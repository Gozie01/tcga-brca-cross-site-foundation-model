from __future__ import annotations

import argparse
import gc
import hashlib
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image
import openslide
import timm
import torch
from timm.data import create_transform, resolve_model_data_config


BASE = Path(os.environ.get("TCGA_BRCA_PIPELINE_ROOT", Path.cwd())).expanduser().resolve()
WSI_BASE = Path(os.environ.get("TCGA_BRCA_WSI_ROOT", "data/tcga_brca_wsi")).expanduser().resolve()
SITES = ("A2", "C8", "E2", "BH")
EXPECTED = {
    "A2": {"rows": 10987, "slides": 11, "patients": 11},
    "C8": {"rows": 45999, "slides": 46, "patients": 46},
    "E2": {"rows": 70998, "slides": 71, "patients": 71},
    "BH": {"rows": 125180, "slides": 126, "patients": 126},
}
MODEL_ID = "hf-hub:MahmoodLab/uni"
EMBED_DIM = 1024
DEFAULT_BATCH_SIZE = 32
REQUIRED_COLUMNS = (
    "patient_id", "file_name", "patch_id", "x_level0", "y_level0",
    "read_size_level0",
)


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def canonical_slide_hash(frame: pd.DataFrame) -> str:
    payload = frame[list(REQUIRED_COLUMNS)].to_csv(
        index=False, lineterminator="\n"
    ).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False, encoding="utf-8-sig")
    os.replace(tmp, path)


def atomic_npy(path: Path, array: np.ndarray) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as f:
        np.save(f, array, allow_pickle=False)
    os.replace(tmp, path)


def index_svs(root: Path) -> dict[str, Path]:
    if not root.exists():
        raise FileNotFoundError(f"WSI directory does not exist: {root}")
    result: dict[str, Path] = {}
    for path in root.rglob("*.svs"):
        if path.name in result:
            raise RuntimeError(f"Duplicate SVS filename under {root}: {path.name}")
        result[path.name] = path
    return result


def read_patch(slide: openslide.OpenSlide, row: pd.Series) -> Image.Image:
    size = int(row["read_size_level0"])
    patch = slide.read_region(
        (int(row["x_level0"]), int(row["y_level0"])), 0, (size, size)
    ).convert("RGB")
    if patch.size != (256, 256):
        patch = patch.resize((256, 256), Image.Resampling.LANCZOS)
    return patch


def checkpoint_paths(checkpoint_dir: Path, number: int, file_name: str):
    token = hashlib.sha256(file_name.encode("utf-8")).hexdigest()[:12]
    stem = f"{number:03d}_{token}"
    return (
        checkpoint_dir / f"{stem}_embeddings.npy",
        checkpoint_dir / f"{stem}_index.csv",
        checkpoint_dir / f"{stem}_metadata.json",
    )


def verified_checkpoint(paths, rows: pd.DataFrame, file_name: str) -> bool:
    npy_path, index_path, meta_path = paths
    if not all(p.exists() for p in paths):
        return False
    try:
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        arr = np.load(npy_path, mmap_mode="r", allow_pickle=False)
        idx = pd.read_csv(index_path)
        expected_shape = (len(rows), EMBED_DIM)
        checks = [
            meta.get("status") == "CHECKPOINT_COMPLETE",
            meta.get("file_name") == file_name,
            meta.get("slide_manifest_sha256") == canonical_slide_hash(rows),
            meta.get("embedding_shape") == list(expected_shape),
            tuple(arr.shape) == expected_shape,
            arr.dtype == np.float32,
            len(idx) == len(rows),
            sha256_file(npy_path) == meta.get("embeddings_sha256"),
            sha256_file(index_path) == meta.get("index_sha256"),
            np.isfinite(arr).all(),
            list(idx["patch_id"].astype(int)) == list(rows["patch_id"].astype(int)),
        ]
        return all(checks)
    except Exception:
        return False


def extract_slide(model, transform, slide_path: Path, rows: pd.DataFrame,
                  batch_size: int) -> np.ndarray:
    output = np.empty((len(rows), EMBED_DIM), dtype=np.float32)
    slide = openslide.OpenSlide(str(slide_path))
    try:
        for start in range(0, len(rows), batch_size):
            end = min(start + batch_size, len(rows))
            tensors = [transform(read_patch(slide, row))
                       for _, row in rows.iloc[start:end].iterrows()]
            batch = torch.stack(tensors).to("cuda", non_blocking=True)
            with torch.inference_mode(), torch.autocast(
                device_type="cuda", dtype=torch.float16
            ):
                embeddings = model(batch)
            if embeddings.ndim != 2 or embeddings.shape[1] != EMBED_DIM:
                raise RuntimeError(f"Unexpected UNI output: {tuple(embeddings.shape)}")
            output[start:end] = embeddings.float().cpu().numpy()
            del tensors, batch, embeddings
        torch.cuda.synchronize()
    finally:
        slide.close()
    if not np.isfinite(output).all() or np.any(np.linalg.norm(output, axis=1) == 0):
        raise RuntimeError(f"Invalid embedding values for {slide_path.name}")
    return output


def run(site: str, batch_size: int) -> None:
    if not torch.cuda.is_available():
        raise RuntimeError("CUDA is unavailable; production run aborted.")
    expected = EXPECTED[site]
    site_manifest = (
        BASE / "multisite_patch_coordinates_v1" / site /
        f"{site}_final_retained_patch_manifest_v1.csv"
    )
    freeze_meta_path = (
        BASE / "multisite_patch_coordinates_v1" / site /
        f"{site}_patch_manifest_freeze_metadata_v1.json"
    )
    if not site_manifest.exists() or not freeze_meta_path.exists():
        raise FileNotFoundError(f"Missing frozen Stage 9A inputs for {site}")
    freeze_meta = json.loads(freeze_meta_path.read_text(encoding="utf-8"))
    if freeze_meta.get("status") != "LOCKED":
        raise RuntimeError(f"Stage 9A manifest for {site} is not LOCKED")
    manifest = pd.read_csv(site_manifest)
    missing_cols = sorted(set(REQUIRED_COLUMNS) - set(manifest.columns))
    if missing_cols:
        raise RuntimeError(f"Manifest columns missing: {missing_cols}")
    actual = {
        "rows": len(manifest),
        "slides": manifest["file_name"].nunique(),
        "patients": manifest["patient_id"].nunique(),
    }
    if actual != expected:
        raise RuntimeError(f"Expected {expected} for {site}; found {actual}")
    if manifest.duplicated(
        ["file_name", "patch_id", "x_level0", "y_level0"]
    ).any():
        raise RuntimeError("Duplicate patch coordinates in frozen manifest")

    outdir = BASE / "uni_multisite_embeddings_v1" / site
    checkpoint_dir = outdir / "slide_checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    svs_index = index_svs(WSI_BASE / ("A2_remaining" if site == "A2" else site))
    missing = sorted(set(manifest["file_name"].astype(str)) - set(svs_index))
    if missing:
        raise RuntimeError(f"{len(missing)} manifest WSIs missing for {site}: {missing[:5]}")

    wsi_names = manifest["file_name"].drop_duplicates().astype(str).tolist()
    pending = []
    for number, file_name in enumerate(wsi_names, 1):
        rows = manifest[manifest["file_name"] == file_name].reset_index(drop=True)
        paths = checkpoint_paths(checkpoint_dir, number, file_name)
        if not verified_checkpoint(paths, rows, file_name):
            pending.append((number, file_name, rows, paths))

    print("=" * 110)
    print(f"STAGE 9B UNI | site={site} | rows={len(manifest):,} | slides={len(wsi_names)}")
    print(f"Verified checkpoints={len(wsi_names)-len(pending)} | pending={len(pending)}")
    print(f"GPU={torch.cuda.get_device_name(0)} | batch={batch_size}")
    print("=" * 110)

    model = None
    transform = None
    parameter_count = None
    data_config = None
    if pending:
        model = timm.create_model(
            MODEL_ID, pretrained=True, init_values=1e-5, dynamic_img_size=True
        ).to("cuda").eval()
        parameter_count = sum(p.numel() for p in model.parameters())
        data_config = resolve_model_data_config(model)
        transform = create_transform(**data_config, is_training=False)

    run_start = time.perf_counter()
    for number, file_name, rows, paths in pending:
        started = time.perf_counter()
        print(f"[{number:03d}/{len(wsi_names):03d}] PROCESS {file_name} ({len(rows):,})")
        array = extract_slide(model, transform, svs_index[file_name], rows, batch_size)
        index_df = rows[list(REQUIRED_COLUMNS)].copy()
        index_df.insert(0, "slide_embedding_index", np.arange(len(rows), dtype=np.int64))
        npy_path, index_path, meta_path = paths
        atomic_npy(npy_path, array)
        atomic_csv(index_path, index_df)
        meta = {
            "pipeline_version": "stage9b_multisite_uni_v1",
            "status": "CHECKPOINT_COMPLETE",
            "site": site,
            "file_name": file_name,
            "patient_id": str(rows.iloc[0]["patient_id"]),
            "slide_manifest_sha256": canonical_slide_hash(rows),
            "embedding_shape": list(array.shape),
            "embedding_dtype": "float32",
            "model_identifier": MODEL_ID,
            "inference_dtype": "float16_autocast",
            "batch_size": batch_size,
            "embeddings_sha256": sha256_file(npy_path),
            "index_sha256": sha256_file(index_path),
            "elapsed_seconds": time.perf_counter() - started,
        }
        atomic_json(meta_path, meta)
        if not verified_checkpoint(paths, rows, file_name):
            raise RuntimeError(f"Checkpoint verification failed: {file_name}")
        print(f"  CHECKPOINT PASS | {meta['elapsed_seconds']:.1f}s")
        del array, index_df
        gc.collect()
        torch.cuda.empty_cache()

    if model is not None:
        del model
        gc.collect()
        torch.cuda.empty_cache()

    # Reverify every checkpoint before creating the immutable site artifact.
    records = []
    checkpoint_sets = []
    for number, file_name in enumerate(wsi_names, 1):
        rows = manifest[manifest["file_name"] == file_name].reset_index(drop=True)
        paths = checkpoint_paths(checkpoint_dir, number, file_name)
        if not verified_checkpoint(paths, rows, file_name):
            raise RuntimeError(f"Final checkpoint gate failed: {file_name}")
        checkpoint_sets.append((file_name, rows, paths))

    embeddings_out = outdir / f"{site}_uni_embeddings_v1.npy"
    index_out = outdir / f"{site}_uni_embeddings_index_v1.csv"
    summary_out = outdir / f"{site}_uni_embeddings_wsi_summary_v1.csv"
    metadata_out = outdir / f"{site}_uni_embeddings_metadata_v1.json"
    tmp_matrix = embeddings_out.with_suffix(".npy.tmp")
    matrix = np.lib.format.open_memmap(
        tmp_matrix, mode="w+", dtype=np.float32, shape=(len(manifest), EMBED_DIM)
    )
    index_frames = []
    cursor = 0
    for file_name, rows, paths in checkpoint_sets:
        arr = np.load(paths[0], mmap_mode="r", allow_pickle=False)
        idx = pd.read_csv(paths[1])
        start = cursor
        cursor += len(arr)
        matrix[start:cursor] = arr
        idx.insert(0, "embedding_index", np.arange(start, cursor, dtype=np.int64))
        index_frames.append(idx)
        records.append({
            "patient_id": str(rows.iloc[0]["patient_id"]),
            "file_name": file_name, "patch_count": len(rows),
            "embedding_start_index": start,
            "embedding_end_index_exclusive": cursor,
        })
    matrix.flush()
    del matrix
    os.replace(tmp_matrix, embeddings_out)
    index_df = pd.concat(index_frames, ignore_index=True)
    summary_df = pd.DataFrame(records)
    atomic_csv(index_out, index_df)
    atomic_csv(summary_out, summary_df)
    site_array = np.load(embeddings_out, mmap_mode="r", allow_pickle=False)
    if site_array.shape != (expected["rows"], EMBED_DIM) or not np.isfinite(site_array).all():
        raise RuntimeError("Final site embedding matrix validation failed")
    metadata = {
        "pipeline_version": "stage9b_multisite_uni_v1",
        "status": "LOCKED_COMPLETE",
        "site": site,
        "model": "MahmoodLab/UNI",
        "model_identifier": MODEL_ID,
        "parameter_count": parameter_count or 303350784,
        "embedding_dimension": EMBED_DIM,
        "embedding_dtype": "float32",
        "inference_dtype": "float16_autocast",
        "batch_size": batch_size,
        "gpu": torch.cuda.get_device_name(0),
        "torch_version": torch.__version__,
        "timm_version": timm.__version__,
        "frozen_manifest": str(site_manifest),
        "frozen_manifest_sha256": sha256_file(site_manifest),
        "stage9a_freeze_metadata_sha256": sha256_file(freeze_meta_path),
        "total_embeddings": len(index_df),
        "unique_wsi_count": index_df["file_name"].nunique(),
        "unique_patient_count": index_df["patient_id"].nunique(),
        "embedding_matrix_shape": list(site_array.shape),
        "preprocessing": data_config or {
            "input_size": [3, 224, 224], "interpolation": "bilinear",
            "mean": [0.485, 0.456, 0.406], "std": [0.229, 0.224, 0.225],
            "crop_pct": 1, "crop_mode": "center",
        },
        "patch_extraction": {
            "source": "SVS whole-slide image", "coordinate_space": "level-0",
            "intermediate_output_size_px": 256, "target_mpp": 0.5,
            "permanent_patch_jpeg_corpus": False,
        },
        "embeddings_file": str(embeddings_out),
        "embeddings_sha256": sha256_file(embeddings_out),
        "index_file": str(index_out),
        "index_sha256": sha256_file(index_out),
        "wsi_summary_file": str(summary_out),
        "wsi_summary_sha256": sha256_file(summary_out),
        "run_elapsed_seconds": time.perf_counter() - run_start,
        "validation": {
            "all_slide_checkpoints_verified": True,
            "finite_values": True,
            "duplicate_coordinate_count": int(index_df.duplicated(
                ["file_name", "patch_id", "x_level0", "y_level0"]
            ).sum()),
            "duplicate_embedding_index_count": int(index_df["embedding_index"].duplicated().sum()),
            "contiguous_embedding_indices": bool(np.array_equal(
                index_df["embedding_index"].to_numpy(), np.arange(len(index_df))
            )),
        },
    }
    if any(v not in (0, True) for v in metadata["validation"].values()):
        raise RuntimeError(f"Site validation failed: {metadata['validation']}")
    atomic_json(metadata_out, metadata)
    print("=" * 110)
    print(f"SITE {site} LOCKED_COMPLETE | embeddings={len(index_df):,} | shape={site_array.shape}")
    print(f"Output: {outdir}")
    print("=" * 110)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--site", required=True, choices=SITES)
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE)
    args = parser.parse_args()
    if args.batch_size < 1:
        raise ValueError("--batch-size must be positive")
    run(args.site, args.batch_size)


if __name__ == "__main__":
    main()
