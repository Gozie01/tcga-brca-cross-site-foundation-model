from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


BASE = Path(os.environ.get("TCGA_BRCA_PROJECT_ROOT", Path.cwd())).resolve()
ROOT = BASE / "uni_multisite_embeddings_v1"
SITES = ("A2", "C8", "E2", "BH")
EXPECTED = {
    "A2": (10987, 11, 11), "C8": (45999, 46, 46),
    "E2": (70998, 71, 71), "BH": (125180, 126, 126),
}
TOTAL_ROWS = 253164
TOTAL_SLIDES = 254
TOTAL_PATIENTS = 254
DIM = 1024


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False, encoding="utf-8-sig")
    os.replace(tmp, path)


def main() -> None:
    site_arrays = []
    site_indexes = []
    site_summaries = []
    site_meta = {}
    gates = {}
    for site in SITES:
        site_dir = ROOT / site
        emb_path = site_dir / f"{site}_uni_embeddings_v1.npy"
        idx_path = site_dir / f"{site}_uni_embeddings_index_v1.csv"
        sum_path = site_dir / f"{site}_uni_embeddings_wsi_summary_v1.csv"
        meta_path = site_dir / f"{site}_uni_embeddings_metadata_v1.json"
        for path in (emb_path, idx_path, sum_path, meta_path):
            if not path.exists():
                raise FileNotFoundError(f"Missing Stage 9B artifact: {path}")
        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        arr = np.load(emb_path, mmap_mode="r", allow_pickle=False)
        idx = pd.read_csv(idx_path)
        summary = pd.read_csv(sum_path)
        rows, slides, patients = EXPECTED[site]
        checks = {
            "status_locked": meta.get("status") == "LOCKED_COMPLETE",
            "shape": arr.shape == (rows, DIM),
            "dtype_float32": arr.dtype == np.float32,
            "index_rows": len(idx) == rows,
            "slide_count": idx["file_name"].nunique() == slides == len(summary),
            "patient_count": idx["patient_id"].nunique() == patients,
            "finite": bool(np.isfinite(arr).all()),
            "embedding_hash": sha256_file(emb_path) == meta.get("embeddings_sha256"),
            "index_hash": sha256_file(idx_path) == meta.get("index_sha256"),
            "summary_hash": sha256_file(sum_path) == meta.get("wsi_summary_sha256"),
            "local_indices": bool(np.array_equal(
                idx["embedding_index"].to_numpy(), np.arange(rows)
            )),
            "no_duplicate_coordinates": not idx.duplicated(
                ["file_name", "patch_id", "x_level0", "y_level0"]
            ).any(),
        }
        gates[site] = checks
        failed = [k for k, passed in checks.items() if not passed]
        if failed:
            raise RuntimeError(f"{site} audit failed: {failed}")
        idx.insert(0, "site", site)
        summary.insert(0, "site", site)
        site_arrays.append((site, arr, emb_path))
        site_indexes.append(idx)
        site_summaries.append(summary)
        site_meta[site] = meta

    combined_index = pd.concat(site_indexes, ignore_index=True)
    combined_summary = pd.concat(site_summaries, ignore_index=True)
    combined_index.rename(columns={"embedding_index": "site_embedding_index"}, inplace=True)
    combined_index.insert(0, "embedding_index", np.arange(TOTAL_ROWS, dtype=np.int64))
    global_gates = {
        "expected_rows": len(combined_index) == TOTAL_ROWS,
        "expected_slides": combined_index["file_name"].nunique() == TOTAL_SLIDES,
        "expected_patients": combined_index["patient_id"].nunique() == TOTAL_PATIENTS,
        "no_duplicate_rows": not combined_index.duplicated(
            ["site", "file_name", "patch_id", "x_level0", "y_level0"]
        ).any(),
        "no_duplicate_patch_ids": not combined_index.duplicated(
            ["site", "file_name", "patch_id"]
        ).any(),
        "contiguous_global_indices": bool(np.array_equal(
            combined_index["embedding_index"].to_numpy(), np.arange(TOTAL_ROWS)
        )),
        "expected_summary_rows": len(combined_summary) == TOTAL_SLIDES,
    }
    failed = [k for k, passed in global_gates.items() if not passed]
    if failed:
        raise RuntimeError(f"Global pre-freeze audit failed: {failed}")

    final_dir = ROOT / "multisite_frozen_v1"
    final_dir.mkdir(parents=True, exist_ok=True)
    emb_out = final_dir / "multisite_uni_embeddings_v1.npy"
    idx_out = final_dir / "multisite_uni_embeddings_index_v1.csv"
    summary_out = final_dir / "multisite_uni_embeddings_wsi_summary_v1.csv"
    meta_out = final_dir / "multisite_uni_embeddings_freeze_metadata_v1.json"
    tmp = emb_out.with_suffix(".npy.tmp")
    matrix = np.lib.format.open_memmap(
        tmp, mode="w+", dtype=np.float32, shape=(TOTAL_ROWS, DIM)
    )
    cursor = 0
    for site, arr, _ in site_arrays:
        end = cursor + len(arr)
        matrix[cursor:end] = arr
        cursor = end
    matrix.flush()
    del matrix
    os.replace(tmp, emb_out)
    atomic_csv(idx_out, combined_index)
    atomic_csv(summary_out, combined_summary)

    frozen = np.load(emb_out, mmap_mode="r", allow_pickle=False)
    finite = True
    zero_norm = 0
    for start in range(0, TOTAL_ROWS, 4096):
        block = np.asarray(frozen[start:start + 4096])
        finite = finite and bool(np.isfinite(block).all())
        zero_norm += int(np.sum(np.linalg.norm(block, axis=1) == 0))
    post_gates = {
        "combined_shape": frozen.shape == (TOTAL_ROWS, DIM),
        "combined_dtype_float32": frozen.dtype == np.float32,
        "finite_values": finite,
        "zero_norm_embedding_count": zero_norm == 0,
    }
    failed = [k for k, passed in post_gates.items() if not passed]
    if failed:
        raise RuntimeError(f"Global post-freeze audit failed: {failed}")
    metadata = {
        "pipeline_version": "stage9b_multisite_uni_v1",
        "status": "LOCKED_COMPLETE",
        "site_order": list(SITES),
        "model": "MahmoodLab/UNI",
        "model_identifier": "hf-hub:MahmoodLab/uni",
        "embedding_dimension": DIM,
        "embedding_dtype": "float32",
        "inference_dtype": "float16_autocast",
        "total_embeddings": TOTAL_ROWS,
        "unique_wsi_count": TOTAL_SLIDES,
        "unique_patient_count": TOTAL_PATIENTS,
        "embedding_matrix_shape": list(frozen.shape),
        "site_totals": {
            site: {"embeddings": EXPECTED[site][0], "wsis": EXPECTED[site][1],
                   "patients": EXPECTED[site][2]} for site in SITES
        },
        "embeddings_file": str(emb_out),
        "embeddings_sha256": sha256_file(emb_out),
        "index_file": str(idx_out),
        "index_sha256": sha256_file(idx_out),
        "wsi_summary_file": str(summary_out),
        "wsi_summary_sha256": sha256_file(summary_out),
        "source_site_metadata_sha256": {
            site: sha256_file(ROOT / site / f"{site}_uni_embeddings_metadata_v1.json")
            for site in SITES
        },
        "validation": {
            "site_gates": gates, "global_pre_freeze_gates": global_gates,
            "global_post_freeze_gates": post_gates,
        },
    }
    atomic_json(meta_out, metadata)
    print("=" * 110)
    print("STAGE 9B MULTISITE UNI AUDIT/FREEZE")
    print("Status              : LOCKED_COMPLETE")
    print(f"WSIs / patients     : {TOTAL_SLIDES} / {TOTAL_PATIENTS}")
    print(f"Embeddings          : {TOTAL_ROWS:,}")
    print(f"Matrix shape        : {frozen.shape}")
    for site in SITES:
        print(f"{site:<20}: PASS ({EXPECTED[site][0]:,})")
    for name, passed in {**global_gates, **post_gates}.items():
        print(f"{name:<40}: {'PASS' if passed else 'FAIL'}")
    print(f"Metadata            : {meta_out}")
    print("=" * 110)


if __name__ == "__main__":
    main()
