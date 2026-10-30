"""Exact-coordinate, frozen Hibou-B extraction with slide checkpoints.

The live HF/OpenSlide adapter requires a local WSI inventory and model access.
Pilot output is separate and cannot be used as a complete feature freeze.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

from audit_assets import PACKAGE, sha, verify_inputs
from boundary_migration import migrate_partial, edge_summary

MODEL = "histai/hibou-b"
DIM = 768


def atomic_json(path, payload):
    tmp = Path(str(path) + ".tmp")
    tmp.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    tmp.replace(path)


def atomic_array(path, array):
    tmp = Path(str(path) + ".tmp")
    with tmp.open("wb") as file:
        np.save(file, array, allow_pickle=False)
    tmp.replace(path)


def slide_coordinate_hash(frame):
    keys = ["patient_id", "file_name", "patch_id", "x_level0", "y_level0", "read_size_level0"]
    return hashlib.sha256(frame[keys].to_csv(index=False, lineterminator="\n").encode()).hexdigest()


def stream_slide_mean(slide, rows, encode_batch, batch_size, dimension=DIM):
    """Float32 patch output -> float64 sum -> float32 slide mean, as in UNI."""
    total = np.zeros(dimension, dtype=np.float64)
    count = 0
    for start in range(0, len(rows), batch_size):
        batch = rows.iloc[start:start + batch_size]
        images = []
        for row in batch.itertuples(index=False):
            size = int(row.read_size_level0)
            x, y = int(row.x_level0), int(row.y_level0)
            if size <= 0 or x < 0 or y < 0:
                raise ValueError("Invalid frozen coordinate")
            # Match original Stage9B UNI read_patch exactly: keep the full
            # frozen request and OpenSlide RGBA -> RGB behavior at slide edges.
            region = slide.read_region((x, y), 0, (size, size)).convert("RGB")
            if region.size != (256, 256):
                region = region.resize((256, 256), Image.Resampling.LANCZOS)
            images.append(region)
        values = np.asarray(encode_batch(images), dtype=np.float32)
        if values.shape != (len(images), dimension) or not np.isfinite(values).all():
            raise ValueError("Invalid encoder output")
        if np.any(np.linalg.norm(values, axis=1) == 0):
            raise ValueError("Zero-norm patch features")
        total += values.astype(np.float64).sum(axis=0)
        count += len(images)
    if count != len(rows) or count == 0:
        raise ValueError("Incomplete slide extraction")
    return (total / count).astype(np.float32)


def pool_patients(slide_features, slide_table, patient_index):
    """Mean of float32 slide means; patients never weighted by patch count."""
    result = []
    for patient in patient_index.patient_id:
        positions = np.flatnonzero(slide_table.patient_id.eq(patient))
        if not len(positions):
            raise ValueError(f"Missing patient slide: {patient}")
        result.append(slide_features[positions].astype(np.float64).mean(axis=0).astype(np.float32))
    return np.stack(result)


def load_backend(args):
    import torch
    import openslide
    from huggingface_hub import HfApi, snapshot_download
    from transformers import AutoImageProcessor, AutoModel
    if args.model_snapshot:
        snapshot = args.model_snapshot.resolve()
        revision = "local_snapshot"
    else:
        # Resolves one immutable publisher commit; access must already be granted.
        revision = HfApi().model_info(MODEL, revision=args.model_revision).sha
        snapshot = Path(snapshot_download(MODEL, revision=revision,
            allow_patterns=["*.json", "*.py", "*.safetensors"]))
    files = {p.relative_to(snapshot).as_posix(): sha(p) for p in snapshot.rglob("*")
             if p.is_file() and p.suffix in {".json", ".py", ".safetensors"}}
    if not any(name.endswith(".safetensors") for name in files):
        raise ValueError("A safetensors model snapshot is required")
    processor = AutoImageProcessor.from_pretrained(str(snapshot), trust_remote_code=True, local_files_only=True)
    model = AutoModel.from_pretrained(str(snapshot), trust_remote_code=True,
                                     local_files_only=True, use_safetensors=True)
    device = args.device
    if device == "auto":
        device = "cuda" if torch.cuda.is_available() else "cpu"
    if device == "cuda" and not torch.cuda.is_available():
        raise ValueError("CUDA requested but unavailable")
    model = model.to(device).eval()
    torch.manual_seed(20261019)
    torch.backends.cudnn.benchmark = False
    if device == "cuda":
        torch.backends.cuda.matmul.allow_tf32 = False
        torch.backends.cudnn.allow_tf32 = False
    def encode(images):
        batch = processor(images=images, return_tensors="pt").to(device)
        with torch.inference_mode():
            output = model(**batch).pooler_output
        return output.float().cpu().numpy()
    identity = {"model": MODEL, "dimension": DIM, "revision": revision, "snapshot_file_hashes": files,
        "processor": processor.to_dict(), "device": device, "precision": "float32_no_autocast",
        "batch_size": args.batch_size, "torch": torch.__version__,
        "gpu": torch.cuda.get_device_name(0) if device == "cuda" else None}
    return encode, openslide.OpenSlide, identity


def extract(args):
    lock = verify_inputs()
    inventory_dir = args.audit.resolve()
    audit = json.loads((inventory_dir / "asset_audit.json").read_text())
    if audit["status"] != "READY_FOR_PILOT":
        raise ValueError("Resolve missing/ambiguous slides and scan errors before extraction")
    if audit["coordinate_manifest_sha256"] != lock["files"]["patch_coordinates.csv.gz"]:
        raise ValueError("Audit used a different coordinate manifest")
    inventory = pd.read_csv(inventory_dir / "slide_inventory.csv")
    template = pd.read_csv(PACKAGE / "data/slide_manifest.csv")
    if inventory.file_name.duplicated().any():
        raise ValueError("Duplicate inventory identity")
    if set(inventory.file_name) != set(template.file_name) or not inventory.status.eq("FOUND").all():
        raise ValueError("Inventory does not cover all locked slides")
    coordinates = pd.read_csv(PACKAGE / "data/patch_coordinates.csv.gz")
    encode, open_slide, model_identity = load_backend(args)
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    protocol = {"version": "brca_revision2_hibou_v2_edgefix", "code_sha256": sha(__file__),
        "coordinate_sha256": lock["files"]["patch_coordinates.csv.gz"], "model": model_identity,
        "edge_handling": "Original Stage9B read_region full frozen request -> RGB -> 256 Lanczos; no clipping, shifting or patch exclusion",
        "boundary_fix_code_sha256": sha(Path(__file__).with_name("boundary_migration.py")),
        "pooling": "float32 patch -> float64 sum -> float32 slide mean -> equal slide mean per patient",
        "pilot_only": args.pilot_patches > 0, "pilot_patches": args.pilot_patches}
    if args.migrate_existing:
        migrate_partial(out, protocol, template, coordinates,
            inventory.set_index("file_name").path.to_dict(), open_slide, sha,
            slide_coordinate_hash, atomic_json)
    protocol_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    protocol_path = out / "extraction_protocol.json"
    if protocol_path.exists() and json.loads(protocol_path.read_text()) != protocol:
        raise ValueError("Extraction protocol changed; choose a fresh output directory")
    atomic_json(protocol_path, protocol)
    checkpoint_dir = out / "slide_checkpoints"
    checkpoint_dir.mkdir(exist_ok=True)
    paths = inventory.set_index("file_name").path.to_dict()
    selected = template.iloc[:1] if args.pilot_patches else template
    vectors, slide_records, edge_records, edge_details = [], [], [], []
    for i, entry in enumerate(selected.itertuples(index=False)):
        rows = coordinates[coordinates.file_name.eq(entry.file_name)].copy()
        if args.pilot_patches:
            rows = rows.iloc[:args.pilot_patches]
        print(f"CHECK SLIDE {i+1}/{len(selected)} {entry.file_name}", flush=True)
        token = hashlib.sha256(entry.file_name.encode()).hexdigest()[:20]
        array_path = checkpoint_dir / f"{token}.npy"
        meta_path = checkpoint_dir / f"{token}.json"
        expected_hash = slide_coordinate_hash(rows)
        wsi_path = Path(paths[entry.file_name])
        # Full image hashing is performed only here, immediately before inference.
        wsi_hash = sha(wsi_path)
        with open_slide(str(wsi_path)) as edge_slide:
            edge_record, edge_detail = edge_summary(edge_slide, rows, entry.file_name)
        edge_records.append(edge_record)
        edge_details.append(edge_detail)
        if edge_record["edge_crossing_patches"]:
            print(f"ORIGINAL UNI EDGE HANDLING: {edge_record['edge_crossing_patches']} frozen patches", flush=True)
        if meta_path.exists():
            meta = json.loads(meta_path.read_text())
            if (meta["protocol_sha256"] != protocol_hash or meta["coordinates_sha256"] != expected_hash
                    or meta["wsi_sha256"] != wsi_hash or meta["mean_sha256"] != sha(array_path)):
                raise ValueError(f"Changed slide checkpoint: {entry.file_name}")
            mean = np.load(array_path, allow_pickle=False)
        else:
            started = time.perf_counter()
            with open_slide(str(wsi_path)) as slide:
                mean = stream_slide_mean(slide, rows, encode, args.batch_size)
            atomic_array(array_path, mean)
            meta = {"file_name": entry.file_name, "patient_id": entry.patient_id,
                "patches": len(rows), "coordinates_sha256": expected_hash, "wsi_sha256": wsi_hash,
                "mean_sha256": sha(array_path), "protocol_sha256": protocol_hash,
                "elapsed_seconds": time.perf_counter() - started}
            atomic_json(meta_path, meta)
        if mean.shape != (DIM,) or not np.isfinite(mean).all() or np.linalg.norm(mean) == 0:
            raise ValueError("Invalid slide mean")
        vectors.append(mean)
        slide_records.append({"patient_id": entry.patient_id, "file_name": entry.file_name,
                              "patches": len(rows), "wsi_sha256": wsi_hash})
        pd.DataFrame(edge_records).to_csv(out / "slide_edge_audit.csv", index=False)
        pd.concat(edge_details, ignore_index=True).to_csv(out / "edge_patch_coordinates.csv", index=False)
        print(f"SLIDE {i+1}/{len(selected)} COMPLETE {entry.file_name}", flush=True)
    if args.pilot_patches:
        atomic_json(out / "pilot_complete.json", {"status": "PILOT_ONLY", "patches": len(rows),
            "protocol_sha256": protocol_hash, "no_complete_patient_feature_matrix": True})
        print("PILOT_ONLY: no complete feature freeze was created.")
        return
    slide_table = pd.DataFrame(slide_records)
    if slide_table.patches.sum() != 345937 or len(slide_table) != 347:
        raise ValueError("Production patch/slide counts changed")
    patient_index = pd.read_csv(PACKAGE / "data/stage10a_uni_patient_features_index_v1.csv")
    patient = pool_patients(np.stack(vectors), slide_table, patient_index)
    if patient.shape != (345, DIM) or not np.isfinite(patient).all():
        raise ValueError("Patient feature gate failed")
    feature_path = out / "hibou_patient_features.npy"
    index_path = out / "hibou_patient_index.csv"
    atomic_array(feature_path, patient)
    patient_index.to_csv(index_path, index=False)
    slide_table.to_csv(out / "slide_pooling_audit.csv", index=False)
    atomic_json(out / "feature_freeze.json", {"status": "LOCKED_COMPLETE", "model": MODEL,
        "patients": 345, "slides": 347, "patches": 345937, "dimension": DIM,
        "features_sha256": sha(feature_path), "index_sha256": sha(index_path),
        "protocol_sha256": protocol_hash, "coordinate_sha256": protocol["coordinate_sha256"],
        "extraction_protocol_sha256": sha(protocol_path),
        "edge_crossing_patches": sum(r["edge_crossing_patches"] for r in edge_records),
        "edge_crossing_slides": sum(r["edge_crossing_patches"] > 0 for r in edge_records),
        "edge_audit_sha256": sha(out / "slide_edge_audit.csv"),
        "edge_coordinates_sha256": sha(out / "edge_patch_coordinates.csv"),
        "migration_sha256": sha(out / "boundary_checkpoint_migration.json") if (out / "boundary_checkpoint_migration.json").exists() else None})
    print("LOCKED_COMPLETE: 345 x 768 patient features; original patch selection preserved.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--audit", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--batch-size", type=int, default=32)
    parser.add_argument("--device", choices=["auto", "cuda", "cpu"], default="auto")
    parser.add_argument("--model-revision", default="main")
    parser.add_argument("--model-snapshot", type=Path)
    parser.add_argument("--pilot-patches", type=int, default=0)
    parser.add_argument("--migrate-existing", action="store_true")
    args = parser.parse_args()
    if args.batch_size < 1 or args.pilot_patches < 0:
        parser.error("Invalid batch or pilot size")
    extract(args)
