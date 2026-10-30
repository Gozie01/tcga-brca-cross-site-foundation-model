"""Validate and migrate partial v1 checkpoints without changing any vectors."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np

LEGACY_CODE_SHA256 = "be1ab441766c1e9e6fd06108a0655bfd7b73019ebaa5ff5ffd631df62331af78"


def fingerprint(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode()).hexdigest()


def edge_summary(slide, rows, file_name):
    width, height = map(int, slide.dimensions)
    x = rows.x_level0.to_numpy(dtype=np.int64)
    y = rows.y_level0.to_numpy(dtype=np.int64)
    size = rows.read_size_level0.to_numpy(dtype=np.int64)
    if np.any(size <= 0) or np.any(x < 0) or np.any(y < 0):
        raise ValueError(f"Invalid frozen coordinate: {file_name}")
    right, bottom = np.maximum(x + size - width, 0), np.maximum(y + size - height, 0)
    crossing = (right > 0) | (bottom > 0)
    detail = rows.loc[crossing, ["patient_id", "file_name", "patch_id", "x_level0",
                               "y_level0", "read_size_level0"]].copy()
    detail["slide_width"] = width
    detail["slide_height"] = height
    detail["right_overflow_pixels"] = right[crossing]
    detail["bottom_overflow_pixels"] = bottom[crossing]
    summary = {"file_name": file_name, "width": width, "height": height,
               "patches": len(rows), "edge_crossing_patches": int(crossing.sum()),
               "fully_outside_patches": int(((x >= width) | (y >= height)).sum()),
               "max_right_overflow_pixels": int(right.max(initial=0)),
               "max_bottom_overflow_pixels": int(bottom.max(initial=0))}
    return summary, detail


def migrate_partial(out, new_protocol, template, coordinates, paths, open_slide,
                    sha, coordinate_hash, atomic_json):
    """Read-only preflight, then backed-up, idempotent metadata migration.

    Every inherited vector must have been computed entirely inside its slide.
    WSI bytes, vector bytes, coordinate rows and all fitting settings must match.
    Vector files are never written by this function.
    """
    out = Path(out)
    protocol_path = out / "extraction_protocol.json"
    if not protocol_path.exists():
        return
    marker_path = out / "boundary_checkpoint_migration.json"
    current = json.loads(protocol_path.read_text())
    marker = json.loads(marker_path.read_text()) if marker_path.exists() else None
    if current == new_protocol and (marker is None or marker.get("status") == "COMPLETE"):
        return
    if (out / "feature_freeze.json").exists():
        raise ValueError("Migration is limited to incomplete extraction outputs")
    old = dict(new_protocol)
    for key in ["edge_handling", "boundary_fix_code_sha256"]:
        old.pop(key)
    old["version"] = "brca_revision2_hibou_v1"
    old["code_sha256"] = LEGACY_CODE_SHA256
    if old.get("pilot_only"):
        raise ValueError("Pilot checkpoints cannot be migrated into production")
    if current != old and current != new_protocol:
        raise ValueError("Existing model, processor, code, coordinates or extraction settings differ")
    old_hash, new_hash = fingerprint(old), fingerprint(new_protocol)
    if marker is not None and (marker.get("old_protocol_sha256") != old_hash or
                               marker.get("new_protocol_sha256") != new_hash):
        raise ValueError("Migration journal belongs to another protocol")
    checkpoint_dir = out / "slide_checkpoints"
    known = {hashlib.sha256(row.file_name.encode()).hexdigest()[:20]: row
             for row in template.itertuples(index=False)}
    meta_files = sorted(checkpoint_dir.glob("*.json"))
    vector_files = {p.stem for p in checkpoint_dir.glob("*.npy")}
    if vector_files != {p.stem for p in meta_files}:
        raise ValueError("Unpaired checkpoint metadata/vector; no files were changed")
    backup = out / "boundary_v1_backup"
    if (backup / "extraction_protocol.json").exists():
        if json.loads((backup / "extraction_protocol.json").read_text()) != old:
            raise ValueError("Legacy protocol backup mismatch")
    validated = []
    for i, meta_path in enumerate(meta_files):
        entry = known.get(meta_path.stem)
        if entry is None:
            raise ValueError("Checkpoint not in locked slide cohort")
        meta = json.loads(meta_path.read_text())
        array_path = meta_path.with_suffix(".npy")
        rows = coordinates[coordinates.file_name.eq(entry.file_name)]
        if (meta.get("file_name") != entry.file_name or meta.get("patient_id") != entry.patient_id
                or meta.get("patches") != len(rows) or len(rows) != int(entry.patch_count)
                or meta.get("protocol_sha256") not in {old_hash, new_hash}
                or meta.get("coordinates_sha256") != coordinate_hash(rows)
                or meta.get("mean_sha256") != sha(array_path)):
            raise ValueError(f"Checkpoint validation failed: {entry.file_name}")
        vector = np.load(array_path, allow_pickle=False)
        if vector.shape != (768,) or vector.dtype != np.float32 or not np.isfinite(vector).all() or not np.linalg.norm(vector):
            raise ValueError("Invalid inherited slide vector")
        path = Path(paths[entry.file_name])
        print(f"VERIFY LEGACY CHECKPOINT {i+1}/{len(meta_files)} {entry.file_name}", flush=True)
        if sha(path) != meta.get("wsi_sha256"):
            raise ValueError(f"Source slide changed: {entry.file_name}")
        with open_slide(str(path)) as slide:
            edges, _ = edge_summary(slide, rows, entry.file_name)
        if edges["edge_crossing_patches"]:
            raise ValueError("Inherited checkpoint has edge patches; cannot certify v1 equivalence")
        legacy_meta = {**meta, "protocol_sha256": old_hash}
        backup_path = backup / meta_path.name
        if backup_path.exists() and json.loads(backup_path.read_text()) != legacy_meta:
            raise ValueError("Legacy checkpoint backup mismatch")
        validated.append((meta_path, legacy_meta))
        if (i + 1) % 20 == 0 or i + 1 == len(meta_files):
            print(f"CHECKPOINT PREFLIGHT {i+1}/{len(meta_files)} VERIFIED", flush=True)
    journal = {"status": "PREPARED", "old_protocol_sha256": old_hash,
               "new_protocol_sha256": new_hash, "reused_slide_count": len(validated),
               "compatibility": "All reused frozen patches fully inside slide; pixel path unchanged",
               "vector_files_unchanged": True,
               "checkpoints": [{"metadata": path.name, "file_name": meta["file_name"],
                                 "mean_sha256": meta["mean_sha256"],
                                 "wsi_sha256": meta["wsi_sha256"],
                                 "coordinates_sha256": meta["coordinates_sha256"]}
                                for path, meta in validated]}
    if marker is not None and marker.get("checkpoints") != journal["checkpoints"]:
        raise ValueError("Checkpoint set changed during an incomplete migration")
    # All source/vector/coordinate/protocol checks above precede any mutation.
    backup.mkdir(exist_ok=True)
    atomic_json(backup / "extraction_protocol.json", old)
    for path, meta in validated:
        atomic_json(backup / path.name, meta)
    atomic_json(marker_path, journal)
    for path, meta in validated:
        atomic_json(path, {**meta, "protocol_sha256": new_hash})
    atomic_json(protocol_path, new_protocol)
    atomic_json(marker_path, {**journal, "status": "COMPLETE"})
    print(f"MIGRATION COMPLETE: {len(validated)} verified vectors retained unchanged", flush=True)
