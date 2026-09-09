from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd


VERSION = "stage9a_multisite_patch_coordinates_v1"
BASE = Path(os.environ.get("TCGA_BRCA_PIPELINE_ROOT", Path.cwd())).expanduser().resolve()
OUT_ROOT = BASE / "multisite_patch_coordinates_v1"
EXPECTED = {"A2": 11, "C8": 46, "E2": 71, "BH": 126}
TARGET_MPP = 0.50
PATCH_OUT = 256
MAX_PATCHES = 1000


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def normalize_bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    out = series.astype(str).str.strip().str.lower().map(
        {"true": True, "false": False, "1": True, "0": False}
    )
    if out.isna().any():
        raise RuntimeError("Boolean normalization failed")
    return out.astype(bool)


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def main() -> None:
    parser = argparse.ArgumentParser(description="Audit and freeze one Stage 9A site")
    parser.add_argument("--site", required=True, choices=EXPECTED.keys())
    args = parser.parse_args()
    site = args.site
    site_dir = OUT_ROOT / site
    checkpoint_dir = site_dir / "slide_checkpoints"
    meta_paths = sorted(checkpoint_dir.glob("*.json"))
    metas = [json.loads(p.read_text(encoding="utf-8")) for p in meta_paths]
    complete = [m for m in metas if m.get("status") == "COMPLETE"]
    errors = [m for m in metas if m.get("status") != "COMPLETE"]
    csv_paths = [Path(m["checkpoint_csv"]) for m in complete]

    if len(complete) != EXPECTED[site] or errors:
        raise RuntimeError(
            f"{site}: expected {EXPECTED[site]} complete checkpoints; "
            f"found complete={len(complete)}, errors={len(errors)}"
        )
    for m, p in zip(complete, csv_paths):
        if not p.exists() or sha256_file(p) != m["checkpoint_sha256"]:
            raise RuntimeError(f"Checkpoint integrity failure: {p}")

    frames = [pd.read_csv(p) for p in csv_paths]
    combined = pd.concat(frames, ignore_index=True)
    bool_cols = [
        "reject_locked_color", "reject_locked_structural",
        "reject_new_locked_union", "retain_after_new_locked_qc",
    ]
    for col in bool_cols:
        combined[col] = normalize_bool(combined[col])

    duplicate_rows = int(combined.duplicated().sum())
    duplicate_coords = int(combined.duplicated(["file_name", "x_level0", "y_level0"]).sum())
    duplicate_patch_ids = int(combined.duplicated(["file_name", "patch_id"]).sum())
    malformed_patients = int(
        (~combined["patient_id"].fillna("").str.match(r"^TCGA-[A-Z0-9]{2}-[A-Z0-9]{4}$")).sum()
    )
    union_expected = combined["reject_locked_color"] | combined["reject_locked_structural"]
    union_mismatch = int((combined["reject_new_locked_union"] != union_expected).sum())
    retain_mismatch = int((combined["retain_after_new_locked_qc"] != ~union_expected).sum())
    threshold_mismatch = int((
        (combined["reject_locked_color"] != (combined["bg_to_red_ratio"] >= 1.10))
        | (combined["reject_locked_structural"] != (combined["mean_local_std"] < 0.50))
    ).sum())
    initial_gate_violations = int((
        (combined["tissue_fraction_hr"] < 0.35)
        | (combined["white_fraction"] > 0.60)
        | (combined["artifact_fraction"] > 0.01)
        | (combined["tissue_fraction_low"] < 0.50)
        | (combined["green_fraction_low"] > 0.05)
    ).sum())
    bad_geometry = int((
        (combined["x_level0"] < 0)
        | (combined["y_level0"] < 0)
        | (combined["output_size"] != PATCH_OUT)
        | (~np.isclose(combined["target_mpp"], TARGET_MPP, atol=1e-6))
    ).sum())

    retained = combined[combined["retain_after_new_locked_qc"]].copy()
    rejected = combined[combined["reject_new_locked_union"]].copy()
    per_slide = combined.groupby(["site", "patient_id", "file_name"], as_index=False).agg(
        initial_coordinates=("patch_id", "count"),
        color_rejects=("reject_locked_color", "sum"),
        structural_rejects=("reject_locked_structural", "sum"),
        union_rejects=("reject_new_locked_union", "sum"),
        final_retained_coordinates=("retain_after_new_locked_qc", "sum"),
        min_mean_local_std=("mean_local_std", "min"),
        max_bg_to_red_ratio=("bg_to_red_ratio", "max"),
    )
    zero_patch_slides = EXPECTED[site] - combined["file_name"].nunique()
    gate = all([
        len(errors) == 0,
        combined["file_name"].nunique() == EXPECTED[site],
        combined["patient_id"].nunique() == EXPECTED[site],
        duplicate_rows == 0,
        duplicate_coords == 0,
        duplicate_patch_ids == 0,
        malformed_patients == 0,
        union_mismatch == 0,
        retain_mismatch == 0,
        threshold_mismatch == 0,
        initial_gate_violations == 0,
        bad_geometry == 0,
        zero_patch_slides == 0,
    ])

    combined_path = site_dir / f"{site}_combined_locked_qc_manifest_v1.csv"
    retained_path = site_dir / f"{site}_final_retained_patch_manifest_v1.csv"
    rejected_path = site_dir / f"{site}_final_rejected_patch_manifest_v1.csv"
    summary_path = site_dir / f"{site}_patch_wsi_summary_v1.csv"
    combined.to_csv(combined_path, index=False, encoding="utf-8-sig")
    retained.to_csv(retained_path, index=False, encoding="utf-8-sig")
    rejected.to_csv(rejected_path, index=False, encoding="utf-8-sig")
    per_slide.to_csv(summary_path, index=False, encoding="utf-8-sig")

    metadata = {
        "version": VERSION,
        "status": "LOCKED" if gate else "REVIEW_REQUIRED",
        "site": site,
        "expected_slides": EXPECTED[site],
        "unique_wsi_count": int(combined["file_name"].nunique()),
        "unique_patient_count": int(combined["patient_id"].nunique()),
        "initial_coordinate_count": int(len(combined)),
        "final_retained_count": int(len(retained)),
        "final_rejected_count": int(len(rejected)),
        "retention_fraction": float(len(retained) / len(combined)),
        "locked_qc": {
            "target_mpp": 0.50,
            "output_size_px": 256,
            "lowres_tissue_min": 0.50,
            "lowres_green_max": 0.05,
            "highres_tissue_min": 0.35,
            "white_fraction_max": 0.60,
            "artifact_fraction_max": 0.01,
            "bg_to_red_ratio_reject_if_ge": 1.10,
            "mean_local_std_reject_if_lt": 0.50,
            "focus_policy": "audit_only",
            "max_initial_coordinates_per_slide": 1000,
            "random_seed": 2026,
            "final_gate_backfill": False,
        },
        "checks": {
            "no_checkpoint_errors": len(errors) == 0,
            "expected_slides": combined["file_name"].nunique() == EXPECTED[site],
            "one_patient_per_target_slide": combined["patient_id"].nunique() == EXPECTED[site],
            "no_duplicate_rows": duplicate_rows == 0,
            "no_duplicate_coordinates": duplicate_coords == 0,
            "no_duplicate_patch_ids": duplicate_patch_ids == 0,
            "patient_ids_valid": malformed_patients == 0,
            "union_decisions_consistent": union_mismatch == 0,
            "retention_decisions_consistent": retain_mismatch == 0,
            "final_thresholds_consistent": threshold_mismatch == 0,
            "initial_locked_gates_valid": initial_gate_violations == 0,
            "geometry_valid": bad_geometry == 0,
            "no_zero_patch_slides": zero_patch_slides == 0,
        },
        "outputs": {
            "combined_manifest": str(combined_path),
            "combined_manifest_sha256": sha256_file(combined_path),
            "retained_manifest": str(retained_path),
            "retained_manifest_sha256": sha256_file(retained_path),
            "rejected_manifest": str(rejected_path),
            "rejected_manifest_sha256": sha256_file(rejected_path),
            "wsi_summary": str(summary_path),
            "wsi_summary_sha256": sha256_file(summary_path),
        },
    }
    meta_path = site_dir / f"{site}_patch_manifest_freeze_metadata_v1.json"
    atomic_json(meta_path, metadata)

    print("=" * 100)
    print(f"STAGE 9A SITE AUDIT | {site}")
    print(f"Status              : {metadata['status']}")
    print(f"WSIs / patients     : {metadata['unique_wsi_count']} / {metadata['unique_patient_count']}")
    print(f"Initial coordinates : {len(combined):,}")
    print(f"Final retained      : {len(retained):,}")
    print(f"Final rejected      : {len(rejected):,}")
    print(f"Retention           : {100*len(retained)/len(combined):.4f}%")
    for name, passed in metadata["checks"].items():
        print(f"{name:36s}: {'PASS' if passed else 'FAIL'}")
    print(f"Metadata            : {meta_path}")
    print("=" * 100)
    if not gate:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
