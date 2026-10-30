from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
import random
import time
from pathlib import Path

import cv2
import numpy as np
import openslide
from PIL import Image


VERSION = "stage9a_multisite_patch_coordinates_v1"
BASE = Path(r"C:\Users\Gozie\Desktop\BioMedical\TCGA_BRCA_Mutation_Audit")
DATA_ROOT = Path(r"D:\TCGA_BRCA_WSI\targeted_sites_v1")
OUT_ROOT = BASE / "multisite_patch_coordinates_v1"

SITE_CONFIG = {
    "A2": {"root": DATA_ROOT / "A2_remaining", "expected": 11},
    "C8": {"root": DATA_ROOT / "C8", "expected": 46},
    "E2": {"root": DATA_ROOT / "E2", "expected": 71},
    "BH": {"root": DATA_ROOT / "BH", "expected": 126},
}

# Frozen Batch01 extraction policy.
TARGET_MPP = 0.50
PATCH_OUT = 256
LOWRES_TISSUE_MIN = 0.50
LOWRES_GREEN_MAX = 0.05
HIGHRES_TISSUE_MIN = 0.35
WHITE_MAX = 0.60
ARTIFACT_MAX = 0.01
MAX_PATCHES_PER_SLIDE = 1000
RANDOM_SEED = 2026
FOCUS_POLICY = "audit_only"

# Frozen final QC. These are applied after the initial cap is filled; rejected
# patches are not backfilled, matching the locked Batch01 workflow.
BG_TO_RED_MAX = 1.10
MEAN_LOCAL_STD_MIN = 0.50


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", encoding="utf-8") as f:
        json.dump(obj, f, indent=2)
    os.replace(tmp, path)


def atomic_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    os.replace(tmp, path)


def rgb_to_hsv(arr: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    x = arr.astype(np.float32) / 255.0
    mx = x.max(axis=2)
    mn = x.min(axis=2)
    diff = mx - mn
    saturation = np.zeros_like(mx)
    nz = mx > 0
    saturation[nz] = diff[nz] / mx[nz]
    return saturation, mx


def majority_filter(mask: np.ndarray, radius: int) -> np.ndarray:
    h, w = mask.shape
    padded = np.pad(mask.astype(np.uint8), radius, mode="edge")
    score = np.zeros((h, w), dtype=np.uint16)
    window = 2 * radius + 1
    for dy in range(window):
        for dx in range(window):
            score += padded[dy:dy + h, dx:dx + w]
    return score >= ((window * window) // 2)


def build_tissue_mask(arr: np.ndarray) -> np.ndarray:
    saturation, value = rgb_to_hsv(arr)
    stained = (saturation > 0.045) & (value < 0.985)
    pale_tissue = (value < 0.93) & (saturation > 0.015)
    dark_tissue = value < 0.80
    mask = stained | pale_tissue | dark_tissue
    near_white = (
        (arr[..., 0] > 245)
        & (arr[..., 1] > 245)
        & (arr[..., 2] > 245)
    )
    mask &= ~near_white
    mask = majority_filter(mask, radius=2)
    return majority_filter(mask, radius=3)


def build_green_mask(arr: np.ndarray) -> np.ndarray:
    r = arr[..., 0].astype(np.int16)
    g = arr[..., 1].astype(np.int16)
    b = arr[..., 2].astype(np.int16)
    return (g > 70) & (g > r + 20) & (g > b + 15)


def initial_qc_metrics(arr: np.ndarray) -> dict:
    r = arr[..., 0].astype(np.int16)
    g = arr[..., 1].astype(np.int16)
    b = arr[..., 2].astype(np.int16)
    saturation, value = rgb_to_hsv(arr)
    white = (r > 235) & (g > 235) & (b > 235)
    tissue = ((saturation > 0.04) & (value < 0.985)) | (value < 0.80)
    green = (g > 70) & (g > r + 20) & (g > b + 15)
    cyan_blue = (b > 80) & (b > r + 20) & (g > r + 10)
    gray = (
        0.299 * arr[..., 0]
        + 0.587 * arr[..., 1]
        + 0.114 * arr[..., 2]
    ).astype(np.float32)
    center = -4.0 * gray
    up = np.zeros_like(gray)
    down = np.zeros_like(gray)
    left = np.zeros_like(gray)
    right = np.zeros_like(gray)
    up[1:, :] = gray[:-1, :]
    down[:-1, :] = gray[1:, :]
    left[:, 1:] = gray[:, :-1]
    right[:, :-1] = gray[:, 1:]
    focus = float((center + up + down + left + right).var())
    return {
        "tissue_fraction_hr": float(tissue.mean()),
        "white_fraction": float(white.mean()),
        "artifact_fraction": float((green | cyan_blue).mean()),
        "focus_score": focus,
    }


def final_locked_metrics(arr: np.ndarray) -> dict:
    # Exact color formula recovered from the locked Batch01 audit.
    rgb_float = arr.astype(np.float32)
    mean_r = float(np.mean(rgb_float[..., 0]))
    mean_g = float(np.mean(rgb_float[..., 1]))
    mean_b = float(np.mean(rgb_float[..., 2]))
    bg_to_red_ratio = float(((mean_g + mean_b) / 2.0) / max(mean_r, 1e-6))

    # Exact structural formula recovered from the locked Batch01 audit.
    gray = cv2.cvtColor(arr, cv2.COLOR_RGB2GRAY).astype(np.float32)
    local_mean = cv2.blur(gray, (9, 9))
    local_sq_mean = cv2.blur(gray * gray, (9, 9))
    local_var = np.maximum(local_sq_mean - local_mean * local_mean, 0.0)
    mean_local_std = float(np.mean(np.sqrt(local_var)))
    return {
        "mean_r": mean_r,
        "mean_g": mean_g,
        "mean_b": mean_b,
        "bg_to_red_ratio": bg_to_red_ratio,
        "mean_local_std": mean_local_std,
    }


def patient_id(file_name: str) -> str:
    parts = file_name.split(".svs")[0].split("-")
    return "-".join(parts[:3]) if len(parts) >= 3 else ""


def tss(file_name: str) -> str:
    parts = file_name.split("-")
    return parts[1] if len(parts) >= 2 else ""


def slide_seed(file_name: str) -> int:
    digest = hashlib.sha256(f"{RANDOM_SEED}|{file_name}".encode()).hexdigest()
    return int(digest[:16], 16)


FIELDS = [
    "site", "site_slide_index", "patient_id", "tss", "file_name",
    "source_path", "patch_id", "x_level0", "y_level0",
    "read_size_level0", "output_size", "target_mpp", "source_mpp_x",
    "source_mpp_y", "mask_level", "mask_downsample",
    "tissue_fraction_low", "green_fraction_low", "tissue_fraction_hr",
    "white_fraction", "artifact_fraction", "focus_score", "qc_pass",
    "mean_r", "mean_g", "mean_b", "bg_to_red_ratio", "mean_local_std",
    "reject_locked_color", "reject_locked_structural",
    "reject_new_locked_union", "retain_after_new_locked_qc",
]


def process_slide(site: str, index: int, total: int, path: Path, site_out: Path) -> dict:
    token = hashlib.sha256(path.name.encode("utf-8")).hexdigest()[:16]
    checkpoint_dir = site_out / "slide_checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    csv_path = checkpoint_dir / f"{token}.csv"
    meta_path = checkpoint_dir / f"{token}.json"

    if csv_path.exists() and meta_path.exists():
        try:
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if (
                meta.get("status") == "COMPLETE"
                and meta.get("source_size") == path.stat().st_size
                and meta.get("checkpoint_sha256") == sha256_file(csv_path)
            ):
                print(f"[{index:03d}/{total:03d}] RESUME-SKIP {path.name}")
                return meta
        except Exception:
            pass

    print(f"[{index:03d}/{total:03d}] PROCESS {path.name}")
    start = time.time()
    rows: list[dict] = []
    slide = None
    try:
        slide = openslide.OpenSlide(str(path))
        props = slide.properties
        mpp_x_raw = props.get(openslide.PROPERTY_NAME_MPP_X, "")
        mpp_y_raw = props.get(openslide.PROPERTY_NAME_MPP_Y, "")
        if mpp_x_raw == "":
            raise RuntimeError("MPP_X missing from slide metadata")
        mpp_x = float(mpp_x_raw)
        mpp_y = float(mpp_y_raw) if mpp_y_raw != "" else mpp_x
        read_size = round(PATCH_OUT * TARGET_MPP / mpp_x)
        if read_size <= 0:
            raise RuntimeError(f"Invalid read_size={read_size}")

        level = slide.level_count - 1
        level_dimensions = slide.level_dimensions[level]
        ds = float(slide.level_downsamples[level])
        low_arr = np.asarray(
            slide.read_region((0, 0), level, level_dimensions).convert("RGB")
        )
        tissue_mask = build_tissue_mask(low_arr)
        green_mask = build_green_mask(low_arr)
        h, w = tissue_mask.shape
        low_patch = max(1, round(read_size / ds))
        step = low_patch
        candidates = []
        for y in range(0, h - low_patch + 1, step):
            for x in range(0, w - low_patch + 1, step):
                tf = float(tissue_mask[y:y + low_patch, x:x + low_patch].mean())
                gf = float(green_mask[y:y + low_patch, x:x + low_patch].mean())
                if tf >= LOWRES_TISSUE_MIN and gf <= LOWRES_GREEN_MAX:
                    candidates.append((x, y, tf, gf))

        rng = random.Random(slide_seed(path.name))
        rng.shuffle(candidates)
        tested = rejected_initial = 0
        reject_low_tissue = reject_background = reject_artifact = 0

        for lx, ly, tissue_low, green_low in candidates:
            if len(rows) >= MAX_PATCHES_PER_SLIDE:
                break
            tested += 1
            x0, y0 = round(lx * ds), round(ly * ds)
            patch = slide.read_region((x0, y0), 0, (read_size, read_size)).convert("RGB")
            patch = patch.resize((PATCH_OUT, PATCH_OUT), Image.Resampling.LANCZOS)
            arr = np.asarray(patch, dtype=np.uint8)
            initial = initial_qc_metrics(arr)
            low_fail = initial["tissue_fraction_hr"] < HIGHRES_TISSUE_MIN
            bg_fail = initial["white_fraction"] > WHITE_MAX
            artifact_fail = initial["artifact_fraction"] > ARTIFACT_MAX
            if low_fail or bg_fail or artifact_fail:
                rejected_initial += 1
                reject_low_tissue += int(low_fail)
                reject_background += int(bg_fail)
                reject_artifact += int(artifact_fail)
                continue

            final = final_locked_metrics(arr)
            reject_color = final["bg_to_red_ratio"] >= BG_TO_RED_MAX
            reject_structural = final["mean_local_std"] < MEAN_LOCAL_STD_MIN
            reject_union = reject_color or reject_structural
            patch_number = len(rows) + 1
            rows.append({
                "site": site,
                "site_slide_index": index,
                "patient_id": patient_id(path.name),
                "tss": tss(path.name),
                "file_name": path.name,
                "source_path": str(path),
                "patch_id": patch_number,
                "x_level0": x0,
                "y_level0": y0,
                "read_size_level0": read_size,
                "output_size": PATCH_OUT,
                "target_mpp": TARGET_MPP,
                "source_mpp_x": round(mpp_x, 6),
                "source_mpp_y": round(mpp_y, 6),
                "mask_level": level,
                "mask_downsample": round(ds, 6),
                "tissue_fraction_low": round(tissue_low, 6),
                "green_fraction_low": round(green_low, 6),
                "tissue_fraction_hr": round(initial["tissue_fraction_hr"], 6),
                "white_fraction": round(initial["white_fraction"], 6),
                "artifact_fraction": round(initial["artifact_fraction"], 6),
                "focus_score": round(initial["focus_score"], 3),
                "qc_pass": "YES",
                "mean_r": round(final["mean_r"], 6),
                "mean_g": round(final["mean_g"], 6),
                "mean_b": round(final["mean_b"], 6),
                "bg_to_red_ratio": round(final["bg_to_red_ratio"], 9),
                "mean_local_std": round(final["mean_local_std"], 9),
                "reject_locked_color": reject_color,
                "reject_locked_structural": reject_structural,
                "reject_new_locked_union": reject_union,
                "retain_after_new_locked_qc": not reject_union,
            })

        atomic_csv(csv_path, rows, FIELDS)
        retained = sum(str(r["retain_after_new_locked_qc"]).lower() == "true" for r in rows)
        meta = {
            "version": VERSION,
            "status": "COMPLETE",
            "site": site,
            "site_slide_index": index,
            "file_name": path.name,
            "patient_id": patient_id(path.name),
            "source_path": str(path),
            "source_size": path.stat().st_size,
            "coarse_candidates": len(candidates),
            "hr_tested": tested,
            "initial_accepted": len(rows),
            "initial_rejected": rejected_initial,
            "reject_low_tissue": reject_low_tissue,
            "reject_high_background": reject_background,
            "reject_artifact": reject_artifact,
            "final_color_rejects": sum(bool(r["reject_locked_color"]) for r in rows),
            "final_structural_rejects": sum(bool(r["reject_locked_structural"]) for r in rows),
            "final_union_rejects": len(rows) - retained,
            "final_retained": retained,
            "cap_reached": len(rows) == MAX_PATCHES_PER_SLIDE,
            "elapsed_seconds": round(time.time() - start, 2),
            "checkpoint_csv": str(csv_path),
            "checkpoint_sha256": sha256_file(csv_path),
        }
        atomic_json(meta_path, meta)
        print(
            f"  initial={len(rows)} retained={retained} "
            f"final_reject={len(rows)-retained} elapsed={meta['elapsed_seconds']:.1f}s"
        )
        return meta
    except Exception as exc:
        meta = {
            "version": VERSION,
            "status": "ERROR",
            "site": site,
            "site_slide_index": index,
            "file_name": path.name,
            "patient_id": patient_id(path.name),
            "source_path": str(path),
            "source_size": path.stat().st_size,
            "elapsed_seconds": round(time.time() - start, 2),
            "error": str(exc),
        }
        atomic_json(meta_path, meta)
        print(f"  ERROR: {exc}")
        return meta
    finally:
        if slide is not None:
            slide.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Resumable locked-QC multisite patch extraction")
    parser.add_argument("--site", required=True, choices=SITE_CONFIG.keys())
    parser.add_argument("--limit", type=int, default=None, help="Optional diagnostic slide limit")
    args = parser.parse_args()
    site = args.site
    root = SITE_CONFIG[site]["root"]
    expected = SITE_CONFIG[site]["expected"]
    if not root.exists():
        raise FileNotFoundError(f"Site root not found: {root}")
    slides = sorted(root.rglob("*.svs"))
    if len(slides) != expected:
        raise RuntimeError(f"{site}: expected {expected} SVS files, found {len(slides)}")
    selected = slides[: args.limit] if args.limit else slides
    site_out = OUT_ROOT / site
    site_out.mkdir(parents=True, exist_ok=True)
    print("=" * 100)
    print(f"{VERSION} | site={site} | slides={len(selected)}/{len(slides)}")
    print("Locked policy: 0.50 MPP, 256px, max 1000/slide, seed 2026")
    print("Final gates: bg_to_red < 1.10; mean_local_std >= 0.50; no backfill")
    print("=" * 100)
    results = [
        process_slide(site, i, len(selected), path, site_out)
        for i, path in enumerate(selected, start=1)
    ]
    errors = [r for r in results if r.get("status") != "COMPLETE"]
    print("=" * 100)
    print(f"Completed={len(results)-len(errors)} Errors={len(errors)}")
    print(f"Output: {site_out}")
    print("Run the audit/freeze script only after Errors=0.")
    print("=" * 100)
    if errors:
        raise SystemExit(2)


if __name__ == "__main__":
    main()
