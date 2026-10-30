"""Offline, standard-library inventory for the exact locked image cohort."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
import os
import platform
import subprocess
import sys
import time
import zipfile
from pathlib import Path

PACKAGE = Path(__file__).resolve().parent.parent


def sha(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_csv(path):
    with Path(path).open(encoding="utf-8-sig", newline="") as f:
        return list(csv.DictReader(f))


def write_csv(path, rows, columns):
    with Path(path).open("w", encoding="utf-8-sig", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)


def verify_inputs(package=PACKAGE):
    data = package / "data"
    lock = json.loads((data / "input_lock.json").read_text())
    for name, expected in lock["files"].items():
        if sha(data / name) != expected:
            raise ValueError(f"Bundled input changed: {name}")
    return lock


def index_slides(roots, wanted):
    """Inspect filenames only; never hash whole WSI files during this audit."""
    found = {name: [] for name in wanted}
    errors, visited = [], set()
    for root in roots:
        root = Path(root).resolve()
        if not root.is_dir():
            errors.append({"path": str(root), "error": "directory_missing"})
            continue
        def onerror(error):
            errors.append({"path": str(error.filename), "error": type(error).__name__})
        for directory, folders, files in os.walk(root, followlinks=False, onerror=onerror):
            folders[:] = [f for f in folders if f not in {
                ".git", "node_modules", ".venv", ".venv-brca-revision", "__pycache__",
            } and not f.startswith(".venv-")]
            for name in files:
                if name not in found:
                    continue
                path = Path(directory, name).resolve()
                if path in visited:
                    continue
                visited.add(path)
                found[name].append(str(path))
    return found, errors


def hardware_info():
    packages = {}
    for name in ["numpy", "pandas", "scikit-learn", "torch", "torchvision",
                 "transformers", "huggingface-hub", "openslide-python", "openslide-bin", "Pillow"]:
        try:
            packages[name] = importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError:
            packages[name] = None
    gpu = {"probe": "not_available"}
    if packages["torch"]:
        command = "import torch,json; print(json.dumps({'cuda_available':torch.cuda.is_available(),'cuda_runtime':torch.version.cuda,'gpu_name':torch.cuda.get_device_name(0) if torch.cuda.is_available() else None}))"
        try:
            result = subprocess.run([sys.executable, "-c", command], capture_output=True,
                                    text=True, timeout=30, check=True)
            gpu = json.loads(result.stdout.strip().splitlines()[-1])
        except Exception as error:
            gpu = {"probe": "failed", "error_type": type(error).__name__}
    return {"python": platform.python_version(), "platform": platform.system(),
            "packages": packages, "gpu": gpu}


def audit(roots, out, package=PACKAGE, probe_hardware=True):
    lock = verify_inputs(package)
    slides = read_csv(package / "data/slide_manifest.csv")
    if len(slides) != lock["slides"]:
        raise ValueError("Slide manifest count changed")
    found, errors = index_slides(roots, {r["file_name"] for r in slides})
    rows = []
    for slide in slides:
        paths = found[slide["file_name"]]
        rows.append({**slide, "status": "FOUND" if len(paths) == 1 else "MISSING" if not paths else "AMBIGUOUS",
                     "path": paths[0] if len(paths) == 1 else "", "candidate_paths": " | ".join(paths)})
    out = Path(out).resolve()
    out.mkdir(parents=True, exist_ok=False)
    columns = ["source", "patient_id", "file_name", "patch_count", "status", "path", "candidate_paths"]
    write_csv(out / "slide_inventory.csv", rows, columns)
    write_csv(out / "missing_or_ambiguous_slides.csv", [r for r in rows if r["status"] != "FOUND"], columns)
    counts = {state: sum(r["status"] == state for r in rows) for state in ["FOUND", "MISSING", "AMBIGUOUS"]}
    result = {"status": "READY_FOR_PILOT" if counts["FOUND"] == len(slides) and not errors else "ASSET_CHECK_REQUIRED",
              "counts": counts, "patients": lock["patients"], "patches": lock["patches"],
              "searched_roots": [str(Path(r).resolve()) for r in roots], "search_errors": errors,
              "coordinate_manifest_sha256": lock["files"]["patch_coordinates.csv.gz"],
              "hardware": hardware_info() if probe_hardware else {"probe": "skipped_for_test"},
              "note": "Offline inventory only. No weights, images, packages, or data were downloaded."}
    (out / "asset_audit.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    with zipfile.ZipFile(out.with_suffix(".zip"), "w", zipfile.ZIP_DEFLATED) as archive:
        for path in sorted(out.iterdir()):
            archive.write(path, path.name)
    print(json.dumps({"status": result["status"], "counts": counts,
                      "report": str(out / "asset_audit.json"), "return_zip": str(out.with_suffix('.zip'))}, indent=2))
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", action="append", type=Path, required=True)
    parser.add_argument("--out", type=Path, default=Path.cwd() / ("brca_revision2_asset_audit_" + time.strftime("%Y%m%d_%H%M%S")))
    args = parser.parse_args()
    audit(args.root, args.out)
