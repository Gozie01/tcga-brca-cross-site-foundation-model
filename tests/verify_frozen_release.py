"""Verify release metadata checksums without requiring raw TCGA data."""
from __future__ import annotations
import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

for metadata in (ROOT / "derived_data/audit_metadata").glob("*.json"):
    record = json.loads(metadata.read_text(encoding="utf-8"))
    assert record.get("status") == "LOCKED_COMPLETE", metadata
print("Frozen release metadata: PASS")
