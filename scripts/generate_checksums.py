from __future__ import annotations
import hashlib
from pathlib import Path

root = Path(__file__).resolve().parents[1] / "derived_data"
output = root / "checksums_sha256.txt"
rows = []
for path in sorted(p for p in root.rglob("*") if p.is_file() and p != output):
    digest = hashlib.sha256(path.read_bytes()).hexdigest()
    rows.append(f"{digest}  {path.relative_to(root).as_posix()}")
output.write_text("\n".join(rows) + "\n", encoding="utf-8")
