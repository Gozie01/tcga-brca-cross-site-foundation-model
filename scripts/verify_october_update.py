from pathlib import Path
import ast, csv, hashlib, json
R = Path(__file__).resolve().parents[1]
def digest(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    manifest = json.loads((R / "UPDATE_FILES_SHA256.json").read_text())
    for rel, expected in manifest.items():
        p = R / rel
        if not p.is_file() or digest(p) != expected:
            raise ValueError("Missing or changed update file: " + rel)
        if p.suffix.lower() == ".py":
            ast.parse(p.read_text(encoding="utf-8-sig"), filename=rel)
    uni = R / "experiments/brca_external_frozen_v1"
    audit = json.loads((uni / "summary.json").read_text())
    for rel, expected in audit["inputs_sha256"].items():
        if digest(uni / rel) != expected:
            raise ValueError("UNI audit mismatch: " + rel)
    matched = R / "experiments/brca_matched_external_v1"
    with (matched / "tcga_fit_audit.csv").open(newline="", encoding="utf-8-sig") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 20
    for row in rows:
        p = matched / "models" / (row["endpoint"] + "_" + row["site"] + "_hibou.joblib")
        if digest(p) != row["model_sha256"]:
            raise ValueError("Hibou audit mismatch: " + p.name)
    print("PASS: update checksums, Python syntax, UNI frozen inputs and all 40 classifiers.")
if __name__ == "__main__": main()
