"""Fetch public CPTAC clinical metadata, without downloading WSIs or deriving labels."""
from __future__ import annotations
import argparse
import hashlib
import json
import re
import urllib.parse
import urllib.request
from collections import Counter
from pathlib import Path

STUDY = "brca_cptac_2020"
BASE = "https://www.cbioportal.org/api"


def fetch(path, params=None):
    url = BASE + path + ("?" + urllib.parse.urlencode(params) if params else "")
    request = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": "BRCA-Research-Metadata-Audit/1.0"})
    with urllib.request.urlopen(request, timeout=45) as response:
        raw = response.read()
    return json.loads(raw), {"url": url, "sha256": hashlib.sha256(raw).hexdigest()}


def paginate(path, params=None):
    collected, provenance = [], []
    for page in range(100):
        rows, prov = fetch(path, {**(params or {}), "pageSize": 500, "pageNumber": page, "projection": "DETAILED"})
        if not isinstance(rows, list): raise ValueError("Expected a metadata list")
        collected.extend(rows); provenance.append(prov)
        if len(rows) < 500: return collected, provenance
    raise ValueError("Unexpected metadata pagination length")


def candidate_attributes(attributes):
    rules = {"ER_binary": r"(^|[^A-Z])ER([^A-Z]|$)|ESTROGEN",
             "PR_binary": r"(^|[^A-Z])PR([^A-Z]|$)|PROGESTERONE",
             "HER2_binary": r"HER.?2|ERBB2", "PAM50_primary_4class": r"PAM.?50|INTRINSIC",
             "PIK3CA_binary": r"PIK3CA"}
    result = {}
    for endpoint, pattern in rules.items():
        result[endpoint] = [a for a in attributes if re.search(pattern,
            str(a.get("clinicalAttributeId", "")) + " " + str(a.get("displayName", "")), re.I)]
    return result


def audit(out):
    out = Path(out).resolve(); out.mkdir(parents=True, exist_ok=False)
    study, study_source = fetch(f"/studies/{STUDY}")
    attributes, attr_source = paginate(f"/studies/{STUDY}/clinical-attributes")
    samples, sample_source = paginate(f"/studies/{STUDY}/samples")
    patient_data, patient_source = paginate(f"/studies/{STUDY}/clinical-data", {"clinicalDataType": "PATIENT"})
    sample_data, clinical_source = paginate(f"/studies/{STUDY}/clinical-data", {"clinicalDataType": "SAMPLE"})
    profiles, profile_source = paginate(f"/studies/{STUDY}/molecular-profiles")
    for name, value in [("study", study), ("clinical_attributes", attributes), ("samples", samples),
                        ("patient_clinical_data", patient_data), ("sample_clinical_data", sample_data),
                        ("molecular_profiles", profiles)]:
        (out / f"{name}.json").write_text(json.dumps(value, indent=2), encoding="utf-8")
    candidates = candidate_attributes(attributes)
    distributions = {}
    for endpoint, attrs in candidates.items():
        distributions[endpoint] = {}
        for attribute in attrs:
            key = attribute["clinicalAttributeId"]
            values = [r.get("value") for r in patient_data + sample_data if r.get("clinicalAttributeId") == key]
            distributions[endpoint][key] = dict(Counter(values))
    report = {"status": "METADATA_ONLY_NOT_VALIDATION", "study": STUDY,
        "patients_in_samples": len({r["patientId"] for r in samples}), "samples": len(samples),
        "endpoint_candidates": candidates, "candidate_value_counts": distributions,
        "provenance": [study_source] + attr_source + sample_source + patient_source + clinical_source + profile_source,
        "required_before_validation": [
            "Match actual TCIA slide patient identifiers to molecular records with a documented one-to-one mapping",
            "Verify this is non-TCGA CPTAC, not the older TCGA-derived CPTAC breast resource",
            "Freeze exact ER/PR/HER2 assay-status fields; do not substitute proteomic or PAM50 proxy labels",
            "Verify mutation-profiled denominator before assigning any PIK3CA wild-type labels",
            "Freeze supplied molecular PAM50 classes and exclusions; do not infer PAM50 from receptor status",
            "Select and quality-control diagnostic tumor slides without examining model outcomes",
            "Freeze TCGA-trained models/calibrators before external scores; no external-outcome recalibration",
            "Audit encoder pretraining overlap where documentation permits; unknown overlap must remain a limitation"],
        "no_WSI_download": True, "no_endpoint_labels_assigned": True, "no_external_metrics_computed": True}
    (out / "external_metadata_audit.json").write_text(json.dumps(report, indent=2), encoding="utf-8")
    print(json.dumps({"status": report["status"], "patients": report["patients_in_samples"],
                      "report": str(out / "external_metadata_audit.json")}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()
    audit(args.out)
