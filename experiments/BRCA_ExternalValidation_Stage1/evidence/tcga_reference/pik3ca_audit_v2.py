from pathlib import Path
import csv
from collections import Counter, defaultdict

BASE = Path(r"C:\Users\Gozie\Desktop\BioMedical\TCGA_BRCA_Mutation_Audit")

PROFILED = BASE / "mutation_profiled_patients_v1.csv"
PIK3CA_ROWS = BASE / "pik3ca_mutation_rows_v1.csv"

OUT_LABELS = BASE / "pik3ca_patient_labels_v2.csv"
OUT_SUMMARY = BASE / "pik3ca_label_summary_v2.csv"
OUT_EXCLUDED = BASE / "pik3ca_nonqualifying_rows_v2.csv"


QUALIFYING_CLASSES = {
    "Missense_Mutation",
    "In_Frame_Del",
    "In_Frame_Ins",
    "Nonsense_Mutation",
    "Frame_Shift_Del",
    "Frame_Shift_Ins",
    "Splice_Site",
    "Translation_Start_Site",
    "Nonstop_Mutation",
}


with PROFILED.open("r", encoding="utf-8-sig") as f:
    profiled = list(csv.DictReader(f))

with PIK3CA_ROWS.open("r", encoding="utf-8-sig") as f:
    pik3ca_rows = list(csv.DictReader(f))


rows_by_patient = defaultdict(list)

for r in pik3ca_rows:
    rows_by_patient[r["patient_id"]].append(r)


patient_rows = []
excluded_rows = []

for r in profiled:

    patient_id = r["patient_id"]

    muts = rows_by_patient.get(patient_id, [])

    qualifying = [
        m for m in muts
        if m.get("variant_classification", "").strip()
        in QUALIFYING_CLASSES
    ]

    nonqualifying = [
        m for m in muts
        if m.get("variant_classification", "").strip()
        not in QUALIFYING_CLASSES
    ]

    for m in nonqualifying:
        excluded_rows.append(m)

    label = "MUT" if qualifying else "WT"

    patient_rows.append({
        "patient_id": patient_id,
        "pik3ca_status": label,
        "qualifying_mutation_count": len(qualifying),
        "all_pik3ca_row_count": len(muts),
        "qualifying_classes": " | ".join(sorted({
            m["variant_classification"]
            for m in qualifying
        })),
        "qualifying_hgvsp": " | ".join(sorted({
            m["hgvsp_short"]
            for m in qualifying
            if m.get("hgvsp_short", "").strip()
        })),
    })


status_counts = Counter(
    r["pik3ca_status"]
    for r in patient_rows
)

summary_rows = [
    {
        "metric": "profiled_patients",
        "value": len(patient_rows)
    },
    {
        "metric": "PIK3CA_MUT",
        "value": status_counts.get("MUT", 0)
    },
    {
        "metric": "PIK3CA_WT",
        "value": status_counts.get("WT", 0)
    },
    {
        "metric": "mutation_prevalence_pct",
        "value": round(
            100 * status_counts.get("MUT", 0) / len(patient_rows),
            2
        )
    },
]


with OUT_LABELS.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=list(patient_rows[0].keys())
    )
    writer.writeheader()
    writer.writerows(patient_rows)


with OUT_SUMMARY.open("w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(
        f,
        fieldnames=["metric", "value"]
    )
    writer.writeheader()
    writer.writerows(summary_rows)


if excluded_rows:
    with OUT_EXCLUDED.open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=list(excluded_rows[0].keys())
        )
        writer.writeheader()
        writer.writerows(excluded_rows)


print("=" * 72)
print("TCGA-BRCA PIK3CA PATIENT LABEL AUDIT V2")
print("=" * 72)

for r in summary_rows:
    print(f"{r['metric']:<30} {r['value']}")

print()
print("QUALIFYING VARIANT CLASSES OBSERVED")
print("-" * 72)

observed = Counter()

for r in patient_rows:
    if r["pik3ca_status"] == "MUT":
        for cls in r["qualifying_classes"].split(" | "):
            if cls:
                observed[cls] += 1

for cls, n in observed.most_common():
    print(f"{cls:<35} {n}")

print()
print("Created:")
print(f"  {OUT_LABELS.name}")
print(f"  {OUT_SUMMARY.name}")
print(f"  {OUT_EXCLUDED.name}")