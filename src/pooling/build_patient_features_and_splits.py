from __future__ import annotations

import hashlib
import json
import os
import platform
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(os.environ.get("TCGA_BRCA_PIPELINE_ROOT", Path.cwd())).expanduser().resolve()
OUT = ROOT / "stage10a_cohort_split_freeze_v1"
SEED = 20260905
DIM = 1024
HOLDOUT_SITES = ("A2", "C8", "E2", "BH")
EXPECTED = {
    "patches": 345937, "slides": 347, "patients": 345,
    "batch_patches": 92773, "batch_slides": 93, "batch_patients": 91,
    "multisite_patches": 253164, "multisite_slides": 254,
    "multisite_patients": 254,
}
ENDPOINTS = {
    "ER_binary": ("er_evaluable", "er_status", {"Negative": 0, "Positive": 1}),
    "PR_binary": ("pr_evaluable", "pr_status", {"Negative": 0, "Positive": 1}),
    "HER2_binary": ("her2_evaluable", "her2_status", {"Negative": 0, "Positive": 1}),
    "PIK3CA_binary": ("pik3ca_evaluable", "pik3ca_status", {"WT": 0, "MUT": 1}),
    "PAM50_primary_4class": (
        "pam50_primary4_evaluable", "pam50_primary_4class_label",
        {"BRCA_LumA": 0, "BRCA_LumB": 1, "BRCA_Basal": 2, "BRCA_Her2": 3},
    ),
}


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def atomic_json(path: Path, obj: dict) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    tmp.write_text(json.dumps(obj, indent=2), encoding="utf-8")
    os.replace(tmp, path)


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(tmp, index=False, encoding="utf-8-sig")
    os.replace(tmp, path)


def atomic_npy(path: Path, array: np.ndarray) -> None:
    tmp = path.with_suffix(path.suffix + ".tmp")
    with tmp.open("wb") as f:
        np.save(f, array, allow_pickle=False)
    os.replace(tmp, path)


def stable_score(*parts: object) -> str:
    return hashlib.sha256("|".join(map(str, parts)).encode("utf-8")).hexdigest()


def validate_embedding_artifact(embedding_path: Path, index_path: Path,
                                metadata_path: Path, rows: int,
                                slides: int, patients: int):
    for path in (embedding_path, index_path, metadata_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing required artifact: {path}")
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("status") != "LOCKED_COMPLETE":
        raise RuntimeError(f"Embedding artifact is not locked: {metadata_path}")
    matrix = np.load(embedding_path, mmap_mode="r", allow_pickle=False)
    index = pd.read_csv(index_path, low_memory=False)
    checks = {
        "matrix_shape": matrix.shape == (rows, DIM),
        "matrix_dtype": matrix.dtype == np.float32,
        "index_rows": len(index) == rows,
        "slide_count": index["file_name"].nunique() == slides,
        "patient_count": index["patient_id"].nunique() == patients,
        "embedding_hash": sha256_file(embedding_path) == metadata.get("embeddings_sha256"),
        "index_hash": sha256_file(index_path) == metadata.get("index_sha256"),
    }
    failed = [k for k, value in checks.items() if not value]
    if failed:
        raise RuntimeError(f"Embedding validation failed ({metadata_path.name}): {failed}")
    index["embedding_index"] = pd.to_numeric(
        index["embedding_index"], errors="raise"
    ).astype(np.int64)
    index = index.sort_values("embedding_index").reset_index(drop=True)
    if not np.array_equal(index["embedding_index"].to_numpy(), np.arange(rows)):
        raise RuntimeError("Embedding index is not contiguous and aligned")
    return matrix, index, metadata, checks


def aggregate_matrix(matrix, codes: np.ndarray, group_count: int,
                     chunk_size: int = 4096):
    sums = np.zeros((group_count, matrix.shape[1]), dtype=np.float64)
    counts = np.zeros(group_count, dtype=np.int64)
    for start in range(0, matrix.shape[0], chunk_size):
        stop = min(start + chunk_size, matrix.shape[0])
        np.add.at(sums, codes[start:stop], np.asarray(matrix[start:stop], dtype=np.float64))
        np.add.at(counts, codes[start:stop], 1)
    if np.any(counts == 0):
        raise RuntimeError("An aggregation group contains no embeddings")
    return (sums / counts[:, None]).astype(np.float32), counts


def pool_source(matrix, index: pd.DataFrame, source: str):
    required = {"embedding_index", "patient_id", "file_name"}
    missing = sorted(required - set(index.columns))
    if missing:
        raise RuntimeError(f"Missing embedding-index columns: {missing}")
    index = index.copy()
    index["patient_id"] = index["patient_id"].astype(str).str.strip()
    index["file_name"] = index["file_name"].astype(str).str.strip()
    slide_keys = index["patient_id"] + "||" + index["file_name"]
    slide_codes, unique_slides = pd.factorize(slide_keys, sort=True)
    slide_features, slide_patch_counts = aggregate_matrix(
        matrix, slide_codes.astype(np.int64), len(unique_slides)
    )
    slide_patients = np.array([x.split("||", 1)[0] for x in unique_slides])
    slide_names = np.array([x.split("||", 1)[1] for x in unique_slides])
    slide_table = pd.DataFrame({
        "source": source,
        "source_slide_feature_row": np.arange(len(unique_slides)),
        "patient_id": slide_patients,
        "file_name": slide_names,
        "patch_count": slide_patch_counts,
    })
    patient_codes, unique_patients = pd.factorize(slide_table["patient_id"], sort=True)
    patient_features, patient_slide_counts = aggregate_matrix(
        slide_features, patient_codes.astype(np.int64), len(unique_patients)
    )
    patch_counts = (
        slide_table.groupby("patient_id")["patch_count"].sum()
        .reindex(unique_patients).to_numpy(dtype=np.int64)
    )
    patient_table = pd.DataFrame({
        "source": source,
        "source_patient_feature_row": np.arange(len(unique_patients)),
        "patient_id": unique_patients,
        "slide_count": patient_slide_counts,
        "patch_count": patch_counts,
    })
    return patient_features, patient_table, slide_table


def endpoint_table(master: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    eval_col, label_col, mapping = ENDPOINTS[endpoint]
    keep = pd.to_numeric(master[eval_col], errors="coerce").eq(1)
    table = master.loc[keep].copy()
    table["endpoint"] = endpoint
    table["label_name"] = table[label_col].astype(str)
    table["label"] = table["label_name"].map(mapping)
    if table["label"].isna().any():
        bad = sorted(table.loc[table["label"].isna(), "label_name"].unique())
        raise RuntimeError(f"Unrecognized {endpoint} labels: {bad}")
    table["label"] = table["label"].astype(int)
    return table


def assign_internal_folds(table: pd.DataFrame, endpoint: str) -> pd.DataFrame:
    output = table[[
        "patient_id", "feature_row", "source", "tss", "label", "label_name"
    ]].copy()
    output["endpoint"] = endpoint
    output["internal_fold"] = -1
    for _, group in output.groupby("label", sort=True):
        ordered = sorted(
            group.index,
            key=lambda i: stable_score(SEED, endpoint, output.at[i, "patient_id"]),
        )
        for position, index in enumerate(ordered):
            output.at[index, "internal_fold"] = position % 5
    if (output["internal_fold"] < 0).any():
        raise RuntimeError(f"Internal fold assignment failed for {endpoint}")
    return output


def assign_site_holdout(table: pd.DataFrame, endpoint: str,
                        holdout: str) -> pd.DataFrame:
    output = table[[
        "patient_id", "feature_row", "source", "tss", "label", "label_name"
    ]].copy()
    output["endpoint"] = endpoint
    output["holdout_site"] = holdout
    output["partition"] = np.where(output["tss"].eq(holdout), "test", "train")
    development = output[output["partition"].eq("train")]
    for (site, label), group in development.groupby(["tss", "label"], sort=True):
        n = len(group)
        n_validation = 0 if n < 2 else min(n - 1, max(1, int(round(0.20 * n))))
        ordered = sorted(
            group.index,
            key=lambda i: stable_score(
                SEED, endpoint, holdout, site, label, output.at[i, "patient_id"]
            ),
        )
        output.loc[ordered[:n_validation], "partition"] = "validation"
    expected_labels = set(table["label"].unique())
    for partition in ("train", "validation", "test"):
        present = set(output.loc[output["partition"].eq(partition), "label"].unique())
        if present != expected_labels:
            raise RuntimeError(
                f"{endpoint}/{holdout}/{partition} lacks labels: "
                f"expected {expected_labels}, found {present}"
            )
    return output


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    batch_dir = ROOT / "uni_batch01_full_embeddings_v1"
    multi_dir = ROOT / "uni_multisite_embeddings_v1" / "multisite_frozen_v1"
    b_matrix, b_index, b_meta, b_gates = validate_embedding_artifact(
        batch_dir / "batch01_uni_embeddings_v1.npy",
        batch_dir / "batch01_uni_embeddings_index_v1.csv",
        batch_dir / "batch01_uni_embeddings_metadata_v1.json",
        EXPECTED["batch_patches"], EXPECTED["batch_slides"], EXPECTED["batch_patients"],
    )
    m_matrix, m_index, m_meta, m_gates = validate_embedding_artifact(
        multi_dir / "multisite_uni_embeddings_v1.npy",
        multi_dir / "multisite_uni_embeddings_index_v1.csv",
        multi_dir / "multisite_uni_embeddings_freeze_metadata_v1.json",
        EXPECTED["multisite_patches"], EXPECTED["multisite_slides"],
        EXPECTED["multisite_patients"],
    )
    print("Pooling Batch01: patch mean -> slide mean -> patient mean")
    bf, bp, bs = pool_source(b_matrix, b_index, "batch01")
    print("Pooling multisite: patch mean -> slide mean -> patient mean")
    mf, mp, ms = pool_source(m_matrix, m_index, "multisite")
    overlap_patients = sorted(set(bp["patient_id"]) & set(mp["patient_id"]))
    overlap_slides = sorted(set(bs["file_name"]) & set(ms["file_name"]))
    if overlap_patients or overlap_slides:
        raise RuntimeError(
            f"Cross-source overlap: patients={len(overlap_patients)}, slides={len(overlap_slides)}"
        )

    feature_lookup = {}
    membership_frames = []
    for features, patients in ((bf, bp), (mf, mp)):
        for row in patients.itertuples(index=False):
            feature_lookup[row.patient_id] = features[row.source_patient_feature_row]
        membership_frames.append(patients)
    membership = pd.concat(membership_frames, ignore_index=True)
    patient_order = sorted(feature_lookup)
    patient_features = np.stack([feature_lookup[p] for p in patient_order]).astype(np.float32)
    feature_index = pd.DataFrame({
        "feature_row": np.arange(len(patient_order), dtype=np.int64),
        "patient_id": patient_order,
    }).merge(membership, on="patient_id", how="left", validate="one_to_one")
    feature_index["tss"] = feature_index["patient_id"].str.split("-").str[1]

    master_path = (
        ROOT / "fullcohort_site_endpoint_design_v1" /
        "fullcohort_site_endpoint_patient_master_v1.csv"
    )
    design_meta_path = (
        ROOT / "fullcohort_site_endpoint_design_v1" /
        "fullcohort_site_endpoint_design_summary_v1.json"
    )
    master = pd.read_csv(master_path, low_memory=False)
    design_meta = json.loads(design_meta_path.read_text(encoding="utf-8"))
    expected_master_hash = design_meta["outputs"]["patient_master_sha256"]
    if sha256_file(master_path) != expected_master_hash:
        raise RuntimeError("Full-cohort endpoint master hash does not match its audit metadata")
    if master["patient_id"].duplicated().any():
        raise RuntimeError("Duplicate patients in full-cohort endpoint master")
    cohort = feature_index.merge(master, on="patient_id", how="left", validate="one_to_one",
                                 suffixes=("", "_endpoint"))
    if cohort["er_evaluable"].isna().any():
        missing = cohort.loc[cohort["er_evaluable"].isna(), "patient_id"].tolist()
        raise RuntimeError(f"Embedding patients missing from endpoint master: {missing}")
    if not cohort["tss"].eq(cohort["tss_endpoint"]).all():
        raise RuntimeError("TSS mismatch between embeddings and endpoint master")
    cohort.drop(columns=["tss_endpoint"], inplace=True)

    slide_audit = pd.concat([bs, ms], ignore_index=True).sort_values(
        ["source", "patient_id", "file_name"]
    ).reset_index(drop=True)
    internal_frames = []
    site_frames = []
    eligibility_rows = []
    for endpoint in ENDPOINTS:
        eligible = endpoint_table(cohort, endpoint)
        internal = assign_internal_folds(eligible, endpoint)
        internal_frames.append(internal)
        counts = eligible["label_name"].value_counts().sort_index().to_dict()
        eligibility_rows.append({
            "endpoint": endpoint, "evaluable_n": len(eligible),
            "class_count": eligible["label"].nunique(),
            "class_counts_json": json.dumps(counts, sort_keys=True),
        })
        for holdout in HOLDOUT_SITES:
            site_frames.append(assign_site_holdout(eligible, endpoint, holdout))

    internal_splits = pd.concat(internal_frames, ignore_index=True)
    site_splits = pd.concat(site_frames, ignore_index=True)
    eligibility = pd.DataFrame(eligibility_rows)
    split_summary = (
        site_splits.groupby(["endpoint", "holdout_site", "partition", "label_name"])
        .size().rename("patient_n").reset_index()
    )
    internal_summary = (
        internal_splits.groupby(["endpoint", "internal_fold", "label_name"])
        .size().rename("patient_n").reset_index()
    )

    features_out = OUT / "stage10a_uni_patient_features_v1.npy"
    index_out = OUT / "stage10a_uni_patient_features_index_v1.csv"
    cohort_out = OUT / "stage10a_patient_endpoint_master_v1.csv"
    slides_out = OUT / "stage10a_slide_pooling_audit_v1.csv"
    eligibility_out = OUT / "stage10a_endpoint_eligibility_summary_v1.csv"
    internal_out = OUT / "stage10a_internal_5fold_manifest_v1.csv"
    internal_summary_out = OUT / "stage10a_internal_5fold_summary_v1.csv"
    site_out = OUT / "stage10a_site_holdout_manifest_v1.csv"
    site_summary_out = OUT / "stage10a_site_holdout_summary_v1.csv"
    metadata_out = OUT / "stage10a_cohort_split_freeze_metadata_v1.json"
    atomic_npy(features_out, patient_features)
    atomic_csv(index_out, feature_index)
    atomic_csv(cohort_out, cohort)
    atomic_csv(slides_out, slide_audit)
    atomic_csv(eligibility_out, eligibility)
    atomic_csv(internal_out, internal_splits)
    atomic_csv(internal_summary_out, internal_summary)
    atomic_csv(site_out, site_splits)
    atomic_csv(site_summary_out, split_summary)

    expected_endpoint_counts = {
        "ER_binary": 337, "PR_binary": 337, "HER2_binary": 314,
        "PIK3CA_binary": 345, "PAM50_primary_4class": 335,
    }
    gates = {
        "expected_patch_count": len(b_index) + len(m_index) == EXPECTED["patches"],
        "expected_slide_count": len(slide_audit) == EXPECTED["slides"],
        "expected_patient_count": len(cohort) == EXPECTED["patients"],
        "zero_cross_source_patient_overlap": len(overlap_patients) == 0,
        "zero_cross_source_slide_overlap": len(overlap_slides) == 0,
        "all_embedding_patients_endpoint_matched": len(cohort) == EXPECTED["patients"],
        "patient_feature_shape": patient_features.shape == (EXPECTED["patients"], DIM),
        "patient_features_finite": bool(np.isfinite(patient_features).all()),
        "patient_features_nonzero": bool(np.all(np.linalg.norm(patient_features, axis=1) > 0)),
        "patch_count_conservation": int(feature_index["patch_count"].sum()) == EXPECTED["patches"],
        "slide_count_conservation": int(feature_index["slide_count"].sum()) == EXPECTED["slides"],
        "endpoint_counts": dict(zip(eligibility["endpoint"], eligibility["evaluable_n"]))
                           == expected_endpoint_counts,
        "four_site_holdouts_per_endpoint": site_splits.groupby("endpoint")[
            "holdout_site"
        ].nunique().eq(4).all(),
        "no_patient_partition_overlap": site_splits.groupby(
            ["endpoint", "holdout_site", "patient_id"]
        )["partition"].nunique().eq(1).all(),
        "five_internal_folds_per_endpoint": internal_splits.groupby("endpoint")[
            "internal_fold"
        ].nunique().eq(5).all(),
    }
    failed = [k for k, value in gates.items() if not bool(value)]
    if failed:
        raise RuntimeError(f"Stage 10A freeze gates failed: {failed}")
    outputs = [
        features_out, index_out, cohort_out, slides_out, eligibility_out,
        internal_out, internal_summary_out, site_out, site_summary_out,
    ]
    metadata = {
        "pipeline_version": "stage10a_cohort_split_freeze_v1",
        "status": "LOCKED_COMPLETE",
        "seed": SEED,
        "modeling_unit": "TCGA patient",
        "aggregation": "patch mean to slide mean to patient mean",
        "cohort": {
            "patch_embeddings": EXPECTED["patches"], "wsis": EXPECTED["slides"],
            "patients": EXPECTED["patients"], "embedding_dimension": DIM,
            "batch01_patients": EXPECTED["batch_patients"],
            "multisite_patients": EXPECTED["multisite_patients"],
            "cross_source_patient_overlap": 0,
        },
        "endpoint_evaluable_counts": expected_endpoint_counts,
        "holdout_sites": list(HOLDOUT_SITES),
        "split_policies": {
            "internal": "endpoint-specific deterministic stratified 5-fold patient split",
            "site_holdout": (
                "endpoint-specific external test site; deterministic 20% validation "
                "sampling within development-site/class strata; remaining patients train"
            ),
            "leakage_control": "all slides and patches from a patient remain in one partition",
        },
        "software": {
            "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "upstream": {
            "batch_embeddings_sha256": b_meta["embeddings_sha256"],
            "batch_index_sha256": b_meta["index_sha256"],
            "multisite_embeddings_sha256": m_meta["embeddings_sha256"],
            "multisite_index_sha256": m_meta["index_sha256"],
            "endpoint_master_sha256": expected_master_hash,
        },
        "validation": {
            "batch_embedding_gates": b_gates,
            "multisite_embedding_gates": m_gates,
            "freeze_gates": {k: bool(v) for k, v in gates.items()},
        },
        "outputs": {
            path.name: {"path": str(path), "sha256": sha256_file(path)} for path in outputs
        },
    }
    atomic_json(metadata_out, metadata)
    print("=" * 110)
    print("STAGE 10A COHORT AND SPLIT FREEZE")
    print("Status                         : LOCKED_COMPLETE")
    print(f"Patches / WSIs / patients      : {EXPECTED['patches']:,} / {EXPECTED['slides']} / {EXPECTED['patients']}")
    print(f"Patient feature matrix         : {patient_features.shape}")
    print("Cross-source patient overlap   : 0")
    for row in eligibility.itertuples(index=False):
        print(f"{row.endpoint:<31}: {row.evaluable_n:>3} evaluable")
    for name, passed in gates.items():
        print(f"{name:<47}: {'PASS' if passed else 'FAIL'}")
    print(f"Metadata                       : {metadata_out}")
    print("=" * 110)


if __name__ == "__main__":
    main()
