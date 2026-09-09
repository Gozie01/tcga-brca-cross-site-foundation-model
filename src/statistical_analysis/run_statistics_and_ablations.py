from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import os
import platform
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler


VERSION = "stage11_statistical_analysis_ablations_v1"
SEED = 20260905
C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
ENDPOINTS = (
    "ER_binary", "PR_binary", "HER2_binary", "PIK3CA_binary",
    "PAM50_primary_4class",
)
SITES = ("A2", "C8", "E2", "BH")
ABLATIONS = ("l2_normalized", "pca_32", "pca_64", "pca_128")
METRICS = ("auroc", "auprc", "brier", "log_loss", "ece_10bin")


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def atomic_csv(path: Path, frame: pd.DataFrame) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    frame.to_csv(temporary, index=False, encoding="utf-8-sig")
    os.replace(temporary, path)


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    os.replace(temporary, path)


def load_stage10b_common(root: Path):
    path = root / "Stage10B_Internal_Prediction_Calibration_v1" / \
        "stage10b_run_internal_prediction_calibration_v1.py"
    if not path.is_file():
        raise FileNotFoundError(f"Missing Stage 10B support script: {path}")
    spec = importlib.util.spec_from_file_location("stage10b_common", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def verify_frozen_output(folder: Path, metadata_name: str,
                         required: tuple[str, ...]) -> dict:
    metadata_path = folder / metadata_name
    if not metadata_path.is_file():
        raise FileNotFoundError(metadata_path)
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("status") != "LOCKED_COMPLETE":
        raise RuntimeError(f"Input is not LOCKED_COMPLETE: {metadata_path}")
    outputs = metadata.get("outputs", {})
    for name in required:
        path = folder / name
        if not path.is_file():
            raise FileNotFoundError(path)
        recorded = outputs.get(name, {}).get("sha256")
        if recorded and sha256_file(path) != recorded:
            raise RuntimeError(f"Frozen-input hash mismatch: {path}")
    return metadata


def probability_matrix(frame: pd.DataFrame, prefix: str,
                       classes: np.ndarray) -> np.ndarray:
    return np.column_stack([
        frame[f"{prefix}{label}"].to_numpy(float) for label in classes
    ])


def paired_bootstrap(first_y: np.ndarray, first_p: np.ndarray,
                     second_p: np.ndarray, classes: np.ndarray, reps: int,
                     seed: int, common) -> list[dict]:
    first_point = common.calculate_metrics(first_y, first_p, classes)
    second_point = common.calculate_metrics(first_y, second_p, classes)
    samples = {metric: [] for metric in METRICS}
    rng = np.random.default_rng(seed)
    for _ in range(reps):
        indices = rng.integers(0, len(first_y), len(first_y))
        if len(np.unique(first_y[indices])) != len(classes):
            continue
        try:
            a = common.calculate_metrics(
                first_y[indices], first_p[indices], classes
            )
            b = common.calculate_metrics(
                first_y[indices], second_p[indices], classes
            )
        except ValueError:
            continue
        for metric in METRICS:
            delta = b[metric] - a[metric]
            if np.isfinite(delta):
                samples[metric].append(float(delta))
    rows = []
    for metric in METRICS:
        values = np.asarray(samples[metric], dtype=float)
        if values.size == 0:
            raise RuntimeError(f"No valid bootstrap replicates for {metric}")
        p_value = min(1.0, 2.0 * min(
            (np.sum(values <= 0) + 1) / (len(values) + 1),
            (np.sum(values >= 0) + 1) / (len(values) + 1),
        ))
        rows.append({
            "metric": metric,
            "first_estimate": first_point[metric],
            "second_estimate": second_point[metric],
            "delta_second_minus_first": second_point[metric] - first_point[metric],
            "delta_ci_lower_95": float(np.percentile(values, 2.5)),
            "delta_ci_upper_95": float(np.percentile(values, 97.5)),
            "paired_bootstrap_p_two_sided": float(p_value),
            "bootstrap_valid_replicates": int(len(values)),
        })
    return rows


def estimator(ablation: str, seed: int) -> Pipeline:
    classifier = LogisticRegression(
        class_weight="balanced", solver="lbfgs", max_iter=5000,
        random_state=seed,
    )
    if ablation == "l2_normalized":
        return Pipeline([
            ("normalize", Normalizer(norm="l2")),
            ("classifier", classifier),
        ])
    components = int(ablation.split("_")[1])
    return Pipeline([
        ("scale", StandardScaler()),
        ("pca", PCA(n_components=components, random_state=seed)),
        ("classifier", classifier),
    ])


def checkpoint_valid(paths: dict[str, Path], fingerprint: str) -> bool:
    if not all(path.is_file() for path in paths.values()):
        return False
    try:
        metadata = json.loads(paths["metadata"].read_text(encoding="utf-8"))
        return (
            metadata.get("fingerprint") == fingerprint
            and sha256_file(paths["predictions"]) == metadata["predictions_sha256"]
            and sha256_file(paths["metrics"]) == metadata["metrics_sha256"]
            and sha256_file(paths["model"]) == metadata["model_sha256"]
        )
    except Exception:
        return False


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--bootstrap-reps", type=int, default=1000)
    parser.add_argument("--n-jobs", type=int, default=-1)
    args = parser.parse_args()
    root = args.root
    common = load_stage10b_common(root)

    a_dir = root / "stage10a_cohort_split_freeze_v1"
    b_dir = root / "stage10b_internal_prediction_calibration_v1"
    c_dir = root / "stage10c_site_holdout_domain_shift_v1"
    feature_name = "stage10a_uni_patient_features_v1.npy"
    split_name = "stage10a_internal_5fold_manifest_v1.csv"
    b_prediction_name = "stage10b_internal_oof_predictions_v1.csv"
    c_prediction_name = "stage10c_site_holdout_predictions_v1.csv"
    c_metric_name = "stage10c_site_holdout_metrics_v1.csv"
    a_meta = verify_frozen_output(
        a_dir, "stage10a_cohort_split_freeze_metadata_v1.json",
        (feature_name, split_name),
    )
    b_meta = verify_frozen_output(
        b_dir, "stage10b_internal_prediction_calibration_metadata_v1.json",
        (b_prediction_name, "stage10b_internal_fold_metrics_v1.csv"),
    )
    c_meta = verify_frozen_output(
        c_dir, "stage10c_site_holdout_metadata_v1.json",
        (c_prediction_name, c_metric_name, "stage10c_domain_gap_v1.csv"),
    )

    X = np.load(a_dir / feature_name, mmap_mode="r", allow_pickle=False)
    splits = pd.read_csv(a_dir / split_name, low_memory=False)
    internal = pd.read_csv(b_dir / b_prediction_name, low_memory=False)
    external = pd.read_csv(c_dir / c_prediction_name, low_memory=False)
    external_metrics = pd.read_csv(c_dir / c_metric_name, low_memory=False)
    out = root / VERSION
    checkpoints = out / "ablation_checkpoints"
    checkpoints.mkdir(parents=True, exist_ok=True)

    print("=" * 110)
    print("STAGE 11: PAIRED STATISTICS, DOMAIN SHIFT, AND REPRESENTATION ABLATIONS")
    print("=" * 110)

    calibration_rows = []
    for endpoint in ENDPOINTS:
        table = internal[internal.endpoint.eq(endpoint)].copy()
        y = table.label.astype(int).to_numpy()
        classes = np.sort(np.unique(y))
        raw = probability_matrix(table, "prob_raw_class_", classes)
        calibrated = probability_matrix(table, "prob_calibrated_class_", classes)
        rows = paired_bootstrap(
            y, raw, calibrated, classes, args.bootstrap_reps,
            SEED + sum(map(ord, endpoint + "calibration")), common,
        )
        for row in rows:
            calibration_rows.append({"endpoint": endpoint, **row})
        print(f"Calibration paired statistics: {endpoint} PASS")
    calibration_stats = pd.DataFrame(calibration_rows)

    domain_rows = []
    for endpoint in ENDPOINTS:
        internal_endpoint = internal[internal.endpoint.eq(endpoint)].copy()
        for site in SITES:
            external_site = external[
                external.endpoint.eq(endpoint) & external.holdout_site.eq(site)
            ].copy()
            paired = external_site.merge(
                internal_endpoint, on=["endpoint", "patient_id", "label"],
                suffixes=("_external", "_internal"), validate="one_to_one",
            )
            y = paired.label.astype(int).to_numpy()
            classes = np.sort(np.unique(y))
            for state, prefix in (
                ("raw", "prob_raw_class_"),
                ("calibrated", "prob_calibrated_class_"),
            ):
                external_p = np.column_stack([
                    paired[f"{prefix}{label}_external"].to_numpy(float)
                    for label in classes
                ])
                internal_p = np.column_stack([
                    paired[f"{prefix}{label}_internal"].to_numpy(float)
                    for label in classes
                ])
                rows = paired_bootstrap(
                    y, internal_p, external_p, classes, args.bootstrap_reps,
                    SEED + sum(map(ord, endpoint + site + state)), common,
                )
                for row in rows:
                    domain_rows.append({
                        "endpoint": endpoint, "holdout_site": site,
                        "calibration_state": state, "paired_patients": len(paired),
                        **row,
                    })
            print(f"Paired internal/external statistics: {endpoint} {site} PASS")
    domain_stats = pd.DataFrame(domain_rows)

    heterogeneity = external_metrics.groupby(
        ["endpoint", "calibration"], as_index=False
    ).agg(
        site_mean_auroc=("auroc", "mean"),
        site_sd_auroc=("auroc", "std"),
        worst_site_auroc=("auroc", "min"),
        best_site_auroc=("auroc", "max"),
        site_mean_auprc=("auprc", "mean"),
        worst_site_auprc=("auprc", "min"),
        site_mean_ece=("ece_10bin", "mean"),
        worst_site_ece=("ece_10bin", "max"),
    )
    heterogeneity["auroc_site_range"] = (
        heterogeneity.best_site_auroc - heterogeneity.worst_site_auroc
    )

    ablation_predictions = []
    ablation_metrics = []
    feature_hash = sha256_file(a_dir / feature_name)
    split_hash = sha256_file(a_dir / split_name)
    for endpoint in ENDPOINTS:
        endpoint_split = splits[splits.endpoint.eq(endpoint)].copy()
        classes = np.sort(endpoint_split.label.unique().astype(int))
        scoring = "average_precision" if len(classes) == 2 else "neg_log_loss"
        for fold in range(5):
            train = endpoint_split[endpoint_split.internal_fold.ne(fold)]
            test = endpoint_split[endpoint_split.internal_fold.eq(fold)]
            train_rows = train.feature_row.astype(int).to_numpy()
            test_rows = test.feature_row.astype(int).to_numpy()
            y_train = train.label.astype(int).to_numpy()
            y_test = test.label.astype(int).to_numpy()
            cv_n = min(4, int(pd.Series(y_train).value_counts().min()))
            if cv_n < 2:
                raise RuntimeError(f"Insufficient inner-CV support: {endpoint}/{fold}")
            for ablation in ABLATIONS:
                stem = f"{endpoint}_fold{fold}_{ablation}"
                paths = {
                    "predictions": checkpoints / f"{stem}_predictions.csv",
                    "metrics": checkpoints / f"{stem}_metrics.csv",
                    "model": checkpoints / f"{stem}_model.joblib",
                    "metadata": checkpoints / f"{stem}_metadata.json",
                }
                fingerprint = hashlib.sha256(
                    f"{VERSION}|{feature_hash}|{split_hash}|{endpoint}|{fold}|{ablation}".encode()
                ).hexdigest()
                if checkpoint_valid(paths, fingerprint):
                    prediction = pd.read_csv(paths["predictions"])
                    metric = pd.read_csv(paths["metrics"])
                    status = "VERIFIED SKIP"
                else:
                    started = time.perf_counter()
                    search = GridSearchCV(
                        estimator(ablation, SEED + fold),
                        {"classifier__C": list(C_GRID)}, scoring=scoring,
                        cv=StratifiedKFold(
                            cv_n, shuffle=True, random_state=SEED + fold
                        ),
                        n_jobs=args.n_jobs, error_score="raise",
                    )
                    search.fit(X[train_rows], y_train)
                    probabilities = search.best_estimator_.predict_proba(X[test_rows])
                    values = common.calculate_metrics(y_test, probabilities, classes)
                    prediction = test[[
                        "endpoint", "patient_id", "feature_row", "source", "tss",
                        "label", "label_name", "internal_fold",
                    ]].copy().reset_index(drop=True)
                    prediction["ablation"] = ablation
                    for position, label in enumerate(classes):
                        prediction[f"prob_class_{label}"] = probabilities[:, position]
                    metric = pd.DataFrame([{
                        "endpoint": endpoint, "fold": fold,
                        "ablation": ablation, "train_n": len(train),
                        "test_n": len(test), "inner_cv_splits": cv_n,
                        "best_C": search.best_params_["classifier__C"],
                        "elapsed_seconds": time.perf_counter() - started, **values,
                    }])
                    atomic_csv(paths["predictions"], prediction)
                    atomic_csv(paths["metrics"], metric)
                    temporary_model = paths["model"].with_suffix(".tmp")
                    joblib.dump(search.best_estimator_, temporary_model)
                    os.replace(temporary_model, paths["model"])
                    atomic_json(paths["metadata"], {
                        "status": "CHECKPOINT_COMPLETE", "fingerprint": fingerprint,
                        "predictions_sha256": sha256_file(paths["predictions"]),
                        "metrics_sha256": sha256_file(paths["metrics"]),
                        "model_sha256": sha256_file(paths["model"]),
                    })
                    status = "CHECKPOINT PASS"
                ablation_predictions.append(prediction)
                ablation_metrics.append(metric)
                print(
                    f"{endpoint:<24} fold={fold} {ablation:<16} "
                    f"AUROC={metric.iloc[0].auroc:.4f} | {status}"
                )

    ablation_prediction_table = pd.concat(ablation_predictions, ignore_index=True)
    ablation_metric_table = pd.concat(ablation_metrics, ignore_index=True)
    ablation_summary = ablation_metric_table.groupby(
        ["endpoint", "ablation"], as_index=False
    ).agg(
        mean_auroc=("auroc", "mean"), sd_auroc=("auroc", "std"),
        mean_auprc=("auprc", "mean"), sd_auprc=("auprc", "std"),
        mean_brier=("brier", "mean"), mean_ece=("ece_10bin", "mean"),
    )

    ablation_comparison_rows = []
    for endpoint in ENDPOINTS:
        baseline = internal[internal.endpoint.eq(endpoint)].copy()
        y = baseline.label.astype(int).to_numpy()
        classes = np.sort(np.unique(y))
        baseline_p = probability_matrix(baseline, "prob_raw_class_", classes)
        baseline_key = baseline[["patient_id", "label"]].reset_index(drop=True)
        for ablation in ABLATIONS:
            alternative = ablation_prediction_table[
                ablation_prediction_table.endpoint.eq(endpoint)
                & ablation_prediction_table.ablation.eq(ablation)
            ].copy().reset_index(drop=True)
            joined = baseline_key.merge(
                alternative, on=["patient_id", "label"], validate="one_to_one"
            )
            if len(joined) != len(baseline):
                raise RuntimeError(f"Ablation pairing incomplete: {endpoint}/{ablation}")
            order = baseline_key.merge(
                alternative[["patient_id", "label"]].reset_index(
                    names="alternative_row"
                ),
                on=["patient_id", "label"], validate="one_to_one",
            ).alternative_row.to_numpy(int)
            alternative_p = np.column_stack([
                alternative.iloc[order][f"prob_class_{label}"].to_numpy(float)
                for label in classes
            ])
            rows = paired_bootstrap(
                y, baseline_p, alternative_p, classes, args.bootstrap_reps,
                SEED + sum(map(ord, endpoint + ablation)), common,
            )
            for row in rows:
                ablation_comparison_rows.append({
                    "endpoint": endpoint, "reference": "full_1024_stage10b_raw",
                    "ablation": ablation, **row,
                })
    ablation_comparisons = pd.DataFrame(ablation_comparison_rows)

    outputs = {
        "calibration": out / "stage11_paired_calibration_statistics_v1.csv",
        "domain": out / "stage11_paired_domain_shift_statistics_v1.csv",
        "heterogeneity": out / "stage11_site_heterogeneity_summary_v1.csv",
        "ablation_predictions": out / "stage11_ablation_oof_predictions_v1.csv",
        "ablation_metrics": out / "stage11_ablation_fold_metrics_v1.csv",
        "ablation_summary": out / "stage11_ablation_summary_v1.csv",
        "ablation_comparisons": out / "stage11_paired_ablation_statistics_v1.csv",
    }
    frames = {
        "calibration": calibration_stats, "domain": domain_stats,
        "heterogeneity": heterogeneity,
        "ablation_predictions": ablation_prediction_table,
        "ablation_metrics": ablation_metric_table,
        "ablation_summary": ablation_summary,
        "ablation_comparisons": ablation_comparisons,
    }
    for name, path in outputs.items():
        atomic_csv(path, frames[name])

    expected_ablation_predictions = len(splits) * len(ABLATIONS)
    gates = {
        "frozen_stage10_inputs": all(
            item.get("status") == "LOCKED_COMPLETE"
            for item in (a_meta, b_meta, c_meta)
        ),
        "calibration_statistics_complete": len(calibration_stats) == 25,
        "domain_statistics_complete": len(domain_stats) == 200,
        "twenty_site_models_summarized": heterogeneity.groupby(
            "endpoint"
        ).size().eq(2).all(),
        "one_hundred_ablation_checkpoints": len(list(
            checkpoints.glob("*_metadata.json")
        )) == 100,
        "ablation_predictions_complete": (
            len(ablation_prediction_table) == expected_ablation_predictions
        ),
        "no_ablation_prediction_duplicates": not ablation_prediction_table.duplicated(
            ["endpoint", "ablation", "patient_id"]
        ).any(),
        "paired_ablation_statistics_complete": len(ablation_comparisons) == 100,
        "all_primary_results_finite": (
            np.isfinite(calibration_stats[[
                "first_estimate", "second_estimate", "delta_second_minus_first",
                "delta_ci_lower_95", "delta_ci_upper_95",
                "paired_bootstrap_p_two_sided",
            ]].to_numpy()).all()
            and np.isfinite(domain_stats[[
                "first_estimate", "second_estimate", "delta_second_minus_first",
                "delta_ci_lower_95", "delta_ci_upper_95",
                "paired_bootstrap_p_two_sided",
            ]].to_numpy()).all()
            and np.isfinite(heterogeneity.select_dtypes(
                include=[np.number]
            ).to_numpy()).all()
            and np.isfinite(ablation_metric_table[list(METRICS)].to_numpy()).all()
            and np.isfinite(ablation_summary.select_dtypes(
                include=[np.number]
            ).to_numpy()).all()
            and np.isfinite(ablation_comparisons[[
                "first_estimate", "second_estimate", "delta_second_minus_first",
                "delta_ci_lower_95", "delta_ci_upper_95",
                "paired_bootstrap_p_two_sided",
            ]].to_numpy()).all()
        ),
    }
    if not all(bool(value) for value in gates.values()):
        raise RuntimeError(f"Stage 11 audit failed: {gates}")

    metadata_path = out / "stage11_statistical_analysis_ablations_metadata_v1.json"
    atomic_json(metadata_path, {
        "pipeline_version": VERSION, "status": "LOCKED_COMPLETE",
        "seed": SEED, "bootstrap_reps": args.bootstrap_reps,
        "paired_statistics": (
            "patient-level nonparametric paired bootstrap with percentile 95% CIs "
            "and finite-sample corrected two-sided bootstrap p-values"
        ),
        "ablations": {
            "reference": "Stage 10B full 1024-dimensional raw OOF predictions",
            "variants": list(ABLATIONS),
            "leakage_control": (
                "normalization, scaling, PCA, and model selection fitted only "
                "within each outer-training fold"
            ),
        },
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__, "pandas": pd.__version__,
            "scikit_learn": sklearn.__version__,
        },
        "validation": {name: bool(value) for name, value in gates.items()},
        "outputs": {
            path.name: {"path": str(path), "sha256": sha256_file(path)}
            for path in outputs.values()
        },
    })
    print("=" * 110)
    print("STAGE 11 STATISTICAL ANALYSIS / ABLATION AUDIT")
    print("Status                         : LOCKED_COMPLETE")
    for name in gates:
        print(f"{name:<47}: PASS")
    print(f"Metadata                       : {metadata_path}")
    print("=" * 110)


if __name__ == "__main__":
    main()
