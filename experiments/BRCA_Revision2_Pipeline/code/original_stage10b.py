from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import time
from pathlib import Path

import joblib
import numpy as np
import pandas as pd
import sklearn
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (
    average_precision_score, balanced_accuracy_score, brier_score_loss,
    confusion_matrix, f1_score, log_loss, matthews_corrcoef, recall_score,
    roc_auc_score,
)
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler, label_binarize


VERSION = "stage10b_internal_prediction_calibration_v1"
SEED = 20260905
C_GRID = (0.001, 0.01, 0.1, 1.0, 10.0, 100.0)
EXPECTED_PATIENTS = 345
EXPECTED_DIM = 1024
EXPECTED_ENDPOINT_N = {
    "ER_binary": 337,
    "PR_binary": 337,
    "HER2_binary": 314,
    "PIK3CA_binary": 345,
    "PAM50_primary_4class": 335,
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


def expected_calibration_error(y: np.ndarray, probabilities: np.ndarray,
                               bins: int = 10) -> float:
    if probabilities.ndim == 1 or probabilities.shape[1] == 2:
        p = probabilities if probabilities.ndim == 1 else probabilities[:, 1]
        confidence = np.maximum(p, 1.0 - p)
        predictions = (p >= 0.5).astype(int)
    else:
        confidence = probabilities.max(axis=1)
        predictions = probabilities.argmax(axis=1)
    correct = predictions == y
    edges = np.linspace(0.0, 1.0, bins + 1)
    ece = 0.0
    for i in range(bins):
        if i == bins - 1:
            mask = (confidence >= edges[i]) & (confidence <= edges[i + 1])
        else:
            mask = (confidence >= edges[i]) & (confidence < edges[i + 1])
        if mask.any():
            ece += mask.mean() * abs(correct[mask].mean() - confidence[mask].mean())
    return float(ece)


def calculate_metrics(y: np.ndarray, probabilities: np.ndarray,
                      classes: np.ndarray) -> dict:
    n_classes = len(classes)
    predicted_positions = probabilities.argmax(axis=1)
    predictions = classes[predicted_positions]
    result = {
        "auroc": np.nan, "auprc": np.nan, "balanced_accuracy": np.nan,
        "macro_f1": np.nan, "mcc": np.nan, "sensitivity": np.nan,
        "specificity": np.nan, "brier": np.nan, "log_loss": np.nan,
        "ece_10bin": np.nan,
    }
    if n_classes == 2:
        positive_position = int(np.where(classes == 1)[0][0])
        p = probabilities[:, positive_position]
        predictions = (p >= 0.5).astype(int)
        tn, fp, fn, tp = confusion_matrix(y, predictions, labels=[0, 1]).ravel()
        result.update({
            "auroc": roc_auc_score(y, p),
            "auprc": average_precision_score(y, p),
            "balanced_accuracy": balanced_accuracy_score(y, predictions),
            "macro_f1": f1_score(y, predictions, average="macro", zero_division=0),
            "mcc": matthews_corrcoef(y, predictions),
            "sensitivity": recall_score(y, predictions, pos_label=1, zero_division=0),
            "specificity": tn / (tn + fp) if tn + fp else np.nan,
            "brier": brier_score_loss(y, p),
            "log_loss": log_loss(y, probabilities, labels=classes),
            "ece_10bin": expected_calibration_error(y, probabilities),
            "tn": int(tn), "fp": int(fp), "fn": int(fn), "tp": int(tp),
        })
    else:
        onehot = label_binarize(y, classes=classes)
        result.update({
            "auroc": roc_auc_score(y, probabilities, labels=classes,
                                    multi_class="ovr", average="macro"),
            "auprc": average_precision_score(onehot, probabilities, average="macro"),
            "balanced_accuracy": balanced_accuracy_score(y, predictions),
            "macro_f1": f1_score(y, predictions, average="macro", zero_division=0),
            "mcc": matthews_corrcoef(y, predictions),
            "brier": float(np.mean(np.sum((probabilities - onehot) ** 2, axis=1))),
            "log_loss": log_loss(y, probabilities, labels=classes),
            "ece_10bin": expected_calibration_error(y, probabilities),
        })
    return {k: (float(v) if isinstance(v, (np.floating, float)) else v)
            for k, v in result.items()}


def inner_splits(y: np.ndarray, maximum: int = 4) -> int:
    minimum_class = int(pd.Series(y).value_counts().min())
    value = min(maximum, minimum_class)
    if value < 2:
        raise RuntimeError("Insufficient class support for inner cross-validation")
    return value


def build_estimator(seed: int) -> Pipeline:
    return Pipeline([
        ("scaler", StandardScaler()),
        ("classifier", LogisticRegression(
            C=1.0, class_weight="balanced", solver="lbfgs",
            max_iter=5000, random_state=seed,
        )),
    ])


def task_paths(checkpoint_dir: Path, endpoint: str, fold: int):
    stem = f"{endpoint}_fold{fold}"
    return {
        "prediction": checkpoint_dir / f"{stem}_predictions.csv",
        "metrics": checkpoint_dir / f"{stem}_metrics.csv",
        "raw_model": checkpoint_dir / f"{stem}_raw_model.joblib",
        "calibrated_model": checkpoint_dir / f"{stem}_calibrated_model.joblib",
        "metadata": checkpoint_dir / f"{stem}_checkpoint.json",
    }


def checkpoint_valid(paths: dict[str, Path], fingerprint: str) -> bool:
    if not all(path.is_file() for path in paths.values()):
        return False
    try:
        meta = json.loads(paths["metadata"].read_text(encoding="utf-8"))
        return (
            meta.get("status") == "CHECKPOINT_COMPLETE"
            and meta.get("fingerprint") == fingerprint
            and sha256_file(paths["prediction"]) == meta.get("prediction_sha256")
            and sha256_file(paths["metrics"]) == meta.get("metrics_sha256")
            and sha256_file(paths["raw_model"]) == meta.get("raw_model_sha256")
            and sha256_file(paths["calibrated_model"]) == meta.get("calibrated_model_sha256")
        )
    except Exception:
        return False


def fit_fold(X: np.ndarray, split: pd.DataFrame, endpoint: str, fold: int,
             checkpoint_dir: Path, fingerprint: str, n_jobs: int):
    paths = task_paths(checkpoint_dir, endpoint, fold)
    if checkpoint_valid(paths, fingerprint):
        return pd.read_csv(paths["prediction"]), pd.read_csv(paths["metrics"]), True
    train = split[split["internal_fold"].ne(fold)].copy()
    test = split[split["internal_fold"].eq(fold)].copy()
    train_rows = train["feature_row"].astype(int).to_numpy()
    test_rows = test["feature_row"].astype(int).to_numpy()
    y_train = train["label"].astype(int).to_numpy()
    y_test = test["label"].astype(int).to_numpy()
    classes = np.sort(split["label"].astype(int).unique())
    seed = SEED + fold
    cv_n = inner_splits(y_train)
    scoring = "average_precision" if len(classes) == 2 else "neg_log_loss"
    search = GridSearchCV(
        build_estimator(seed),
        {"classifier__C": list(C_GRID)},
        scoring=scoring,
        cv=StratifiedKFold(cv_n, shuffle=True, random_state=seed),
        n_jobs=n_jobs,
        refit=True,
        error_score="raise",
    )
    started = time.perf_counter()
    search.fit(X[train_rows], y_train)
    raw_model = search.best_estimator_
    calibration_cv = StratifiedKFold(cv_n, shuffle=True, random_state=seed + 1000)
    calibrated_model = CalibratedClassifierCV(
        estimator=clone(raw_model), method="sigmoid", cv=calibration_cv,
        n_jobs=n_jobs,
    )
    calibrated_model.fit(X[train_rows], y_train)
    raw_p = raw_model.predict_proba(X[test_rows])
    calibrated_p = calibrated_model.predict_proba(X[test_rows])
    raw_metrics = calculate_metrics(y_test, raw_p, classes)
    cal_metrics = calculate_metrics(y_test, calibrated_p, classes)
    elapsed = time.perf_counter() - started
    prediction = test[[
        "endpoint", "patient_id", "feature_row", "source", "tss",
        "label", "label_name", "internal_fold",
    ]].copy().reset_index(drop=True)
    for position, label in enumerate(classes):
        prediction[f"prob_raw_class_{label}"] = raw_p[:, position]
        prediction[f"prob_calibrated_class_{label}"] = calibrated_p[:, position]
    metric_rows = []
    for calibration, values in (("raw", raw_metrics), ("sigmoid_cv", cal_metrics)):
        metric_rows.append({
            "endpoint": endpoint, "fold": fold, "calibration": calibration,
            "train_n": len(train), "test_n": len(test),
            "inner_cv_splits": cv_n, "best_C": search.best_params_["classifier__C"],
            "elapsed_seconds": elapsed, **values,
        })
    metrics = pd.DataFrame(metric_rows)
    atomic_csv(paths["prediction"], prediction)
    atomic_csv(paths["metrics"], metrics)
    raw_tmp = paths["raw_model"].with_suffix(".joblib.tmp")
    cal_tmp = paths["calibrated_model"].with_suffix(".joblib.tmp")
    joblib.dump(raw_model, raw_tmp)
    joblib.dump(calibrated_model, cal_tmp)
    os.replace(raw_tmp, paths["raw_model"])
    os.replace(cal_tmp, paths["calibrated_model"])
    meta = {
        "version": VERSION, "status": "CHECKPOINT_COMPLETE",
        "fingerprint": fingerprint, "endpoint": endpoint, "fold": fold,
        "prediction_sha256": sha256_file(paths["prediction"]),
        "metrics_sha256": sha256_file(paths["metrics"]),
        "raw_model_sha256": sha256_file(paths["raw_model"]),
        "calibrated_model_sha256": sha256_file(paths["calibrated_model"]),
    }
    atomic_json(paths["metadata"], meta)
    if not checkpoint_valid(paths, fingerprint):
        raise RuntimeError(f"Checkpoint verification failed: {endpoint}/fold{fold}")
    return prediction, metrics, False


def bootstrap_intervals(predictions: pd.DataFrame, endpoint: str,
                        reps: int) -> pd.DataFrame:
    table = predictions[predictions["endpoint"].eq(endpoint)].copy()
    y = table["label"].astype(int).to_numpy()
    classes = np.sort(np.unique(y))
    records = []
    for calibration, prefix in (("raw", "prob_raw_class_"),
                                ("sigmoid_cv", "prob_calibrated_class_")):
        probabilities = np.column_stack([
            table[f"{prefix}{label}"].to_numpy(float) for label in classes
        ])
        point = calculate_metrics(y, probabilities, classes)
        rng = np.random.default_rng(SEED + sum(map(ord, endpoint + calibration)))
        samples = {name: [] for name in (
            "auroc", "auprc", "balanced_accuracy", "macro_f1", "mcc",
            "sensitivity", "specificity", "brier", "log_loss", "ece_10bin",
        )}
        for _ in range(reps):
            indices = rng.integers(0, len(y), len(y))
            if len(np.unique(y[indices])) != len(classes):
                continue
            try:
                values = calculate_metrics(y[indices], probabilities[indices], classes)
            except ValueError:
                continue
            for name in samples:
                value = values.get(name, np.nan)
                if np.isfinite(value):
                    samples[name].append(value)
        for name, values in samples.items():
            if not values or not np.isfinite(point.get(name, np.nan)):
                continue
            records.append({
                "endpoint": endpoint, "calibration": calibration, "metric": name,
                "estimate": point[name],
                "ci_lower_95": float(np.percentile(values, 2.5)),
                "ci_upper_95": float(np.percentile(values, 97.5)),
                "bootstrap_valid_replicates": len(values),
            })
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=Path.cwd())
    parser.add_argument("--bootstrap-reps", type=int, default=1000)
    parser.add_argument("--n-jobs", type=int, default=-1)
    args = parser.parse_args()
    root = args.root
    stage10a = root / "stage10a_cohort_split_freeze_v1"
    features_path = stage10a / "stage10a_uni_patient_features_v1.npy"
    feature_index_path = stage10a / "stage10a_uni_patient_features_index_v1.csv"
    split_path = stage10a / "stage10a_internal_5fold_manifest_v1.csv"
    freeze_meta_path = stage10a / "stage10a_cohort_split_freeze_metadata_v1.json"
    for path in (features_path, feature_index_path, split_path, freeze_meta_path):
        if not path.is_file():
            raise FileNotFoundError(f"Missing frozen Stage 10A input: {path}")
    freeze_meta = json.loads(freeze_meta_path.read_text(encoding="utf-8"))
    if freeze_meta.get("status") != "LOCKED_COMPLETE":
        raise RuntimeError("Stage 10A is not LOCKED_COMPLETE")
    upstream = freeze_meta.get("outputs", {})
    for path in (features_path, feature_index_path, split_path):
        expected_hash = upstream.get(path.name, {}).get("sha256")
        if not expected_hash or sha256_file(path) != expected_hash:
            raise RuntimeError(f"Stage 10A hash mismatch: {path.name}")
    X = np.load(features_path, mmap_mode="r", allow_pickle=False)
    feature_index = pd.read_csv(feature_index_path)
    splits = pd.read_csv(split_path)
    gates = {
        "feature_shape": X.shape == (EXPECTED_PATIENTS, EXPECTED_DIM),
        "feature_dtype_float32": X.dtype == np.float32,
        "feature_index_rows": len(feature_index) == EXPECTED_PATIENTS,
        "five_folds_per_endpoint": splits.groupby("endpoint")["internal_fold"].nunique().eq(5).all(),
        "endpoint_counts": splits.groupby("endpoint")["patient_id"].nunique().to_dict()
                           == EXPECTED_ENDPOINT_N,
        "unique_endpoint_patients": not splits.duplicated(["endpoint", "patient_id"]).any(),
    }
    failed = [name for name, passed in gates.items() if not bool(passed)]
    if failed:
        raise RuntimeError(f"Stage 10B input gates failed: {failed}")

    out = root / "stage10b_internal_prediction_calibration_v1"
    checkpoint_dir = out / "split_checkpoints"
    checkpoint_dir.mkdir(parents=True, exist_ok=True)
    fingerprint_base = hashlib.sha256(
        (sha256_file(features_path) + sha256_file(split_path) + VERSION).encode("utf-8")
    ).hexdigest()
    predictions_all = []
    metrics_all = []
    print("=" * 110)
    print("STAGE 10B INTERNAL PREDICTION AND CALIBRATION")
    print(f"Patients={EXPECTED_PATIENTS} | features={X.shape} | bootstrap={args.bootstrap_reps}")
    print("=" * 110)
    for endpoint in EXPECTED_ENDPOINT_N:
        endpoint_split = splits[splits["endpoint"].eq(endpoint)].copy()
        print(f"\n{endpoint} | evaluable={len(endpoint_split)}")
        for fold in range(5):
            fingerprint = hashlib.sha256(
                f"{fingerprint_base}|{endpoint}|{fold}".encode("utf-8")
            ).hexdigest()
            prediction, metrics, skipped = fit_fold(
                X, endpoint_split, endpoint, fold, checkpoint_dir,
                fingerprint, args.n_jobs,
            )
            predictions_all.append(prediction)
            metrics_all.append(metrics)
            status = "VERIFIED SKIP" if skipped else "CHECKPOINT PASS"
            raw = metrics[metrics["calibration"].eq("raw")].iloc[0]
            cal = metrics[metrics["calibration"].eq("sigmoid_cv")].iloc[0]
            print(f"  fold {fold}: {status} | raw AUROC={raw.auroc:.4f} "
                  f"cal AUROC={cal.auroc:.4f} | raw ECE={raw.ece_10bin:.4f} "
                  f"cal ECE={cal.ece_10bin:.4f}")

    predictions = pd.concat(predictions_all, ignore_index=True)
    fold_metrics = pd.concat(metrics_all, ignore_index=True)
    if predictions.duplicated(["endpoint", "patient_id"]).any():
        raise RuntimeError("OOF prediction duplication detected")
    if predictions.groupby("endpoint")["patient_id"].nunique().to_dict() != EXPECTED_ENDPOINT_N:
        raise RuntimeError("OOF prediction count mismatch")
    bootstrap = pd.concat([
        bootstrap_intervals(predictions, endpoint, args.bootstrap_reps)
        for endpoint in EXPECTED_ENDPOINT_N
    ], ignore_index=True)
    summary = (
        fold_metrics.groupby(["endpoint", "calibration"])[[
            "auroc", "auprc", "balanced_accuracy", "macro_f1", "mcc",
            "sensitivity", "specificity", "brier", "log_loss", "ece_10bin",
        ]].agg(["mean", "std"]).reset_index()
    )
    summary.columns = [
        "_".join(str(x) for x in col if str(x)) if isinstance(col, tuple) else col
        for col in summary.columns
    ]
    predictions_out = out / "stage10b_internal_oof_predictions_v1.csv"
    fold_metrics_out = out / "stage10b_internal_fold_metrics_v1.csv"
    summary_out = out / "stage10b_internal_summary_metrics_v1.csv"
    bootstrap_out = out / "stage10b_internal_bootstrap_ci_v1.csv"
    metadata_out = out / "stage10b_internal_prediction_calibration_metadata_v1.json"
    atomic_csv(predictions_out, predictions)
    atomic_csv(fold_metrics_out, fold_metrics)
    atomic_csv(summary_out, summary)
    atomic_csv(bootstrap_out, bootstrap)
    final_gates = {
        **{f"input_{k}": bool(v) for k, v in gates.items()},
        "twenty_five_split_checkpoints": len(list(checkpoint_dir.glob("*_checkpoint.json"))) == 25,
        "complete_oof_predictions": predictions.groupby("endpoint")["patient_id"].nunique().to_dict()
                                    == EXPECTED_ENDPOINT_N,
        "no_oof_duplicates": not predictions.duplicated(["endpoint", "patient_id"]).any(),
        "two_calibration_states_per_fold": fold_metrics.groupby(
            ["endpoint", "fold"]
        )["calibration"].nunique().eq(2).all(),
        "all_primary_metrics_finite": np.isfinite(
            fold_metrics[["auroc", "auprc", "balanced_accuracy", "macro_f1",
                          "mcc", "brier", "log_loss", "ece_10bin"]].to_numpy(float)
        ).all(),
        "bootstrap_rows_present": len(bootstrap) > 0,
    }
    failed = [name for name, passed in final_gates.items() if not bool(passed)]
    if failed:
        raise RuntimeError(f"Stage 10B final gates failed: {failed}")
    outputs = [predictions_out, fold_metrics_out, summary_out, bootstrap_out]
    metadata = {
        "pipeline_version": VERSION,
        "status": "LOCKED_COMPLETE",
        "seed": SEED,
        "model": "L2-regularized logistic-regression linear probe",
        "preprocessing": "StandardScaler fitted within each training fold",
        "hyperparameter_selection": {
            "C_grid": list(C_GRID), "inner_cv_max_folds": 4,
            "binary_metric": "average_precision",
            "multiclass_metric": "negative_log_loss",
        },
        "calibration": (
            "sigmoid calibration fitted by stratified CV inside each outer-training set; "
            "outer test folds untouched"
        ),
        "evaluation": "frozen endpoint-specific patient-level five-fold OOF",
        "bootstrap_replicates": args.bootstrap_reps,
        "endpoint_evaluable_counts": EXPECTED_ENDPOINT_N,
        "software": {
            "python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "scikit_learn": sklearn.__version__,
        },
        "stage10a_metadata_sha256": sha256_file(freeze_meta_path),
        "validation": {name: bool(passed) for name, passed in final_gates.items()},
        "outputs": {
            path.name: {"path": str(path), "sha256": sha256_file(path)}
            for path in outputs
        },
    }
    atomic_json(metadata_out, metadata)
    print("\n" + "=" * 110)
    print("STAGE 10B INTERNAL PREDICTION/CALIBRATION AUDIT")
    print("Status                         : LOCKED_COMPLETE")
    print(f"Endpoints                      : {len(EXPECTED_ENDPOINT_N)}")
    print("Outer folds per endpoint       : 5")
    print(f"Unique endpoint-patient OOF    : {len(predictions):,}")
    for name, passed in final_gates.items():
        print(f"{name:<47}: {'PASS' if passed else 'FAIL'}")
    print(f"Metadata                       : {metadata_out}")
    print("=" * 110)


if __name__ == "__main__":
    main()
