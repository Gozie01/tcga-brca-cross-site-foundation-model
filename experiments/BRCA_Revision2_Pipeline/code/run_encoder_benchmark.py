"""Matched UNI/Hibou benchmark; requires a complete real extraction freeze."""
from __future__ import annotations
import argparse
import hashlib
import json
import platform
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import sklearn
from sklearn.exceptions import ConvergenceWarning
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from threadpoolctl import threadpool_limits

from audit_assets import PACKAGE, sha, verify_inputs
from extract_hibou import atomic_json
import run_controls as controls
import original_stage10c as site_original


def fit_site(X_train, y_train, X_validation, y_validation, X_test):
    classes = np.unique(y_train)
    n_cv = controls.original.inner_splits(y_train)
    search = GridSearchCV(controls.estimator("std_balanced", controls.SEED),
        {"classifier__C": list(controls.original.C_GRID)},
        scoring="average_precision" if len(classes) == 2 else "neg_log_loss",
        cv=StratifiedKFold(n_cv, shuffle=True, random_state=controls.SEED),
        n_jobs=1, error_score="raise")
    search.fit(X_train, y_train)
    raw_validation = search.predict_proba(X_validation)
    raw_test = search.predict_proba(X_test)
    calibrated, calibration = site_original.calibrate(raw_validation, y_validation, raw_test, classes)
    return {"full_raw": raw_test, "validation_calibrated": calibrated}, {
        "best_C": float(search.best_params_["classifier__C"]), "calibration": calibration}


def feature_input(folder):
    metadata = json.loads((folder / "feature_freeze.json").read_text())
    if metadata.get("status") != "LOCKED_COMPLETE" or metadata.get("model") != "histai/hibou-b":
        raise ValueError("A complete real Hibou-B feature freeze is required; pilot features are ineligible")
    for file, key in [("hibou_patient_features.npy", "features_sha256"), ("hibou_patient_index.csv", "index_sha256")]:
        if sha(folder / file) != metadata[key]:
            raise ValueError("Changed second-encoder features or index")
    if metadata.get("coordinate_sha256") != verify_inputs()["files"]["patch_coordinates.csv.gz"]:
        raise ValueError("Second encoder used different patch coordinates")
    X = np.load(folder / "hibou_patient_features.npy", allow_pickle=False)
    index = pd.read_csv(folder / "hibou_patient_index.csv")
    original = pd.read_csv(PACKAGE / "data/stage10a_uni_patient_features_index_v1.csv")
    pd.testing.assert_frame_equal(index, original)
    if X.shape != (345, 768) or X.dtype != np.float32 or not np.isfinite(X).all():
        raise ValueError("Incomplete or invalid second-encoder feature matrix")
    if np.any(np.linalg.norm(X, axis=1) == 0):
        raise ValueError("Zero patient representation")
    return X, metadata


def checkpoint(out, tag, fingerprint, fit):
    array_path, meta_path = out / f"{tag}.npz", out / f"{tag}.json"
    if meta_path.exists():
        metadata = json.loads(meta_path.read_text())
        if metadata["fingerprint"] != fingerprint or metadata["predictions_sha256"] != sha(array_path):
            raise ValueError(f"Changed benchmark checkpoint: {tag}")
        with np.load(array_path, allow_pickle=False) as saved:
            prediction = {name: saved[name] for name in saved.files}
    else:
        with warnings.catch_warnings():
            warnings.simplefilter("error", ConvergenceWarning)
            prediction, details = fit()
        np.savez_compressed(array_path, **prediction)
        metadata = {"fingerprint": fingerprint, "predictions_sha256": sha(array_path), **details}
        atomic_json(meta_path, metadata)
    return prediction


def compare_group(new, reference, reps, context):
    new = new.sort_values("patient_id")
    reference = reference.sort_values("patient_id")
    if not np.array_equal(new.patient_id, reference.patient_id) or not np.array_equal(new.label, reference.label):
        raise ValueError("Encoder patient pairing failed")
    y = new.label.to_numpy(int)
    k = len(np.unique(y))
    a, b = [frame[[f"p_{c}" for c in range(k)]].to_numpy(float) for frame in [reference, new]]
    for p in [a, b]: controls.validate_probability(p, len(y), k)
    seed = controls.SEED + sum(map(ord, json.dumps(context, sort_keys=True)))
    rng = np.random.default_rng(seed)
    draws = rng.integers(0, len(y), (reps, len(y)))
    weights = np.vstack([np.ones(len(y)), np.stack([np.bincount(i, minlength=len(y)) for i in draws])]).astype(float)
    distributions = {"UNI": controls.metric_distributions(y, a, weights),
                     "Hibou-B": controls.metric_distributions(y, b, weights)}
    points, contrasts = [], []
    for encoder, dist in distributions.items():
        for metric, values in dist.items():
            valid = values[1:][np.isfinite(values[1:])]
            lo, hi = np.percentile(valid, [2.5, 97.5])
            points.append({**context, "encoder": encoder, "metric": metric, "estimate": float(values[0]),
                           "ci_lower_95": float(lo), "ci_upper_95": float(hi), "n": len(y)})
    for metric in distributions["UNI"]:
        delta = distributions["Hibou-B"][metric] - distributions["UNI"][metric]
        valid = delta[1:][np.isfinite(delta[1:])]
        lo, hi = np.percentile(valid, [2.5, 97.5])
        contrasts.append({**context, "metric": metric, "delta_hibou_minus_uni": float(delta[0]),
                          "ci_lower_95": float(lo), "ci_upper_95": float(hi), "valid_replicates": len(valid)})
    return points, contrasts


def benchmark(args):
    lock = verify_inputs()
    X, freeze = feature_input(args.features.resolve())
    data = PACKAGE / "data"
    internal = pd.read_csv(data / "stage10a_internal_5fold_manifest_v1.csv")
    sites = pd.read_csv(data / "stage10a_site_holdout_manifest_v1.csv")
    reference_internal = pd.read_csv(data / "uni_control_predictions.csv", float_precision="round_trip")
    reference_internal = reference_internal[reference_internal.variant.eq("std_balanced")]
    reference_site = pd.read_csv(data / "stage10c_site_holdout_predictions_v1.csv", float_precision="round_trip")
    out = args.out.resolve(); out.mkdir(parents=True, exist_ok=True)
    protocol = {"version": "brca_revision2_encoder_comparison_v1", "code_sha256": sha(__file__),
        "support_code_sha256": {name: sha(Path(__file__).with_name(name)) for name in
             ["run_controls.py", "original_stage10b.py", "original_stage10c.py", "extract_hibou.py", "audit_assets.py"]},
        "feature_freeze_sha256": sha(args.features / "feature_freeze.json"), "input_lock": lock,
        "bootstrap_reps": args.bootstrap, "python": platform.python_version(), "sklearn": sklearn.__version__,
        "numpy": np.__version__, "pandas": pd.__version__,
        "design": "Same coordinates, pooling, patient folds and training-only C grid; encoder-specific dimensions and official pixel normalization",
        "primary_pipeline": "Standardized balanced logistic regression; original endpoint-specific tuning criteria",
        "inference": "Exploratory fixed-prediction patient bootstrap; no multiplicity adjustment",
        "site_calibration": "Binary validation Platt; PAM50 validation temperature, matching original protocol"}
    fingerprint = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    protocol_path = out / "benchmark_protocol.json"
    if protocol_path.exists() and json.loads(protocol_path.read_text()) != protocol:
        raise ValueError("Benchmark protocol changed; use a new output directory")
    atomic_json(protocol_path, protocol)
    checkpoints = out / "checkpoints"; checkpoints.mkdir(exist_ok=True)
    prediction_rows = []
    for endpoint in controls.ENDPOINTS:
        table = internal[internal.endpoint.eq(endpoint)]
        for fold in range(5):
            train, test = table[table.internal_fold.ne(fold)], table[table.internal_fold.eq(fold)]
            def fit():
                prediction, _, details = controls.fit_control(X[train.feature_row], train.label.to_numpy(int),
                    X[test.feature_row], "std_balanced", controls.SEED + fold)
                return prediction, details
            prediction = checkpoint(checkpoints, f"{endpoint}_internal_{fold}", fingerprint, fit)
            for state, matrix in prediction.items():
                frame = test.copy(); frame["state"] = state; frame["domain"] = "internal"; frame["holdout_site"] = ""
                for c in range(matrix.shape[1]): frame[f"p_{c}"] = matrix[:, c]
                prediction_rows.append(frame)
            print(f"INTERNAL {endpoint} fold {fold} COMPLETE", flush=True)
        for site in ["A2", "C8", "E2", "BH"]:
            group = sites[sites.endpoint.eq(endpoint) & sites.holdout_site.eq(site)]
            train, val, test = [group[group.partition.eq(p)] for p in ["train", "validation", "test"]]
            if set(train.patient_id) & set(val.patient_id) or set(group[group.partition.ne('test')].patient_id) & set(test.patient_id):
                raise ValueError("Site partition separation failed")
            def fit():
                return fit_site(X[train.feature_row], train.label.to_numpy(int), X[val.feature_row],
                                val.label.to_numpy(int), X[test.feature_row])
            prediction = checkpoint(checkpoints, f"{endpoint}_site_{site}", fingerprint, fit)
            for state, matrix in prediction.items():
                frame = test.copy(); frame["state"] = state; frame["domain"] = "site"
                for c in range(matrix.shape[1]): frame[f"p_{c}"] = matrix[:, c]
                prediction_rows.append(frame)
            print(f"SITE {endpoint} {site} COMPLETE", flush=True)
    predictions = pd.concat(prediction_rows, ignore_index=True)
    points, contrasts = [], []
    for endpoint in controls.ENDPOINTS:
        for state in ["full_raw", "ensemble_raw", "ensemble_sigmoid"]:
            new = predictions[predictions.endpoint.eq(endpoint) & predictions.domain.eq("internal") & predictions.state.eq(state)]
            ref = reference_internal[reference_internal.endpoint.eq(endpoint) & reference_internal.state.eq(state)]
            p, c = compare_group(new, ref, args.bootstrap, {"endpoint": endpoint, "domain": "internal", "state": state, "holdout_site": ""})
            points += p; contrasts += c
        for site in ["A2", "C8", "E2", "BH"]:
            for state, prefix in [("full_raw", "prob_raw_class_"), ("validation_calibrated", "prob_calibrated_class_")]:
                new = predictions[predictions.endpoint.eq(endpoint) & predictions.domain.eq("site") & predictions.holdout_site.eq(site) & predictions.state.eq(state)]
                ref = reference_site[reference_site.endpoint.eq(endpoint) & reference_site.holdout_site.eq(site)].copy()
                for c in range(new.label.nunique()): ref[f"p_{c}"] = ref[f"{prefix}{c}"]
                p, c = compare_group(new, ref, args.bootstrap, {"endpoint": endpoint, "domain": "site", "state": state, "holdout_site": site})
                points += p; contrasts += c
    if predictions.duplicated(["endpoint", "domain", "holdout_site", "state", "patient_id"]).any():
        raise ValueError("Duplicate benchmark predictions")
    predictions.to_csv(out / "hibou_predictions.csv", index=False)
    pd.DataFrame(points).to_csv(out / "encoder_metrics.csv", index=False)
    pd.DataFrame(contrasts).to_csv(out / "paired_encoder_contrasts.csv", index=False)
    atomic_json(out / "benchmark_complete.json", {"status": "COMPLETE", "fingerprint": fingerprint,
        "internal_endpoint_fold_fits": 25, "site_fits": 20, "bootstrap_reps": args.bootstrap,
        "outputs": {p.name: sha(p) for p in out.glob("*.csv")}})
    print("COMPLETE: paired encoder metrics and contrasts written.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--features", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    if args.bootstrap < 100: parser.error("Use at least 100 bootstrap replicates")
    with threadpool_limits(limits=1): benchmark(args)
