"""Controlled reanalysis of the released TCGA-BRCA patient representations.

No image extraction, new partitioning, model selection on test patients, or
replacement of original results. See protocol.json and README.md.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scipy
import sklearn
from scipy.optimize import brentq, minimize
from scipy.special import expit, logit
from sklearn.base import clone
from sklearn.calibration import CalibratedClassifierCV
from sklearn.exceptions import ConvergenceWarning
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import Normalizer, StandardScaler
from threadpoolctl import threadpool_limits

import original_stage10b as original

VERSION = "brca_revision1_controls_v1"
SEED = original.SEED
ENDPOINTS = tuple(original.EXPECTED_ENDPOINT_N)
PARITY_TOLERANCE = 5e-5
VARIANTS = (
    "std_balanced", "std_unweighted", "none_balanced",
    "l2_balanced", "l2_std_balanced",
)


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False), encoding="utf-8")


def verify_table(folder, metadata_name, names):
    metadata = json.loads((folder / metadata_name).read_text(encoding="utf-8-sig"))
    if metadata.get("status") != "LOCKED_COMPLETE":
        raise ValueError(f"Upstream is not locked: {folder}")
    hashes = {}
    for name in names:
        actual = digest(folder / name)
        expected = metadata["outputs"][name]["sha256"]
        if actual != expected:
            raise ValueError(f"Frozen input changed: {name}")
        hashes[name] = actual
    hashes[metadata_name] = digest(folder / metadata_name)
    return hashes


def load_inputs(root):
    a = root / "stage10a_cohort_split_freeze_v1"
    b = root / "stage10b_internal_prediction_calibration_v1"
    c = root / "stage10c_site_holdout_domain_shift_v1"
    anames = [p.name for p in a.glob("*.csv")] + ["stage10a_uni_patient_features_v1.npy"]
    hashes = verify_table(a, "stage10a_cohort_split_freeze_metadata_v1.json", anames)
    hashes.update(verify_table(b, "stage10b_internal_prediction_calibration_metadata_v1.json", [
        "stage10b_internal_oof_predictions_v1.csv", "stage10b_internal_fold_metrics_v1.csv",
    ]))
    hashes.update(verify_table(c, "stage10c_site_holdout_metadata_v1.json", [
        "stage10c_site_holdout_predictions_v1.csv", "stage10c_site_holdout_metrics_v1.csv",
    ]))
    X = np.load(a / "stage10a_uni_patient_features_v1.npy", allow_pickle=False)
    index = pd.read_csv(a / "stage10a_uni_patient_features_index_v1.csv")
    splits = pd.read_csv(a / "stage10a_internal_5fold_manifest_v1.csv")
    sites = pd.read_csv(a / "stage10a_site_holdout_manifest_v1.csv")
    if X.shape != (345, 1024) or X.dtype != np.float32 or not np.isfinite(X).all():
        raise ValueError("Unexpected or non-finite feature matrix")
    if index.patient_id.duplicated().any() or index.feature_row.duplicated().any():
        raise ValueError("Duplicate patient features")
    if not np.array_equal(index.feature_row.to_numpy(), np.arange(len(X))):
        raise ValueError("Feature index ordering changed")
    for table, keys in [(splits, ["endpoint", "patient_id"]),
                        (sites, ["endpoint", "holdout_site", "patient_id"])]:
        if table.duplicated(keys).any():
            raise ValueError("Duplicate partition assignments")
        mapped = index.set_index("feature_row").loc[table.feature_row, "patient_id"].to_numpy()
        if not np.array_equal(mapped, table.patient_id.to_numpy()):
            raise ValueError("Patient-to-feature mapping changed")
        if not np.array_equal(index.set_index("feature_row").loc[table.feature_row, "tss"], table.tss):
            raise ValueError("Site mapping changed")
    if splits.groupby("endpoint").size().to_dict() != original.EXPECTED_ENDPOINT_N:
        raise ValueError("Endpoint counts changed")
    for endpoint in ENDPOINTS:
        d = splits[splits.endpoint.eq(endpoint)]
        if set(d.internal_fold) != set(range(5)):
            raise ValueError("Frozen outer folds changed")
        wanted = set(range(4)) if endpoint.startswith("PAM50") else {0, 1}
        if set(d.label) != wanted:
            raise ValueError("Unexpected class encoding")
    old = pd.read_csv(b / "stage10b_internal_oof_predictions_v1.csv")
    old_metrics = pd.read_csv(b / "stage10b_internal_fold_metrics_v1.csv")
    external = pd.read_csv(c / "stage10c_site_holdout_predictions_v1.csv")
    expected = splits.set_index(["endpoint", "patient_id"])[["feature_row", "label", "internal_fold"]].sort_index()
    found = old.set_index(["endpoint", "patient_id"])[expected.columns].sort_index()
    pd.testing.assert_frame_equal(expected, found)
    return X, index, splits, sites, old, old_metrics, external, hashes


def estimator(variant, seed):
    steps = []
    if variant.startswith("l2"):
        steps.append(("normalize", Normalizer(norm="l2")))
    if variant.startswith("std") or variant.startswith("l2_std"):
        steps.append(("scaler", StandardScaler()))
    steps.append(("classifier", LogisticRegression(
        class_weight=None if variant.endswith("unweighted") else "balanced",
        solver="lbfgs", max_iter=5000, random_state=seed,
    )))
    return Pipeline(steps)


def validate_probability(p, n, k):
    if p.shape != (n, k) or not np.isfinite(p).all():
        raise ValueError("Invalid prediction matrix")
    if (p < 0).any() or (p > 1).any() or not np.allclose(p.sum(axis=1), 1, atol=1e-10):
        raise ValueError("Invalid probability simplex")


def fit_control(X_train, y_train, X_test, variant, seed, fixed_C=None):
    """Outer-test labels are deliberately not accepted by this function."""
    classes = np.unique(y_train)
    n_cv = original.inner_splits(y_train)
    model = estimator(variant, seed)
    if fixed_C is None:
        search = GridSearchCV(model, {"classifier__C": list(original.C_GRID)},
            scoring="average_precision" if len(classes) == 2 else "neg_log_loss",
            cv=StratifiedKFold(n_cv, shuffle=True, random_state=seed),
            refit=True, n_jobs=1, error_score="raise")
        search.fit(X_train, y_train)
        raw = search.best_estimator_
        best_C = search.best_params_["classifier__C"]
    else:
        best_C = float(fixed_C)
        raw = model.set_params(classifier__C=best_C).fit(X_train, y_train)
    predictions = {"full_raw": raw.predict_proba(X_test)}
    inner_records, components = [], {}
    if variant in ("std_balanced", "std_unweighted"):
        cv = list(StratifiedKFold(n_cv, shuffle=True, random_state=seed + 1000).split(X_train, y_train))
        calibrated = CalibratedClassifierCV(estimator=clone(raw), method="sigmoid",
                                           cv=cv, ensemble=True, n_jobs=1)
        calibrated.fit(X_train, y_train)
        raw_components, calibrated_components = [], []
        for (train_i, cal_i), pair in zip(cv, calibrated.calibrated_classifiers_):
            if set(train_i) & set(cal_i) or len(train_i) + len(cal_i) != len(y_train):
                raise ValueError("Inner partition separation failed")
            if not np.array_equal(pair.estimator.classes_, classes):
                raise ValueError("Inner model class ordering changed")
            raw_components.append(pair.estimator.predict_proba(X_test))
            calibrated_components.append(pair.predict_proba(X_test))
            inner_records.append({"fit_local_rows": train_i.tolist(), "calibration_local_rows": cal_i.tolist()})
        components = {"raw_components": np.stack(raw_components),
                      "calibrated_components": np.stack(calibrated_components)}
        predictions["ensemble_raw"] = components["raw_components"].mean(axis=0)
        predictions["ensemble_sigmoid"] = components["calibrated_components"].mean(axis=0)
        if not np.allclose(predictions["ensemble_sigmoid"], calibrated.predict_proba(X_test), atol=1e-12):
            raise ValueError("Calibrated ensemble reconstruction failed")
    for p in predictions.values():
        validate_probability(p, len(X_test), len(classes))
    return predictions, components, {"best_C": float(best_C), "inner_splits": n_cv,
                                    "inner_partition_rows": inner_records}


def prevalence_predictions(y_development, n_test, classes):
    counts = np.array([(y_development == c).sum() for c in classes], dtype=float)
    if not counts.sum() or not (counts > 0).all():
        raise ValueError("Prevalence reference lacks class support")
    return np.tile(counts / counts.sum(), (n_test, 1))


def metric_distributions(y, p, weights):
    """Paired ordinary patient bootstrap, with exact weighted tie handling."""
    n, k = p.shape
    onehot = np.eye(k)[y]
    losses = ((p - onehot) ** 2).sum(axis=1)
    if k == 2:
        losses /= 2  # manuscript binary Brier convention
    total = weights.sum(axis=1)
    result = {"brier": weights @ losses / total,
              "log_loss": weights @ (-np.log(np.clip(p[np.arange(n), y],
                  np.finfo(p.dtype).eps, 1 - np.finfo(p.dtype).eps))) / total}
    # Match the manuscript's binary 0.5 tie rule and multiclass argmax.
    predicted = (p[:, 1] >= .5).astype(int) if k == 2 else p.argmax(axis=1)
    confidence = p.max(axis=1)
    # Match the released np.linspace boundaries, including floating-point edges.
    bins = np.clip(np.searchsorted(np.linspace(0., 1., 11), confidence, side="right") - 1, 0, 9)
    result["ece_10bin"] = sum(np.abs(weights @ ((bins == b) * ((predicted == y) - confidence)))
                             for b in range(10)) / total
    aucs, aps = [], []
    for c in ([1] if k == 2 else range(k)):
        order = np.argsort(p[:, c], kind="stable")
        starts = np.r_[0, np.flatnonzero(np.diff(p[order, c])) + 1]
        positive = np.add.reduceat(weights[:, order] * (y[order] == c), starts, axis=1)
        negative = np.add.reduceat(weights[:, order] * (y[order] != c), starts, axis=1)
        tp, tn = positive.sum(axis=1), negative.sum(axis=1)
        with np.errstate(divide="ignore", invalid="ignore"):
            auc = (positive * (np.cumsum(negative, axis=1) - negative / 2)).sum(axis=1) / (tp * tn)
            pos_d, neg_d = positive[:, ::-1], negative[:, ::-1]
            cum_pos = np.cumsum(pos_d, axis=1)
            cum_all = cum_pos + np.cumsum(neg_d, axis=1)
            precision = np.divide(cum_pos, cum_all, out=np.zeros_like(cum_pos), where=cum_all > 0)
            ap = (precision * pos_d).sum(axis=1) / tp
        aucs.append(auc); aps.append(ap)
    result["auroc"] = np.mean(aucs, axis=0)
    result["auprc"] = np.mean(aps, axis=0)
    return result


def calibration_diagnostics(y, p):
    rows = []
    k = p.shape[1]
    for c in ([1] if k == 2 else range(k)):
        target = (y == c).astype(float)
        z = logit(np.clip(p[:, c], 1e-6, 1 - 1e-6))
        alpha = brentq(lambda a: (expit(z + a) - target).sum(), -100, 100)
        def objective(ab):
            eta = ab[0] + ab[1] * z
            error = expit(eta) - target
            return np.logaddexp(0, eta).sum() - target @ eta, np.array([error.sum(), error @ z])
        fitted = minimize(objective, [0., 1.], jac=True, method="L-BFGS-B",
                          options={"maxiter": 2000, "ftol": 1e-13, "gtol": 1e-8})
        if not fitted.success:
            raise ValueError(f"Calibration diagnostic fit failed: {fitted.message}")
        rows.append({"class": c, "calibration_in_the_large": float(alpha),
                     "recalibration_intercept": float(fitted.x[0]), "calibration_slope": float(fitted.x[1])})
    return rows


def run(args):
    out = args.out.resolve()
    out.mkdir(parents=True, exist_ok=True)
    X, index, splits, sites, old, old_metrics, external, hashes = load_inputs(args.root.resolve())
    protocol = {"version": VERSION, "seed": SEED, "C_grid": list(original.C_GRID),
        "outer_partitions": "Unmodified released endpoint-specific patient folds",
        "calibration": "Explicit four-fold ensemble; same component estimators for raw and calibrated means",
        "balanced_C": "Reuse released outer-training-only tuning; verify prediction parity",
        "unweighted_and_preprocessing_C": "Retune in outer training with original grid and selection criterion",
        "preprocessing": "2x2 sample L2 (off/on) by feature standardization (off/on), balanced LR",
        "prevalence_reference": "Outer-training only internally; training+validation development only for site test",
        "primary_comparisons": "Paired changes in Brier and log loss; AUROC, AP and top-label ECE also reported",
        "bootstrap": "1000 ordinary patient resamples, fixed predictions; no retraining uncertainty",
        "multiplicity": "Exploratory, no multiplicity adjustment; percentile intervals, no confirmatory claims",
        "clinical_utility": "No clinical utility, external validity or assay-replacement claim",
        "input_hashes": hashes, "code_sha256": digest(__file__),
        "original_code_sha256": digest(Path(__file__).with_name("original_stage10b.py")),
        "software": {"python": platform.python_version(), "numpy": np.__version__,
            "pandas": pd.__version__, "sklearn": sklearn.__version__, "scipy": scipy.__version__},
        "bootstrap_reps": args.bootstrap, "parity_probability_tolerance": PARITY_TOLERANCE}
    protocol["bootstrap"] = f"{args.bootstrap} ordinary patient resamples, fixed predictions; no retraining uncertainty"
    protocol_hash = hashlib.sha256(json.dumps(protocol, sort_keys=True).encode()).hexdigest()
    protocol_path = out / "protocol.json"
    if protocol_path.exists() and json.loads(protocol_path.read_text()) != protocol:
        raise ValueError("Existing run protocol differs; choose a fresh output directory")
    write_json(protocol_path, protocol)  # frozen before new predictions or outcome summaries
    checkpoints = out / "checkpoints"
    checkpoints.mkdir(exist_ok=True)
    all_rows, parity = [], []
    for endpoint in ENDPOINTS:
        table = splits[splits.endpoint.eq(endpoint)].copy()
        classes = np.sort(table.label.unique())
        for fold in range(5):
            train, test = table[table.internal_fold.ne(fold)], table[table.internal_fold.eq(fold)]
            y_train = train.label.to_numpy(int)
            for variant in VARIANTS:
                tag = f"{endpoint}_fold{fold}_{variant}"
                npz_path, meta_path = checkpoints / f"{tag}.npz", checkpoints / f"{tag}.json"
                if meta_path.exists():
                    metadata = json.loads(meta_path.read_text())
                    if metadata["protocol_sha256"] != protocol_hash or metadata["npz_sha256"] != digest(npz_path):
                        raise ValueError(f"Checkpoint mismatch: {tag}")
                    with np.load(npz_path, allow_pickle=False) as saved:
                        prediction = {k.removeprefix("prediction_"): saved[k] for k in saved.files if k.startswith("prediction_")}
                else:
                    print(f"FIT {tag}", flush=True)
                    before = time.perf_counter()
                    fixed_C = None
                    if variant == "std_balanced":
                        row = old_metrics[old_metrics.endpoint.eq(endpoint) & old_metrics.fold.eq(fold) & old_metrics.calibration.eq("raw")]
                        if len(row) != 1: raise ValueError("Ambiguous released C")
                        fixed_C = row.best_C.iloc[0]
                    with warnings.catch_warnings():
                        warnings.simplefilter("error", ConvergenceWarning)
                        prediction, components, details = fit_control(
                            X[train.feature_row], y_train, X[test.feature_row], variant, SEED + fold, fixed_C)
                    np.savez_compressed(npz_path, **{f"prediction_{k}": v for k, v in prediction.items()}, **components)
                    metadata = {"protocol_sha256": protocol_hash, "npz_sha256": digest(npz_path),
                        "endpoint": endpoint, "fold": fold, "variant": variant,
                        "train_patient_ids": train.patient_id.tolist(), "test_patient_ids": test.patient_id.tolist(),
                        "elapsed_seconds": time.perf_counter() - before, **details}
                    write_json(meta_path, metadata)
                if metadata["train_patient_ids"] != train.patient_id.tolist() or metadata["test_patient_ids"] != test.patient_id.tolist():
                    raise ValueError("Checkpoint patient identity changed")
                if variant == "std_balanced":
                    prior = old[old.endpoint.eq(endpoint)].set_index("patient_id").loc[test.patient_id]
                    for state, prefix in [("full_raw", "prob_raw_class_"), ("ensemble_sigmoid", "prob_calibrated_class_")]:
                        original_p = prior[[f"{prefix}{c}" for c in classes]].to_numpy()
                        delta = float(np.max(np.abs(original_p - prediction[state])))
                        parity.append({"endpoint": endpoint, "fold": fold, "state": state, "max_abs_probability_difference": delta,
                                       "pass": delta <= PARITY_TOLERANCE})
                        if delta > PARITY_TOLERANCE:
                            raise ValueError(f"Released-prediction parity failed: {tag}/{state}: {delta}")
                for state, p in prediction.items():
                    validate_probability(p, len(test), len(classes))
                    rows = test.copy()
                    rows["variant"], rows["state"] = variant, state
                    for j, c in enumerate(classes): rows[f"p_{c}"] = p[:, j]
                    all_rows.append(rows)
            rows = test.copy()
            rows["variant"], rows["state"] = "prevalence", "development_only"
            p = prevalence_predictions(y_train, len(test), classes)
            for j, c in enumerate(classes): rows[f"p_{c}"] = p[:, j]
            all_rows.append(rows)
        print(f"ENDPOINT COMPLETE {endpoint}", flush=True)
    predictions = pd.concat(all_rows, ignore_index=True)
    if predictions.duplicated(["endpoint", "patient_id", "variant", "state"]).any():
        raise ValueError("Duplicate out-of-fold predictions")
    predictions.to_csv(out / "internal_predictions.csv", index=False)
    pd.DataFrame(parity).to_csv(out / "reproduction_parity.csv", index=False)
    summarize_internal(predictions, args.bootstrap, out)
    summarize_sites(sites, external, args.bootstrap, out)
    # Recheck frozen source bytes at completion.
    if load_inputs(args.root.resolve())[-1] != hashes:
        raise ValueError("Inputs changed during analysis")
    files = [p for p in out.rglob("*") if p.is_file() and p.name != "completion.json"]
    write_json(out / "completion.json", {"status": "COMPLETE", "protocol_sha256": protocol_hash,
        "all_50_reproduction_checks_pass": all(r["pass"] for r in parity),
        "internal_prediction_rows": len(predictions), "original_input_hashes_unchanged": True,
        "outputs": {p.relative_to(out).as_posix(): digest(p) for p in files}})
    print("COMPLETE: frozen inputs unchanged; 50 released-prediction parity checks passed.", flush=True)


def summarize_internal(predictions, reps, out):
    points, intervals, contrasts, diagnostics = [], [], [], []
    for endpoint in ENDPOINTS:
        table = predictions[predictions.endpoint.eq(endpoint)]
        reference = table[table.variant.eq("prevalence")].sort_values("patient_id")
        ids, y = reference.patient_id.to_numpy(), reference.label.to_numpy(int)
        k = len(np.unique(y))
        rng = np.random.default_rng(SEED + sum(map(ord, endpoint)))
        indices = rng.integers(0, len(y), (reps, len(y)))
        weights = np.vstack([np.ones(len(y)), np.stack([np.bincount(i, minlength=len(y)) for i in indices])]).astype(float)
        distributions = {}
        for (variant, state), d in table.groupby(["variant", "state"], sort=False):
            d = d.sort_values("patient_id")
            if not np.array_equal(d.patient_id, ids) or not np.array_equal(d.label, y):
                raise ValueError("Patient pairing failed")
            p = d[[f"p_{c}" for c in range(k)]].to_numpy(float)
            dist = metric_distributions(y, p, weights)
            distributions[(variant, state)] = dist
            points.append({"endpoint": endpoint, "variant": variant, "state": state, "n": len(y),
                           **{m: float(v[0]) for m, v in dist.items()}})
            for m, values in dist.items():
                valid = values[1:][np.isfinite(values[1:])]
                lo, hi = np.percentile(valid, [2.5, 97.5])
                intervals.append({"endpoint": endpoint, "variant": variant, "state": state,
                    "metric": m, "estimate": values[0], "ci_lower_95": lo, "ci_upper_95": hi,
                    "valid_replicates": len(valid)})
            for row in calibration_diagnostics(y, p):
                diagnostics.append({"endpoint": endpoint, "variant": variant, "state": state, **row})
        pairs = []
        for variant in ("std_balanced", "std_unweighted"):
            pairs.extend([
                ("matched_sigmoid_effect", (variant, "ensemble_raw"), (variant, "ensemble_sigmoid")),
                ("ensemble_training_effect", (variant, "full_raw"), (variant, "ensemble_raw")),
                ("full_procedure_effect", (variant, "full_raw"), (variant, "ensemble_sigmoid")),
            ])
        for state in ("full_raw", "ensemble_raw", "ensemble_sigmoid"):
            pairs.append(("unweighted_minus_balanced", ("std_balanced", state), ("std_unweighted", state)))
        for a, b in [("none_balanced", "l2_balanced"), ("std_balanced", "l2_std_balanced"),
                     ("none_balanced", "std_balanced"), ("l2_balanced", "l2_std_balanced")]:
            pairs.append(("controlled_preprocessing", (a, "full_raw"), (b, "full_raw")))
        for key in distributions:
            if key[0] != "prevalence":
                pairs.append(("model_minus_prevalence", ("prevalence", "development_only"), key))
        for name, a, b in pairs:
            for metric in distributions[a]:
                delta = distributions[b][metric] - distributions[a][metric]
                valid = delta[1:][np.isfinite(delta[1:])]
                lo, hi = np.percentile(valid, [2.5, 97.5])
                contrasts.append({"endpoint": endpoint, "comparison": name, "first_variant": a[0],
                    "first_state": a[1], "second_variant": b[0], "second_state": b[1], "metric": metric,
                    "delta_second_minus_first": delta[0], "ci_lower_95": lo, "ci_upper_95": hi,
                    "valid_replicates": len(valid)})
    pd.DataFrame(points).to_csv(out / "internal_metrics.csv", index=False)
    pd.DataFrame(intervals).to_csv(out / "internal_metric_intervals.csv", index=False)
    pd.DataFrame(contrasts).to_csv(out / "paired_contrasts.csv", index=False)
    pd.DataFrame(diagnostics).to_csv(out / "calibration_diagnostics.csv", index=False)


def summarize_sites(sites, external, reps, out):
    rows, predictions, contrasts = [], [], []
    for (endpoint, site), table in sites.groupby(["endpoint", "holdout_site"]):
        if set(table.partition) != {"train", "validation", "test"}:
            raise ValueError("Unexpected site partition")
        dev, test = table[table.partition.ne("test")], table[table.partition.eq("test")].sort_values("patient_id")
        if set(dev.patient_id) & set(test.patient_id) or not test.tss.eq(site).all() or dev.tss.eq(site).any():
            raise ValueError("Site separation failed")
        y, classes = test.label.to_numpy(int), np.sort(table.label.unique())
        old = external[external.endpoint.eq(endpoint) & external.holdout_site.eq(site)].sort_values("patient_id")
        if not np.array_equal(old.patient_id, test.patient_id) or not np.array_equal(old.label, y):
            raise ValueError("Site prediction pairing failed")
        baseline = prevalence_predictions(dev.label.to_numpy(int), len(test), classes)
        p = old[[f"prob_calibrated_class_{c}" for c in classes]].to_numpy(float)
        rng = np.random.default_rng(SEED + sum(map(ord, endpoint + site)))
        draws = rng.integers(0, len(y), (reps, len(y)))
        weights = np.vstack([np.ones(len(y)), np.stack([np.bincount(i, minlength=len(y)) for i in draws])]).astype(float)
        states = {"development_prevalence": baseline, "released_calibrated": p}
        dist = {}
        for state, probs in states.items():
            validate_probability(probs, len(test), len(classes))
            dist[state] = metric_distributions(y, probs, weights)
            rows.append({"endpoint": endpoint, "holdout_site": site, "state": state,
                         "test_n": len(test), "development_n": len(dev),
                         **{m: float(v[0]) for m, v in dist[state].items()}})
            d = test.copy(); d["state"] = state
            for j, c in enumerate(classes): d[f"p_{c}"] = probs[:, j]
            predictions.append(d)
        for metric in dist["development_prevalence"]:
            delta = dist["released_calibrated"][metric] - dist["development_prevalence"][metric]
            valid = delta[1:][np.isfinite(delta[1:])]
            lo, hi = np.percentile(valid, [2.5, 97.5])
            contrasts.append({"endpoint": endpoint, "holdout_site": site, "metric": metric,
                "delta_calibrated_minus_development_prevalence": delta[0], "ci_lower_95": lo,
                "ci_upper_95": hi, "valid_replicates": len(valid)})
    pd.DataFrame(rows).to_csv(out / "site_prevalence_metrics.csv", index=False)
    pd.DataFrame(contrasts).to_csv(out / "site_prevalence_contrasts.csv", index=False)
    pd.concat(predictions, ignore_index=True).to_csv(out / "site_prevalence_predictions.csv", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    package = Path(__file__).resolve().parent.parent
    parser.add_argument("--root", type=Path, default=package / "frozen_inputs")
    parser.add_argument("--out", type=Path, default=package / "results")
    parser.add_argument("--bootstrap", type=int, default=1000)
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error("Use at least 100 bootstrap replicates; 1000 is the study default")
    with threadpool_limits(limits=1):
        run(args)
