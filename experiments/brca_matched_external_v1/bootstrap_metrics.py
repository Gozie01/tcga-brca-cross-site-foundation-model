import numpy as np
from scipy.special import expit,logit
from scipy.optimize import brentq,minimize
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

