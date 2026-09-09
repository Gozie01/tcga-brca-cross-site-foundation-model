from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "evidence" / "Stage12_Results_Synthesis_v1" / "tables"
STATS = ROOT / "evidence" / "Stage12_inputs_20260905_030608" / "stage11_statistical_analysis_ablations_v1"
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8,
    "axes.titlesize": 9,
    "axes.labelsize": 8,
    "figure.dpi": 160,
    "savefig.dpi": 400,
    "axes.spines.top": False,
    "axes.spines.right": False,
})

ENDPOINTS = ["ER", "PR", "HER2", "PIK3CA", "PAM50"]
COLORS = {"Raw": "#6b7280", "Calibrated": "#0072B2"}


def parse_ci(value):
    text = str(value).replace("–", "-")
    estimate = float(text.split()[0])
    if "(" not in text:
        return estimate, np.nan, np.nan
    lo, hi = text.split("(", 1)[1].rstrip(")").split("-")
    return estimate, float(lo), float(hi)


def save(fig, name):
    fig.savefig(OUT / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(OUT / f"{name}.png", bbox_inches="tight")
    plt.close(fig)


def figure1():
    fig, ax = plt.subplots(figsize=(7.4, 4.6))
    ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
    navy, blue, green, orange = "#173F5F", "#DCEAF6", "#E5F3EC", "#FFF0D9"

    def stage(y, text, color):
        ax.text(0.04, y, text, ha="left", va="center", fontsize=6.3,
                fontweight="bold", color=color)
        ax.plot([0.04, 0.96], [y-0.026, y-0.026], color=color, lw=0.8, alpha=0.55)

    def box(x, y, w, h, title, body, face):
        ax.add_patch(FancyBboxPatch((x, y), w, h,
            boxstyle="round,pad=0.006,rounding_size=0.010",
            linewidth=0.9, edgecolor=navy, facecolor=face, zorder=2))
        ax.text(x+w/2, y+h*0.67, title, ha="center", va="center", fontsize=6.1,
                fontweight="bold", color=navy, zorder=3)
        ax.text(x+w/2, y+h*0.34, body, ha="center", va="center", fontsize=5.4,
                color="#263B4A", linespacing=1.2, zorder=3)
        return x, y, w, h

    def arrow(a, b, color=navy):
        ax.add_patch(FancyArrowPatch(a, b, arrowstyle="-|>", mutation_scale=8,
            linewidth=0.9, color=color, zorder=4, connectionstyle="arc3,rad=0"))

    def elbow(points, color=navy):
        for a, b in zip(points[:-2], points[1:-1]):
            ax.plot([a[0], b[0]], [a[1], b[1]], color=color, lw=0.9, zorder=4)
        arrow(points[-2], points[-1], color=color)

    ax.text(0.5, 0.96, "Leakage-controlled pathology foundation-model framework",
            ha="center", fontsize=10.3, fontweight="bold", color=navy)
    ax.text(0.5, 0.915, "Frozen representations | patient-level separation | calibration | site-aware evaluation",
            ha="center", fontsize=6.6, color="#526777")

    stage(0.865, "A  DATA AND REPRESENTATION", navy)
    top = [
        box(0.055, 0.675, 0.205, 0.145, "TCGA-BRCA", "347 diagnostic WSIs\n345 patients", blue),
        box(0.285, 0.675, 0.205, 0.145, "Locked patch QC", "345,937 patches | 0.50 μm/px\nmax 1,000/slide; no backfill", blue),
        box(0.515, 0.675, 0.205, 0.145, "Frozen UNI", "224 × 224 input\n1,024-D patch embeddings", blue),
        box(0.745, 0.675, 0.205, 0.145, "Patient representation", "patch mean → slide mean\n→ patient mean", blue),
    ]
    for left, right in zip(top[:-1], top[1:]):
        arrow((left[0]+left[2]+0.004, 0.748), (right[0]-0.004, 0.748))

    stage(0.61, "B  MODELING AND CALIBRATION", "#2E6B57")
    mid = [
        box(0.105, 0.415, 0.24, 0.145, "Five endpoint models", "ER | PR | HER2 | PIK3CA | PAM50\npatient-specific eligibility", green),
        box(0.38, 0.415, 0.24, 0.145, "Linear probes", "StandardScaler + L2 logistic model\ninner-fold regularization selection", green),
        box(0.655, 0.415, 0.24, 0.145, "Development-only calibration", "Platt sigmoid: binary endpoints\ntemperature scaling: PAM50", green),
    ]
    arrow((0.848, 0.675), (0.34, 0.565))
    arrow((0.345+0.004, 0.488), (0.38-0.004, 0.488))
    arrow((0.62+0.004, 0.488), (0.655-0.004, 0.488))

    stage(0.35, "C  VALIDATION AND STATISTICAL INFERENCE", "#9A6500")
    bottom = [
        box(0.075, 0.13, 0.25, 0.155, "Internal validation", "5 patient-separated outer folds\none OOF prediction per patient", orange),
        box(0.375, 0.13, 0.25, 0.155, "Institution holdouts", "A2 | C8 | E2 | BH\nwithin TCGA-BRCA", orange),
        box(0.675, 0.13, 0.25, 0.155, "Statistical evaluation", "AUROC | AUPRC | Brier | ECE\n1,000 patient bootstraps", orange),
    ]
    arrow((0.775, 0.415), (0.30, 0.29))
    arrow((0.775, 0.415), (0.60, 0.29))
    arrow((0.325+0.004, 0.208), (0.675-0.004, 0.208), color="#9A6500")
    arrow((0.625+0.004, 0.208), (0.675-0.004, 0.208), color="#9A6500")

    ax.text(0.5, 0.055, "No patient, slide, preprocessing, or calibration leakage across test partitions",
            ha="center", fontsize=6.1, color="#6B7280", style="italic")
    save(fig, "Figure1_study_design")


def figure2():
    df = pd.read_csv(DATA / "Table2_internal_performance.csv")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.3), sharey=True)
    ybase = np.arange(len(ENDPOINTS))[::-1]
    for ax, metric, title in zip(axes, ["auroc", "auprc"], ["AUROC", "AUPRC"]):
        for state, offset in [("Raw", -0.10), ("Calibrated", 0.10)]:
            rows = df[df.calibration.eq(state)].set_index("endpoint").loc[ENDPOINTS]
            vals = [parse_ci(v) for v in rows[metric]]
            est = np.array([v[0] for v in vals]); lo = np.array([v[1] for v in vals]); hi = np.array([v[2] for v in vals])
            ax.errorbar(est, ybase + offset, xerr=[est-lo, hi-est], fmt="o", ms=4.5, capsize=2,
                        lw=1.2, color=COLORS[state], label=state)
        ax.axvline(0.5, color="#bbbbbb", ls="--", lw=0.8)
        ax.set_xlim(0.25 if metric == "auprc" else 0.45, 1.0)
        ax.set_xlabel(title)
        ax.grid(axis="x", color="#e5e7eb", lw=0.6)
    axes[0].set_yticks(ybase, ENDPOINTS)
    axes[1].legend(frameon=False, loc="lower right")
    fig.suptitle("Internal patient-level out-of-fold performance (95% bootstrap CI)", y=1.01, fontweight="bold")
    fig.tight_layout()
    save(fig, "Figure2_internal_performance")


def figure3():
    df = pd.read_csv(STATS / "stage11_paired_calibration_statistics_v1.csv", encoding="utf-8-sig")
    df["endpoint"] = df.endpoint.str.replace("_binary", "", regex=False).str.replace("_primary_4class", "", regex=False)
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)
    y = np.arange(len(ENDPOINTS))[::-1]
    for ax, metric, title in zip(axes, ["brier", "ece_10bin"], ["Δ Brier score", "Δ ECE (10 bins)"]):
        rows = df[df.metric.eq(metric)].set_index("endpoint").reindex(ENDPOINTS)
        est = rows.delta_second_minus_first.to_numpy(float)
        lo = rows.delta_ci_lower_95.to_numpy(float); hi = rows.delta_ci_upper_95.to_numpy(float)
        sig = rows.paired_bootstrap_p_two_sided.to_numpy(float) < 0.05
        for i in range(len(ENDPOINTS)):
            ax.errorbar(est[i], y[i], xerr=[[est[i]-lo[i]], [hi[i]-est[i]]], fmt="o", capsize=2,
                        color="#0072B2" if sig[i] else "#6b7280", ms=5)
        ax.axvline(0, color="black", lw=0.8)
        ax.set_xlabel(f"{title}\n(calibrated minus raw; negative is better)")
        ax.grid(axis="x", color="#e5e7eb", lw=0.6)
    axes[0].set_yticks(y, ENDPOINTS)
    fig.suptitle("Paired patient-level calibration effects (95% bootstrap CI)", y=1.02, fontweight="bold")
    fig.tight_layout()
    save(fig, "Figure3_calibration_effects")


def figure4():
    df = pd.read_csv(DATA / "Table3_site_holdout_performance.csv")
    df = df[df.calibration_state.eq("Calibrated")]
    sites = ["A2", "C8", "E2", "BH"]
    metrics = [("auroc", "AUROC", 0.5, 1.0), ("auprc", "AUPRC", 0.25, 1.0), ("ece_10bin", "ECE", 0, 0.4)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 4.0), constrained_layout=True)
    for ax, (metric, title, vmin, vmax) in zip(axes, metrics):
        pivot = df.pivot(index="endpoint", columns="holdout_site", values=metric).reindex(index=ENDPOINTS, columns=sites)
        im = ax.imshow(pivot.values, cmap="Blues_r" if metric == "ece_10bin" else "Blues", vmin=vmin, vmax=vmax, aspect="auto")
        for i in range(len(ENDPOINTS)):
            for j in range(len(sites)):
                sub = df[(df.endpoint == ENDPOINTS[i]) & (df.holdout_site == sites[j])].iloc[0]
                ax.text(j, i, f"{pivot.iloc[i,j]:.2f}\n(n={int(sub.test_n)})", ha="center", va="center",
                        fontsize=6.4, color="white" if (metric != "ece_10bin" and pivot.iloc[i,j] > .78) else "black")
        ax.set_xticks(range(len(sites)), sites)
        ax.set_yticks(range(len(ENDPOINTS)), ENDPOINTS if ax is axes[0] else [])
        ax.set_title(title, fontweight="bold")
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.03)
    fig.suptitle("Calibrated institution-held-out performance within TCGA-BRCA", y=1.04, fontweight="bold")
    save(fig, "Figure4_site_holdout")


def figure5():
    df = pd.read_csv(DATA / "Table4_representation_ablation_auroc.csv")
    order = ["full_1024", "l2_normalized", "pca_32", "pca_64", "pca_128"]
    labels = ["Full 1024", "L2 norm", "PCA-32", "PCA-64", "PCA-128"]
    pivot = df.pivot(index="endpoint", columns="ablation", values="mean_auroc").reindex(index=ENDPOINTS, columns=order)
    fig, ax = plt.subplots(figsize=(7.2, 3.4))
    x = np.arange(len(ENDPOINTS)); width = 0.16
    palette = ["#4b5563", "#0072B2", "#56B4E9", "#009E73", "#E69F00"]
    for j, (key, label, color) in enumerate(zip(order, labels, palette)):
        ax.bar(x + (j-2)*width, pivot[key], width=width, label=label, color=color)
    ax.set_xticks(x, ENDPOINTS)
    ax.set_ylim(0.48, 0.91)
    ax.set_ylabel("Mean outer-fold AUROC")
    ax.grid(axis="y", color="#e5e7eb", lw=0.6)
    ax.legend(ncol=5, fontsize=6.5, frameon=False, loc="upper center")
    ax.text(3-0.16, pivot.loc["PIK3CA", "l2_normalized"]+0.012, "*", ha="center", fontsize=13)
    for j in [2, 3, 4]:
        ax.text(4 + (j-2)*width, pivot.iloc[4, j]+0.010, "*", ha="center", fontsize=13)
    ax.set_title("Representation-ablation discrimination", fontweight="bold")
    ax.text(0.99, 0.02, "* paired bootstrap p<0.05 vs full representation", transform=ax.transAxes,
            ha="right", va="bottom", fontsize=6.5)
    fig.tight_layout()
    save(fig, "Figure5_ablation")


if __name__ == "__main__":
    figure1(); figure2(); figure3(); figure4(); figure5()
    print(f"Generated figures in {OUT}")
