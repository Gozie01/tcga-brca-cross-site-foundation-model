"""Publication figure from completed revision-1 controls, never pilot data."""
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from audit_assets import PACKAGE

df = pd.read_csv(PACKAGE / "data/completed_control_contrasts.csv")
endpoints = ["ER_binary", "PR_binary", "HER2_binary", "PIK3CA_binary", "PAM50_primary_4class"]
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.4), sharey=True)
for ax, variant, title, color in zip(axes, ["std_balanced", "std_unweighted"],
        ["Class-weighted probes", "Unweighted probes"], ["#245D91", "#B36527"]):
    d = df[df.comparison.eq("matched_sigmoid_effect") & df.metric.eq("brier") & df.first_variant.eq(variant)].set_index("endpoint").loc[endpoints]
    value, low, high = [d[name].to_numpy() for name in ["delta_second_minus_first", "ci_lower_95", "ci_upper_95"]]
    ax.errorbar(value, np.arange(5), xerr=np.vstack([value-low, high-value]), fmt="o", color=color,
                ecolor=color, capsize=3, markersize=5, linewidth=1.5)
    ax.axvline(0, color="#555555", linestyle="--", linewidth=1)
    ax.set_title(title, pad=12)
    ax.set_yticks(np.arange(5), ["ER", "PR", "HER2", "PIK3CA", "PAM50"])
    ax.set_xlim(-.077, .039)
    ax.grid(axis="x", color="#E3E6E9", linewidth=.7)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
axes[0].invert_yaxis()
fig.supxlabel("Change in Brier score: calibrated minus raw ensemble", y=.06)
fig.text(.5, .005, "Paired patient bootstrap 95% intervals; negative values indicate improvement. PAM50 uses multiclass Brier score.", ha="center", fontsize=8)
fig.tight_layout(rect=(0, .10, 1, 1))
output = PACKAGE / "manuscript_edits/matched_calibration_effects.png"
fig.savefig(output, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(output)
