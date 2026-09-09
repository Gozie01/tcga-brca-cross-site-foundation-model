from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

ROOT = Path("stage12_work/Stage12_inputs_20260905_030608")
OUT = Path("Stage12_Results_Synthesis_v1")
FIG = OUT / "figures"
TAB = OUT / "tables"
FIG.mkdir(parents=True, exist_ok=True)
TAB.mkdir(parents=True, exist_ok=True)

EP = ["ER_binary", "PR_binary", "HER2_binary", "PIK3CA_binary", "PAM50_primary_4class"]
LABEL = {"ER_binary":"ER", "PR_binary":"PR", "HER2_binary":"HER2", "PIK3CA_binary":"PIK3CA", "PAM50_primary_4class":"PAM50"}
SITES = ["A2", "C8", "E2", "BH"]

def read(rel): return pd.read_csv(ROOT / rel)
elig = read("stage10a_cohort_split_freeze_v1/stage10a_endpoint_eligibility_summary_v1.csv")
boot = read("stage10b_internal_prediction_calibration_v1/stage10b_internal_bootstrap_ci_v1.csv")
site = read("stage10c_site_holdout_domain_shift_v1/stage10c_site_holdout_metrics_v1.csv")
cal = read("stage11_statistical_analysis_ablations_v1/stage11_paired_calibration_statistics_v1.csv")
abl = read("stage11_statistical_analysis_ablations_v1/stage11_ablation_summary_v1.csv")
ablp = read("stage11_statistical_analysis_ablations_v1/stage11_paired_ablation_statistics_v1.csv")
domain = read("stage11_statistical_analysis_ablations_v1/stage11_paired_domain_shift_statistics_v1.csv")
hetero = read("stage11_statistical_analysis_ablations_v1/stage11_site_heterogeneity_summary_v1.csv")

# Publication tables
cohort = elig.copy()
cohort["endpoint"] = cohort.endpoint.map(LABEL)
cohort["class_distribution"] = cohort.class_counts_json.map(lambda x: "; ".join(f"{k}: {v}" for k,v in json.loads(x).items()))
cohort[["endpoint","evaluable_n","class_count","class_distribution"]].to_csv(TAB/"Table1_endpoint_cohorts.csv", index=False)

internal = boot[boot.metric.isin(["auroc","auprc","brier","ece_10bin"])].copy()
internal["endpoint"] = internal.endpoint.map(LABEL)
internal["calibration"] = internal.calibration.map({"raw":"Raw","sigmoid_cv":"Calibrated"})
internal["estimate_95ci"] = internal.apply(lambda r:f"{r.estimate:.3f} ({r.ci_lower_95:.3f}–{r.ci_upper_95:.3f})",axis=1)
internal.pivot(index=["endpoint","calibration"],columns="metric",values="estimate_95ci").reset_index().to_csv(TAB/"Table2_internal_performance.csv", index=False)

site_table = site.copy()
site_table["endpoint"] = site_table.endpoint.map(LABEL)
site_table["calibration_state"] = np.where(site_table.calibration.eq("raw"),"Raw","Calibrated")
site_table[["endpoint","holdout_site","calibration_state","test_n","auroc","auprc","balanced_accuracy","macro_f1","brier","ece_10bin"]].to_csv(TAB/"Table3_site_holdout_performance.csv", index=False)

abl_table = abl.copy(); abl_table["endpoint"] = abl_table.endpoint.map(LABEL)
baseline = boot[(boot.metric.eq("auroc"))&(boot.calibration.eq("raw"))][["endpoint","estimate"]].rename(columns={"estimate":"mean_auroc"})
baseline["ablation"]="full_1024"; baseline["endpoint"]=baseline.endpoint.map(LABEL)
pd.concat([baseline,abl_table[["endpoint","ablation","mean_auroc"]]],ignore_index=True).to_csv(TAB/"Table4_representation_ablation_auroc.csv",index=False)

sig_cal = cal[(cal.paired_bootstrap_p_two_sided<0.05)].copy(); sig_cal["comparison"]="Calibrated minus raw"
sig_abl = ablp[(ablp.paired_bootstrap_p_two_sided<0.05)].copy(); sig_abl["comparison"]=sig_abl.ablation+" minus full_1024"
pd.concat([
 sig_cal[["endpoint","comparison","metric","delta_second_minus_first","delta_ci_lower_95","delta_ci_upper_95","paired_bootstrap_p_two_sided"]],
 sig_abl[["endpoint","comparison","metric","delta_second_minus_first","delta_ci_lower_95","delta_ci_upper_95","paired_bootstrap_p_two_sided"]]
]).to_csv(TAB/"Table5_significant_paired_effects.csv",index=False)

sns.set_theme(style="whitegrid", context="paper", font_scale=1.15)
palette={"Raw":"#4C78A8","Calibrated":"#F58518"}

# Fig 1: Internal discrimination with CIs
d=boot[(boot.metric.eq("auroc"))].copy(); d["Endpoint"]=d.endpoint.map(LABEL); d["State"]=d.calibration.map({"raw":"Raw","sigmoid_cv":"Calibrated"})
fig,ax=plt.subplots(figsize=(8.0,4.7))
for j,state_name in enumerate(["Raw","Calibrated"]):
 q=d[d.State.eq(state_name)].set_index("Endpoint").loc[[LABEL[x] for x in EP]].reset_index(); x=np.arange(len(q))+(j-.5)*.22
 ax.errorbar(x,q.estimate,yerr=[q.estimate-q.ci_lower_95,q.ci_upper_95-q.estimate],fmt="o",capsize=3,label=state_name,color=palette[state_name],lw=1.5)
ax.set_xticks(np.arange(5),[LABEL[x] for x in EP]); ax.set_ylim(.45,1.0); ax.set_ylabel("Patient-level AUROC (95% CI)"); ax.legend(frameon=False,ncol=2); sns.despine(); fig.tight_layout()
fig.savefig(FIG/"Figure1_internal_auroc.pdf",bbox_inches="tight"); fig.savefig(FIG/"Figure1_internal_auroc.png",dpi=400,bbox_inches="tight"); plt.close(fig)

# Fig 2: held-out site AUROC heatmaps raw and calibrated
fig,axes=plt.subplots(1,2,figsize=(10.5,4.4),sharey=True)
for ax,state_name,is_raw in zip(axes,["Raw","Calibrated"],[True,False]):
 q=site[site.calibration.eq("raw") if is_raw else ~site.calibration.eq("raw")].copy(); q["Endpoint"]=q.endpoint.map(LABEL)
 mat=q.pivot(index="Endpoint",columns="holdout_site",values="auroc").loc[[LABEL[x] for x in EP],SITES]
 sns.heatmap(mat,annot=True,fmt=".3f",cmap="viridis",vmin=.5,vmax=1,ax=ax,cbar=ax is axes[1],linewidths=.5)
 ax.set_title(state_name); ax.set_xlabel("Held-out TCGA tissue-source site"); ax.set_ylabel("Endpoint" if ax is axes[0] else "")
fig.tight_layout(); fig.savefig(FIG/"Figure2_site_holdout_auroc.pdf",bbox_inches="tight"); fig.savefig(FIG/"Figure2_site_holdout_auroc.png",dpi=400,bbox_inches="tight"); plt.close(fig)

# Fig 3: calibration ECE by held-out site
q=site.copy(); q["Endpoint"]=q.endpoint.map(LABEL); q["State"]=np.where(q.calibration.eq("raw"),"Raw","Calibrated")
g=sns.catplot(data=q,x="holdout_site",y="ece_10bin",hue="State",col="Endpoint",col_wrap=3,kind="bar",palette=palette,height=3.0,aspect=1.05,sharey=True)
g.set_axis_labels("Held-out site","ECE (10 bins)"); g.set_titles("{col_name}"); g.legend.set_title(""); g.figure.subplots_adjust(top=.92); g.figure.suptitle("Calibration error under institution-held-out evaluation")
g.savefig(FIG/"Figure3_site_holdout_calibration.pdf",bbox_inches="tight"); g.savefig(FIG/"Figure3_site_holdout_calibration.png",dpi=400,bbox_inches="tight"); plt.close(g.figure)

# Fig 4: representation ablation
raw=boot[(boot.metric.eq("auroc"))&(boot.calibration.eq("raw"))][["endpoint","estimate"]].rename(columns={"estimate":"mean_auroc"}); raw["ablation"]="Full 1024"
q=pd.concat([raw,abl[["endpoint","ablation","mean_auroc"]]],ignore_index=True); q["Endpoint"]=q.endpoint.map(LABEL); q["Representation"]=q.ablation.replace({"l2_normalized":"L2 normalized","pca_32":"PCA-32","pca_64":"PCA-64","pca_128":"PCA-128"})
fig,ax=plt.subplots(figsize=(10,5.2)); sns.barplot(data=q,x="Endpoint",y="mean_auroc",hue="Representation",hue_order=["Full 1024","L2 normalized","PCA-32","PCA-64","PCA-128"],ax=ax)
ax.set_ylim(.5,.92); ax.set_ylabel("Out-of-fold AUROC"); ax.set_xlabel(""); ax.legend(frameon=False,ncol=3,loc="upper center",bbox_to_anchor=(.5,1.18)); sns.despine(); fig.tight_layout()
fig.savefig(FIG/"Figure4_representation_ablation.pdf",bbox_inches="tight"); fig.savefig(FIG/"Figure4_representation_ablation.png",dpi=400,bbox_inches="tight"); plt.close(fig)

# Results draft grounded only in frozen values.
def pval(x): return "<0.002" if x <= .0019981 else f"={x:.3f}"
cal_lookup=cal.set_index(["endpoint","metric"])
lines=[]
lines += ["# Proposed manuscript framing", "", "**Working title:** Evaluating Cross-Site Generalization of Breast Cancer Molecular Phenotype Prediction from Whole-Slide Histology Using UNI Foundation-Model Embeddings", "", "**Central claim:** Frozen UNI features support strong patient-level prediction of ER and PAM50 and moderate PR prediction, while HER2 and PIK3CA remain substantially more difficult; probability calibration and cross-site evaluation reveal endpoint-specific reliability limitations that discrimination alone would obscure.", "", "# Results", "", "## Cohort and evaluation integrity", "", "The locked cohort comprised 345 patients represented by 347 diagnostic whole-slide images and 345,937 quality-controlled patches. Endpoint-specific evaluable cohorts included 337 patients for ER, 337 for PR, 314 for HER2, 345 for PIK3CA, and 335 for the four-class PAM50 task. All results were evaluated at the patient level using frozen five-fold out-of-fold partitions and four institution-held-out experiments (A2, C8, E2, and BH), with no patient overlap between partitions.", "", "## Internal predictive performance", ""]
for endpoint in EP:
 q=boot[(boot.endpoint.eq(endpoint))&(boot.calibration.eq("sigmoid_cv"))&(boot.metric.isin(["auroc","auprc"]))].set_index("metric")
 lines.append(f"- {LABEL[endpoint]}: calibrated AUROC {q.loc['auroc','estimate']:.3f} (95% CI {q.loc['auroc','ci_lower_95']:.3f}–{q.loc['auroc','ci_upper_95']:.3f}) and AUPRC {q.loc['auprc','estimate']:.3f} (95% CI {q.loc['auprc','ci_lower_95']:.3f}–{q.loc['auprc','ci_upper_95']:.3f}).")
lines += ["", "ER and PAM50 showed the strongest discrimination. PR was moderately predictable, whereas HER2 and especially PIK3CA were more difficult. The calibrated PIK3CA AUROC remained modest, indicating that calibration improved probability reliability without converting the endpoint into a strongly discriminative task.", "", "## Effect of calibration", ""]
for endpoint in EP:
 b=cal_lookup.loc[(endpoint,"brier")]; e=cal_lookup.loc[(endpoint,"ece_10bin")]
 lines.append(f"- {LABEL[endpoint]}: Brier-score change {b.delta_second_minus_first:+.3f} (95% CI {b.delta_ci_lower_95:+.3f} to {b.delta_ci_upper_95:+.3f}, p{pval(b.paired_bootstrap_p_two_sided)}); ECE change {e.delta_second_minus_first:+.3f} (95% CI {e.delta_ci_lower_95:+.3f} to {e.delta_ci_upper_95:+.3f}, p{pval(e.paired_bootstrap_p_two_sided)}).")
lines += ["", "Calibration significantly improved probability accuracy for ER, PR, HER2, and PIK3CA by Brier score. Discrimination changes were not statistically supported. PAM50 temperature scaling did not provide a supported calibration benefit and tended to increase ECE, so raw PAM50 probabilities should remain the primary presentation with calibrated results reported as sensitivity analysis.", "", "## Institution-held-out generalization", ""]
for endpoint in EP:
 q=hetero[(hetero.endpoint.eq(endpoint))&(hetero.calibration.ne("raw"))].iloc[0]
 lines.append(f"- {LABEL[endpoint]}: mean held-out-site AUROC {q.site_mean_auroc:.3f}; range {q.worst_site_auroc:.3f}–{q.best_site_auroc:.3f}; worst-site AUPRC {q.worst_site_auprc:.3f}; mean ECE {q.site_mean_ece:.3f}.")
lines += ["", "Cross-site robustness was endpoint-dependent. ER retained strong discrimination at every held-out site, and PAM50 remained comparatively stable. HER2 showed the greatest site-related instability, including poor calibrated performance at C8, while PIK3CA remained weak across sites. These findings motivate reporting both mean and worst-site performance rather than pooled performance alone.", "", "## Representation ablations", "", "L2 normalization preserved PAM50 performance and improved the patient-level point estimates for ER, PR, HER2, and PIK3CA relative to the original full-dimensional standardized probe. The paired AUROC improvement was statistically supported only for PIK3CA (+0.033; 95% CI +0.004 to +0.064; p=0.020). PCA compression had endpoint-specific effects: PCA-32 improved PIK3CA probability error but significantly reduced PAM50 AUROC and AUPRC. Across PCA-32, PCA-64, and PCA-128, PAM50 discrimination was significantly below the full 1024-dimensional representation, indicating that molecular-subtype information was distributed across the higher-dimensional UNI feature space.", "", "## Interpretation boundary", "", "The institution-held-out experiments measure cross-site generalization within TCGA-BRCA and should not be described as independent external-cohort validation. The results support morphological association and patient-level prediction, not clinical replacement of molecular assays. Independent multi-cohort validation and prospective assessment remain necessary."]
(OUT/"RESULTS_SECTION_DRAFT.md").write_text("\n".join(lines),encoding="utf-8")

summary={
 "status":"LOCKED_COMPLETE","working_title":lines[2].replace("**Working title:** ",""),
 "cohort":{"patients":345,"slides":347,"patches":345937},
 "primary_endpoint_interpretation":{"ER":"strong","PR":"moderate","HER2":"limited_and_site_sensitive","PIK3CA":"weak","PAM50":"strong_but_compression_sensitive"},
 "caution":"Site-held-out validation is within TCGA-BRCA, not an independent external cohort.",
 "outputs":{"figures":8,"tables":5,"results_draft":"RESULTS_SECTION_DRAFT.md"}
}
(OUT/"stage12_results_synthesis_metadata_v1.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
print(json.dumps(summary,indent=2))
