# Reproducibility

All stochastic operations use seed `20260905`; parameters are frozen in `configs/experiment_config.yaml`. Stage 9A applies the policy in `configs/patch_qc_policy.yaml`. Stage 9B extracts 1024-dimensional UNI features. Stage 10A pools patch features to slides then patients and creates endpoint-specific five-fold patient splits. Stage 10B fits L2 logistic-regression probes and performs fold-contained sigmoid calibration. Stage 10C uses Platt calibration for binary tasks and temperature scaling for PAM50 on frozen validation sets. Stage 11 performs 1,000-replicate paired bootstrap intervals, PCA, and regularization ablations.

To reconstruct WSIs, use GDC project `TCGA-BRCA`, diagnostic-slide image files, and the GDC manifest records retained in the upstream audit package. Do not commit downloaded files. The supplied frozen CSVs allow reproduction of manuscript tables and figures without redistributing images.

Run `python figures/generate_manuscript_figures_and_tables.py` after arranging Stage 10-11 outputs under `stage12_work/Stage12_inputs_20260905_030608/`. Table 1 uses endpoint eligibility, Table 2 internal bootstrap metrics, Table 3 held-out metrics, Table 4 ablations, and Table 5 paired comparisons.
