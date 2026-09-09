# Proposed manuscript framing

**Working title:** Evaluating Cross-Site Generalization of Breast Cancer Molecular Phenotype Prediction from Whole-Slide Histology Using UNI Foundation-Model Embeddings

**Central claim:** Frozen UNI features support strong patient-level prediction of ER and PAM50 and moderate PR prediction, while HER2 and PIK3CA remain substantially more difficult; probability calibration and cross-site evaluation reveal endpoint-specific reliability limitations that discrimination alone would obscure.

# Results

## Cohort and evaluation integrity

The locked cohort comprised 345 patients represented by 347 diagnostic whole-slide images and 345,937 quality-controlled patches. Endpoint-specific evaluable cohorts included 337 patients for ER, 337 for PR, 314 for HER2, 345 for PIK3CA, and 335 for the four-class PAM50 task. All results were evaluated at the patient level using frozen five-fold out-of-fold partitions and four institution-held-out experiments (A2, C8, E2, and BH), with no patient overlap between partitions.

## Internal predictive performance

- ER: calibrated AUROC 0.843 (95% CI 0.782–0.900) and AUPRC 0.923 (95% CI 0.882–0.962).
- PR: calibrated AUROC 0.768 (95% CI 0.709–0.825) and AUPRC 0.830 (95% CI 0.774–0.894).
- HER2: calibrated AUROC 0.672 (95% CI 0.597–0.751) and AUPRC 0.353 (95% CI 0.260–0.499).
- PIK3CA: calibrated AUROC 0.589 (95% CI 0.530–0.654) and AUPRC 0.447 (95% CI 0.369–0.550).
- PAM50: calibrated AUROC 0.856 (95% CI 0.827–0.882) and AUPRC 0.633 (95% CI 0.581–0.695).

ER and PAM50 showed the strongest discrimination. PR was moderately predictable, whereas HER2 and especially PIK3CA were more difficult. The calibrated PIK3CA AUROC remained modest, indicating that calibration improved probability reliability without converting the endpoint into a strongly discriminative task.

## Effect of calibration

- ER: Brier-score change -0.019 (95% CI -0.035 to -0.004, p=0.022); ECE change -0.047 (95% CI -0.089 to -0.003, p=0.036).
- PR: Brier-score change -0.039 (95% CI -0.061 to -0.018, p<0.002); ECE change -0.042 (95% CI -0.103 to +0.006, p=0.074).
- HER2: Brier-score change -0.041 (95% CI -0.061 to -0.020, p<0.002); ECE change -0.080 (95% CI -0.111 to -0.026, p=0.006).
- PIK3CA: Brier-score change -0.089 (95% CI -0.120 to -0.058, p<0.002); ECE change -0.186 (95% CI -0.220 to -0.114, p<0.002).
- PAM50: Brier-score change -0.009 (95% CI -0.037 to +0.023, p=0.579); ECE change +0.051 (95% CI -0.014 to +0.101, p=0.120).

Calibration significantly improved probability accuracy for ER, PR, HER2, and PIK3CA by Brier score. Discrimination changes were not statistically supported. PAM50 temperature scaling did not provide a supported calibration benefit and tended to increase ECE, so raw PAM50 probabilities should remain the primary presentation with calibrated results reported as sensitivity analysis.

## Institution-held-out generalization

- ER: mean held-out-site AUROC 0.869; range 0.814–0.931; worst-site AUPRC 0.833; mean ECE 0.093.
- PR: mean held-out-site AUROC 0.781; range 0.730–0.819; worst-site AUPRC 0.768; mean ECE 0.125.
- HER2: mean held-out-site AUROC 0.656; range 0.577–0.720; worst-site AUPRC 0.297; mean ECE 0.199.
- PIK3CA: mean held-out-site AUROC 0.613; range 0.587–0.661; worst-site AUPRC 0.404; mean ECE 0.101.
- PAM50: mean held-out-site AUROC 0.809; range 0.767–0.827; worst-site AUPRC 0.568; mean ECE 0.123.

Cross-site robustness was endpoint-dependent. ER retained strong discrimination at every held-out site, and PAM50 remained comparatively stable. HER2 showed the greatest site-related instability, including poor calibrated performance at C8, while PIK3CA remained weak across sites. These findings motivate reporting both mean and worst-site performance rather than pooled performance alone.

## Representation ablations

L2 normalization preserved PAM50 performance and improved the patient-level point estimates for ER, PR, HER2, and PIK3CA relative to the original full-dimensional standardized probe. The paired AUROC improvement was statistically supported only for PIK3CA (+0.033; 95% CI +0.004 to +0.064; p=0.020). PCA compression had endpoint-specific effects: PCA-32 improved PIK3CA probability error but significantly reduced PAM50 AUROC and AUPRC. Across PCA-32, PCA-64, and PCA-128, PAM50 discrimination was significantly below the full 1024-dimensional representation, indicating that molecular-subtype information was distributed across the higher-dimensional UNI feature space.

## Interpretation boundary

The institution-held-out experiments measure cross-site generalization within TCGA-BRCA and should not be described as independent external-cohort validation. The results support morphological association and patient-level prediction, not clinical replacement of molecular assays. Independent multi-cohort validation and prospective assessment remain necessary.
