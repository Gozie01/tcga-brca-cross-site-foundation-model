# Frozen UNI external validation results

Computed from the supplied TCGA checkpoints and CPTAC batch-32 patient features. No training or external recalibration was performed. All 20 frozen models passed local replay against stored TCGA predictions (maximum difference 3.89e-16).

The existing four site-holdout models are reported separately. The ranges below describe model-to-model variation, not confidence intervals. Per-model 95% patient-bootstrap intervals are in external_metrics_95CI.csv. No model was chosen using external performance.

| Endpoint | Evaluable patients | Calibrated AUROC range | Calibrated Brier range |
|---|---:|---:|---:|
| ER_binary | 118 | 0.856–0.908 | 0.182–0.209 |
| HER2_binary | 95 | 0.436–0.614 | 0.120–0.234 |
| PAM50_primary_4class | 114 | 0.815–0.854 | 0.451–0.578 |
| PIK3CA_binary | 120 | 0.565–0.733 | 0.202–0.232 |
| PR_binary | 112 | 0.690–0.824 | 0.211–0.317 |

ER and PAM50 have the most consistent discrimination across these models. HER2 discrimination is weak and imprecise (13 positives among 95 evaluable patients). PR and PIK3CA vary substantially across source models. Calibration transfer is not uniformly reliable; see both raw and calibrated metrics rather than assuming calibration improves external probabilities. PAM50 Brier is an unnormalised four-class sum and is not directly comparable to binary Brier.

These results remain provisional for manuscript interpretation because all 391 slides have pending human tumor review. Technical QC and extraction-coordinate correspondence are verified, but neither establishes tumor-region content or matching to the molecularly assayed block. Pretraining independence is also unestablished.

The immediate next step is qualified histology review, blinded to outcomes and predictions, with documented decisions. If the eligible cohort changes, rerun the same frozen evaluation and retain this original version. Representation comparison remains a separate analysis requiring matching TCGA-side models/features for each encoder.
