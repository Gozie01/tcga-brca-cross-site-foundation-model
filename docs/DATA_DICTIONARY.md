# Data dictionary

## Identifiers and partitions

| Field | Meaning |
|---|---|
| `patient_id` | Public TCGA case barcode used as the patient-level modeling identifier. |
| `file_name` | TCGA diagnostic-slide filename. Raw slide content is not included. |
| `feature_row` | Zero-based row in `stage10a_uni_patient_features_v1.npy`. |
| `source` | Upstream embedding source (`batch01` or `multisite`). |
| `tss` | TCGA tissue source site code derived from the case barcode. |
| `endpoint` | Molecular phenotype modeling task. |
| `internal_fold` | Frozen patient-level outer fold, numbered 0–4. |
| `holdout_site` | TCGA tissue source site used as the test domain. |
| `partition` | `train`, `validation`, or `test` for a site-holdout experiment. |

## Labels

| Endpoint | Encoded classes |
|---|---|
| `ER_binary` | 0 = negative; 1 = positive. |
| `PR_binary` | 0 = negative; 1 = positive. |
| `HER2_binary` | 0 = negative; 1 = positive. |
| `PIK3CA_binary` | 0 = wild type; 1 = mutated. |
| `PAM50_primary_4class` | 0 = Luminal A; 1 = Luminal B; 2 = Basal-like; 3 = HER2-enriched. |

## Predictions and metrics

| Field pattern | Meaning |
|---|---|
| `prob_raw_class_<k>` | Uncalibrated predicted probability for encoded class `k`. |
| `prob_calibrated_class_<k>` | Validation-calibrated probability for class `k`. |
| `calibration` | Raw or calibrated evaluation state/method. |
| `auroc` | Patient-level area under the receiver operating characteristic curve; macro one-vs-rest for PAM50. |
| `auprc` | Patient-level average precision; macro one-vs-rest for PAM50. |
| `balanced_accuracy` | Mean recall across classes. |
| `macro_f1` | Unweighted mean class-wise F1 score. |
| `mcc` | Matthews correlation coefficient. |
| `brier` | Brier score; multiclass generalization for PAM50. |
| `log_loss` | Cross-entropy/logarithmic loss. |
| `ece_10bin` | Expected calibration error computed with 10 bins. |
| `ci_lower_95`, `ci_upper_95` | Percentile 95% bootstrap confidence limits. |
| `paired_bootstrap_p_two_sided` | Finite-sample-corrected two-sided paired-bootstrap p-value. |

## Feature matrix

`stage10a_uni_patient_features_v1.npy` is a `(345, 1024)` float32 NumPy array. Each row is the hierarchical mean of UNI patch representations: patch mean within slide, followed by slide mean within patient. Row identity is defined exclusively by `stage10a_uni_patient_features_index_v1.csv`.
