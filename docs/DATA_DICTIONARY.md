# Data Dictionary

`cohort_manifest.csv`: `patient_id` is the public TCGA case identifier; `feature_row` indexes the patient feature matrix; `tss` is the TCGA tissue-source site; endpoint and evaluability columns define task labels.

`slide_patient_mapping.csv`: links a non-identifying slide filename to `patient_id`, `tss`, pooling counts, and feature rows.

`frozen_patient_folds.csv`: endpoint-specific `patient_id`, integer `label`, `label_name`, and `internal_fold` (0-4).

`institution_holdout_assignments.csv`: endpoint-specific patient assignments with `holdout_site` and `partition` (`train`, `validation`, or `test`).

Prediction files share `endpoint`, `patient_id`, `label`, site/fold fields, and `prob_raw_class_*` or `prob_calibrated_class_*` probabilities. Metric files report AUROC, AUPRC, Brier score, log loss, ECE, and associated sample counts. Feature-index rows map `feature_row` to the matching feature-matrix row.
