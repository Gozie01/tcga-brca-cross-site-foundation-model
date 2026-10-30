# Frozen UNI external evaluation, version 1

This package evaluates all 20 pre-existing TCGA site-holdout models on CPTAC patient features. No model, scaler, calibrator, threshold or ensemble is fitted or selected using CPTAC. The four overlapping site-holdout models are reported separately and are not independent replications. A single final model trained on all TCGA patients was not assumed.

Run `python run_external.py` from any directory after installing requirements.txt in a separate environment. Inputs and models are bundled; no WSI downloads or GPU are required. Model files originate from the supplied asset audit.

## Outputs

- external_predictions.csv: raw and calibrated class probabilities for evaluable patient-endpoint pairs, for each model.
- external_metrics_95CI.csv: AUROC, average precision, Brier score, log loss and 10-bin ECE with 1,000 patient-bootstrap percentile intervals. Draws missing an endpoint class are retried. Intervals are conditional on the frozen models and observed eligible cohort; they do not include model-training uncertainty.
- endpoint_eligibility.csv and external_label_alignment.csv: endpoint-specific missingness and label source records. Missing/equivocal/unmapped labels are excluded, never recoded negative.
- reliability_bins.csv: positive-class reliability for binary endpoints and top-label reliability for PAM50. PAM50 Brier is the unnormalised sum across four classes.
- local_tcga_replay.json: all 20 model replay checks against stored TCGA predictions.
- tissue_review_status.csv and tissue_documentation.json: technical QC is complete; histological tumor confirmation remains pending. All 391 coordinate hashes match the source QC report.
- summary.json: input hashes and interpretation limits.

## Interpretation limits

These are computational results for the existing technically filtered tissue cohort. Do not describe it as a pathologist-confirmed tumor-region cohort. Tumor assessment requires qualified human histology review of the existing contact sheets or WSIs, blinded to predictions and outcomes, with documented slide decisions. If eligibility changes, rerun the same frozen evaluation on the documented cohort and retain this version.

Image-to-assayed-block correspondence and encoder pretraining independence have not been established. PIK3CA negatives mean no qualifying reported mutation under the frozen source definition, not proof of absence under every mutation assay. The binary clinical receptor endpoints and PAM50 HER2-enriched subtype are distinct.

The report presents every pre-existing site model without selecting the best external score. Any future ensemble must be explicitly specified and evaluated as a new analysis. Calibration diagnostics do not authorize external recalibration on this same evaluation cohort.
