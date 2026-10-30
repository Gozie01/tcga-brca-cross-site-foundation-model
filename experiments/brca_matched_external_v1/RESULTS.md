# Matched external UNI versus Hibou-B comparison

All 20 Hibou-B classifiers were trained using the original TCGA endpoint/site train partitions, with the original training-only C grid, preprocessing, tuning criterion, and seed. Calibration used only the corresponding TCGA validation partition. All models were locked before external evaluation. UNI uses its pre-existing frozen checkpoints and verified predictions. No external tuning, recalibration, model selection, or ensemble was performed.

| Endpoint | External n | UNI calibrated AUROC range | Hibou-B calibrated AUROC range |
|---|---:|---:|---:|
| ER_binary | 118 | 0.856–0.908 | 0.738–0.823 |
| HER2_binary | 95 | 0.436–0.614 | 0.421–0.530 |
| PAM50_primary_4class | 114 | 0.815–0.854 | 0.694–0.724 |
| PIK3CA_binary | 120 | 0.565–0.733 | 0.492–0.552 |
| PR_binary | 112 | 0.690–0.824 | 0.684–0.759 |

Ranges describe the four source models, not confidence intervals. Per-model intervals and paired differences are provided separately. UNI has higher calibrated AUROC point estimates in all 20 comparisons. Paired Hibou-minus-UNI AUROC intervals exclude zero in the UNI-favouring direction for PAM50 in 4/4 source models, ER in 3/4, PIK3CA in 2/4, and PR in 1/4; none for HER2. These are exploratory, unadjusted intervals across overlapping models and endpoints, not independent confirmatory replications. Both encoders perform weakly for HER2. Probability calibration is not uniformly portable.

## Files

- matched_external_metrics.csv: per-encoder, per-model raw/calibrated metrics and paired-bootstrap-derived marginal intervals.
- paired_encoder_contrasts.csv: Hibou minus UNI paired changes with 95% percentile intervals, 1,000 valid patient resamples. Higher AUROC/AP is better; lower Brier/log loss/ECE is better.
- paired_calibration_effects.csv: calibrated minus raw changes within each encoder/model, using identical patient resamples. No external recalibration is applied.
- calibration_diagnostics.csv: descriptive calibration-in-the-large, joint recalibration intercept and slope. These regression fits diagnose predictions; they are not applied to reported probabilities. PAM50 uses one-vs-rest diagnostics.
- hibou_external_predictions.csv and hibou_tcga_site_predictions.csv: auditable probabilities.
- tcga_heldout_context_metrics.csv: source-domain held-out point estimates for these exact models; not a rerun of all internal five-fold comparisons.
- tcga_fit_audit.csv, models_locked.json, protocol.json: model choices, sizes, source hashes and fixed protocol.
- models and inputs: frozen models and the data needed to reproduce this analysis.

## Outstanding evidence

Both encoders use the same 391 external slides, patient ordering, source-file hashes, patch counts, and coordinate hashes; TCGA Hibou feature hashes and patient indexing pass the supplied freeze. The original complete TCGA archive was subsequently recovered. Its extraction protocol SHA matches the feature freeze, and its recorded Hibou revision and all five model/configuration/code hashes match the verified external pilot records. Official processor settings and float32 precision also match the external full-extraction summary, which reports zero pilot-probe difference. Recorded checkpoint provenance is therefore resolved. A separate TCGA Hibou numerical replay across batch sizes has not been performed; TCGA used batch 32 and the external pilot used batch 8. See hibou_provenance_resolution.json and the preserved source protocols. This addendum changes documentation, not the prediction results.

Human tumor-content review remains pending. Matching to the molecularly assayed block and encoder pretraining independence are not proven. These computational results therefore remain conditional on those limitations. Missing/equivocal/unmapped endpoint labels are excluded under the frozen reference definitions. No clinical utility or assay-replacement claim is supported.

Binary positive-class ECE is provided alongside the original manuscript top-label ECE (ece_10bin); they are distinct definitions. Multiclass Brier is the unnormalised four-class sum. Patient bootstrap intervals condition on fixed trained models, preserve patient pairing, and do not include training uncertainty.

## Reproduce

Install requirements.txt in a separate environment. Run python run_matched.py, then python summarize.py. No GPU, slide downloads, or feature extraction are required. A clean run fits Hibou models on TCGA only; an existing package reuses the bundled models. Existing UNI results are supplied and were verified against stored TCGA predictions in the earlier external evaluation.
