# October 2026 experiment update

This update adds archived study code, derived patient-level data, fitted downstream classifiers, predictions and results. It extends the original TCGA release with CPTAC breast-cancer evaluation and a matched UNI/Hibou-B comparison.

## Included analyses

| Path under `experiments/` | Purpose |
| --- | --- |
| `brca_external_frozen_v1/` | Frozen UNI external evaluation: 120 patients, 20 existing site-holdout models, no external fitting or calibration |
| `brca_matched_external_v1/` | Matched external UNI/Hibou-B comparison, TCGA-only fitting/calibration and paired bootstrap contrasts |
| `BRCA_Revision2_Pipeline/` | Revision analyses, frozen TCGA inputs, controls and encoder comparison code |
| `BRCA_ExternalValidation_Stage1/` | CPTAC acquisition, annotation linkage and endpoint freezing |
| `BRCA_ExternalValidation_Stage2_FullQC/` | Technical image QC code and summary |
| `BRCA_Stage3_*`, `BRCA_Stage4_*`, `BRCA_Stage6_*`, `BRCA_Stage7_*` | Extraction/preflight, numerical compatibility and batch-size evidence |
| `BRCA_AI_Tissue_Assessment/` | Exploratory tissue segmentation, overlap audit, technical-QC outputs and single-case sensitivity |
| `brca_stain_provenance*`, `trace_*`, `tcga_qc_reference*`, `source_policy/` | Provenance, patch trace and QC-policy evidence |

Other pilot and audit directories retain historical evidence; they are not additional independent validations. Historical paths, dates, model provenance, protocols and outputs have been preserved. Local absolute paths describe the original execution environment and require remapping for new extraction runs. For legacy programs accepting `--project-root`, the sibling stage folders now sit under `experiments/`; check each program's arguments before running it. This is an archival update, not a claim of a portable one-command WSI pipeline.

## Analysis environments and reruns

Use the per-analysis `requirements.txt` in the external and matched folders in separate Python 3.11 environments. Do not replace these with the original repository root requirements: recorded versions differ. The exploratory tissue run records Python 3.11.9, tiatoolbox 2.1.3 and torch 2.12.1+cu130 in `run_info.json`; a complete environment lock and its generating scripts were not supplied. Extraction also uses transformers, huggingface_hub and platform-compatible OpenSlide/PyTorch installations; the original root environment alone is insufficient.

The fitted UNI and Hibou classifiers are included in their respective `models/` directories, distinct from pretrained encoder weights, which are excluded. Their hashes were compared with `summary.json` and `tcga_fit_audit.csv`: all 40 matched. Joblib files are Python serialization; use only trusted release artifacts in a compatible environment.

To rerun an external evaluation, copy that analysis directory to a scratch location first, install its own requirements and run `python run_external.py` or `python run_matched.py` from the copy. These programs write outputs into their analysis directory; do not overwrite archived results. The matched program can load the supplied fitted classifiers. End-to-end extraction and numerical replay were not rerun while assembling this release.

## Interpretation and remaining gaps

- CPTAC is an external cohort; TCGA TSS holdouts remain within-TCGA domain-shift tests.
- All four overlapping site-holdout models are reported separately; they are not independent replications or an ensemble.
- Expert tumour-region confirmation remains pending in the archived protocols. Automated technical QC and segmentation are not pathologist review.
- The TCGA-5L-AAT0 sensitivity removes evaluation rows only. It does not retrain models; that patient can remain in other folds' training/calibration partitions.
- Generating scripts for the October 29 tissue pilot/inference, technical-QC segmentation and single-case sensitivity were not found in the supplied packages. Their outputs are released as archived evidence, not fully regenerable workflows.
- Full patch-coordinate tables and patch embeddings are not included. Exact end-to-end extraction replay still needs those tables or a separately identified archive. Raw WSIs, patch images, pretrained encoders, downloaded third-party repositories, PDFs and unrelated projects are excluded.
- Assayed-block correspondence and encoder pretraining independence are not established by these files.

## Verification

`UPDATE_FILES_SHA256.json` records every file delivered in this update except itself. Run `python scripts/verify_october_update.py` from the repository root. It checks delivered bytes, parses Python sources without executing them, and verifies the 40 classifier hashes against original audit records. These checks do not establish numerical replay or scientific validity.

The installer removes tracked PDFs from the current Git index while preserving local files. Earlier commits retain their historical PDFs; Git history is not rewritten. New study PNGs were not supplied in the reviewed packages and are not added by this update.
