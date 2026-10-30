# Full-cohort CPTAC technical QC

Extract under C:\Users\Gozie\Desktop\BioMedical. Run Start-BRCAExternalFullQC.ps1 with PowerShell. Default WSI root is D:\CPTAC_BRCA_WSI and the existing .venv-brca-revision environment is used.

Processes all 391 slides from 120 patients, using the exact supplied Stage9A policy. Patient and site identifiers come from the acquisition manifest. No endpoint labels are opened. No encoder features or predictions are computed.

The target MPP, patch cap, tissue/artifact/colour/structural thresholds, seed, Lanczos resizing, and edge handling are unchanged. Blurred patches can pass: focus is audit-only, as in TCGA. Preview review is technical only and does not certify tumor content.

Per-slide coordinate checkpoints have checksums and can be reused on rerun. A run lock binds the runner, policy, frozen manifest and WSI root. Source sizes and MPP are checked; this is not a fresh full-file checksum of each WSI. Do not launch concurrent copies. Four-slide pilot outputs are preserved; this separate cohort run recomputes those slides once. Rerun the same command after interruption. Outputs are in BRCA_ExternalValidation_Stage2_FullQC/full_qc_v1.

At completion, upload the brca_stage2_full_qc_*.zip printed in the console. This includes 16 evenly spaced retained-patch previews per slide (when available), all retained coordinates and slide-level counters. Zero-retained slides remain in the summary and are explicitly flagged. Completion is technical QC only: visual/pathology review, linkage limitations and pretraining overlap remain unresolved. Do not tune QC thresholds to CPTAC labels or scores.

Existing TCGA data, Stage1 acquisition and the pilot are not modified. No model download or package installation occurs.
