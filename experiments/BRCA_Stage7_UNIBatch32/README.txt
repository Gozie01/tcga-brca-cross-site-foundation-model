UNI-only extraction at original TCGA batch size 32.
Extract under Desktop\BioMedical and run Start-UNIBatch32.ps1.
Exactly the same 391 CPTAC slides and 368,985 frozen retained patches.
Pinned UNI weights, transforms, float16 autocast, and equal-slide patient pooling preserved.
Hibou is not loaded or rerun. No downloads, training, outcome selection or molecular predictions.

Outputs go to this NEW folder's features_uni32_v1. Previous Stage4 features_v1 remains untouched.
Allow 3 GiB free on this drive. New UNI patch arrays total approximately 1.51 GB.
Each slide is checkpointed and verified on resume. Rerun the same command if interrupted.
Input slide hashes must match the completed batch8 run's source registry.

The original 128-patch batch8 probe is repeated at batch8 to verify the cached adapter first.
FULL COHORT EXTRACTION then uses batch32, including each slide's final shorter batch.
The bundled Stage6 report documents exact batch32 agreement on 160 original TCGA patches;
this is numerical compatibility evidence, not proof of historical weight-file identity.

The compact report includes complete slide and patient means, checkpoint metadata, and a
patient-feature comparison against completed batch8 outputs. It excludes the large patch arrays.
Differences in features are reported without making an untested clinical-effect claim.
Upload the printed brca_uni32_extraction_report_*.zip; retain the features_uni32_v1 folder.
