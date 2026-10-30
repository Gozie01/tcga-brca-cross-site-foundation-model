Regenerates the first 32 indexed original TCGA patches from Batch01, A2, C8, E2 and BH.
Compares against verified stored UNI embeddings at batch sizes 32 and 8.
This is a deterministic diagnostic subset, not a representative quality or performance evaluation.
Checks original matrix/index hashes, cached weights and preprocessing configuration.
No downloads, installations, model training, CPTAC predictions or full feature extraction.
The report describes numerical agreement without imposing an unvalidated acceptance threshold.
Even exact agreement on 160 patches cannot prove historical weight-file identity or pretraining independence.

Extract under Desktop\BioMedical and run Start-UNICompatibility.ps1.
Original TCGA slides must be under D:\TCGA_BRCA_WSI; override -WSIRoot if relocated.
Original embeddings remain in TCGA_BRCA_Mutation_Audit.
Hashing the original matrices may take several minutes; GPU inference is small.
Upload the printed report ZIP even if an error is recorded.

Bundled tcga_classifier_replay.json records successful independent replay of all 20 original
site-holdout classifiers and calibrators from the uploaded assets. It contains no CPTAC outcomes.
