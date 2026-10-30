Stage5 collects existing artifacts needed to choose and validate the external prediction protocol.
It does not load pickle/joblib models, train, download, extract patches, or generate predictions.

Extract under Desktop\BioMedical; run Start-AssetAudit.ps1. Upload the printed report ZIP.
The scan is limited to BioMedical subfolders starting TCGA_BRCA, BRCA_Revision, or
BRCA_ExternalValidation_Stage1 (including nested project folders).
Includes relevant scripts, classifier files, metadata, splits, small patient features, label tables
and prediction records. It excludes whole slides, patch embedding arrays, environments and cache
directories. Files over 20 MiB and files beyond the 200 MiB collection limit are listed as skipped.
No model is chosen automatically; collected checkpoints may represent different folds/protocols.

Next review will establish original feature aggregation, encoder provenance, classifier and
calibrator availability, class ordering, training split identities, and CPTAC label linkage.
Missing evidence remains explicit; collection success is not a scientific validation result.
