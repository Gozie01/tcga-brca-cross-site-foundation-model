Stage 4: offline full-cohort UNI and Hibou-B patch feature extraction.

Extract this folder under C:\Users\Gozie\Desktop\BioMedical and run Start-FullExtraction.ps1.
Uses the existing .venv-brca-revision environment and original Stage2 QC report.
No installations, downloads, training, molecular labels or predictions.

Scope: exactly 391 slides, 120 patients and 368,985 patches per encoder.
Uses the pilot's pinned weights, batch size 8, preprocessing and precision.
Both encoders repeat all 128 reference pilot patches before extraction; exact equality is required.
Preserves the five boundary-crossing patches and all frozen QC decisions.
Each input slide is hashed before processing for each encoder. This adds disk-reading time.
The run stops on failures; no failed slides are silently excluded.

Restart: run the same command again. Completed slide arrays are validated by input identity,
coordinate checksum, array checksum, dimensions and finite/nonzero values before reuse.
Only the interrupted slide needs extraction again. A changed environment or damaged committed
checkpoint stops the run for review rather than silently mixing results.
Concurrent runs in this folder are refused by an operating-system lock.

Allow at least 5 GiB free on the drive holding this folder. Patch arrays alone total approximately
2.64 GB for both encoders. Retain features_v1; the compact upload ZIP excludes these large arrays.
Keep the computer awake and D:\CPTAC_BRCA_WSI accessible. Extraction runs UNI then Hibou.
Upload the printed brca_full_extraction_report_*.zip, including after a reported failure.

The report includes 391 slide means and 120 equal-slide patient means per encoder.
Patient pooling is a provisional convenience output. Compatibility with the original UNI predictor
must be established before using it for predictions. All patch features remain available for exact
reconstruction of the required original aggregation.
Original TCGA weight identity, human tumor review, specimen-alias linkage and encoder pretraining
overlap remain unresolved. Successful extraction does not resolve these scientific questions.

Local validation covers checkpoint corruption, changed inputs, partial writes, invalid arrays,
patient pooling, the real QC coordinate table and Python syntax. Actual CUDA inference and full
slide reading must run on the user's computer.
