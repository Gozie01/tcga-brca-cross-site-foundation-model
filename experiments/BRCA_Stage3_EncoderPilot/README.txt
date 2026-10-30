Offline diagnostic extraction of up to 128 patches (32 per slide) using UNI and Hibou-B.
Four slides: first alphabetical and lowest patch count at each of the two resolutions.
Uses frozen retained coordinates without QC changes. No labels, predictions or training.
Pinned cached model hashes come from the uploaded preflight report.
UNI follows original timm preprocessing with CUDA fp16 autocast.
Hibou follows original AutoImageProcessor and pooler_output with float32 inference.
Models run sequentially with batch size 8. No installation or model downloads.
Output features are diagnostic subsets, NOT complete slide or patient features.
Recorded extraction speed excludes slide reads/model loading and is not end-to-end throughput.
Matching the original TCGA weight identity and preprocessing/library provenance remains pending.
Extract under Desktop\BioMedical and run Start-EncoderPilot.ps1; upload printed report ZIP even on failure.
