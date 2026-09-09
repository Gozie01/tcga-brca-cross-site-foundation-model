# TCGA-BRCA Cross-Site Foundation-Model Evaluation

Reproducible evaluation of UNI pathology foundation-model representations across TCGA-BRCA collection sites. The pipeline performs WSI QC, patch selection, UNI extraction, patch-to-patient pooling, frozen patient-level splits, linear-probe evaluation, calibration, institution-held-out evaluation, and statistical ablations.

This repository contains no WSI files, GDC-controlled data, credentials, UNI weights, patch images, or incomplete checkpoints. A2, C8, E2, and BH are institution-held-out sites within TCGA-BRCA, not independent clinical cohorts.

## Install

```bash
conda env create -f environment.yml
conda activate tcga-brca-cross-site-foundation-model
```

Obtain TCGA-BRCA WSIs through the GDC using the supplied manifests and obtain UNI weights from the authorized provider. Set `TCGA_BRCA_WSI_ROOT` to the reconstructed WSI directory and `TCGA_BRCA_PIPELINE_ROOT` to a writable pipeline-output directory.

## Reproduce

Run the Stage 9A, 9B, 10A, 10B, 10C, and 11 wrappers in `scripts/` in order. Each accepts `--root` where applicable. Generate manuscript tables and figures with `python figures/generate_manuscript_figures_and_tables.py` after staging the Stage 10-11 outputs under `stage12_work/`.

The frozen results used for release verification are in `derived_data/`. Run `python tests/verify_frozen_release.py` to validate all included CSV and metadata checksums. See `docs/REPRODUCIBILITY.md` and `docs/DATA_DICTIONARY.md`.

Code is MIT-licensed. Original documentation and derived tables are CC BY 4.0; TCGA data remain subject to GDC/TCGA terms, and UNI weights must be obtained under their original license.
