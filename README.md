# Cross-site generalization of pathology foundation-model representations in TCGA-BRCA

This repository contains the portable analysis code and frozen derived outputs supporting the manuscript **“Cross-Site Generalization and Calibration of Pathology Foundation-Model Representations for Breast Cancer Molecular Phenotype Prediction.”**

The study uses 345,937 quality-controlled patches from 347 diagnostic whole-slide images representing 345 patients in TCGA-BRCA. Frozen 1,024-dimensional UNI patch representations are pooled hierarchically (patch to slide to patient) and evaluated for ER, PR, HER2, PIK3CA, and four-class PAM50 prediction. Evaluation comprises patient-separated internal five-fold testing, validation-only probability calibration, four TCGA tissue-source-site holdouts (A2, C8, E2, and BH), patient-level bootstrap inference, and representation ablations.

## Interpretation boundary

The site-held-out experiments evaluate institution-associated domain shift **within TCGA-BRCA**. They are not independent external-cohort validation. The results demonstrate morphological associations and research-level patient prediction; they do not support replacement of clinical molecular assays.

## Repository contents

- `Stage9A_Multisite_Patch_Pipeline/`: locked patch-coordinate extraction and QC audit.
- `Stage9B_Multisite_UNI_Pipeline_v1/`: frozen UNI representation extraction and audit.
- `Stage10A_Cohort_Split_Pipeline_v1/`: hierarchical pooling and split construction.
- `Stage10B_Internal_Prediction_Calibration_v1/`: internal five-fold prediction and calibration.
- `Stage10C_Site_Holdout_Domain_Shift_v1/`: TCGA tissue-source-site holdout evaluation.
- `Stage11_Statistical_Analysis_Ablations_v1/`: paired inference and representation ablations.
- `Stage12_Results_Synthesis_v1/`: manuscript-oriented tables and figures.
- Lower-case `stage10*` and `stage11*` directories: frozen derived outputs used for reproducibility.
- `manuscript_assets/`: source script and vector figures used in the manuscript.
- `docs/`: reproducibility, data dictionary, provenance, and release instructions.

## Quick reproduction of reported analyses

Python 3.11 is recommended. Create an environment and install the analysis dependencies:

```bash
python -m venv .venv
# Windows: .venv\Scripts\Activate.ps1
# Linux/macOS: source .venv/bin/activate
python -m pip install --upgrade pip
python -m pip install -r requirements-analysis.txt
```

The repository includes the frozen patient-level feature matrix and partition manifests needed to rerun Stages 10B, 10C, 11, and 12 without redistributing TCGA whole-slide images:

```powershell
.\Stage10B_Internal_Prediction_Calibration_v1\run_stage10b_internal_prediction_calibration.ps1
.\Stage10C_Site_Holdout_Domain_Shift_v1\run_stage10c_site_holdout.ps1
.\Stage11_Statistical_Analysis_Ablations_v1\run_stage11_statistics_ablations.ps1
python .\Stage12_Results_Synthesis_v1\build_stage12_results_synthesis.py --root .
```

Stages 9A, 9B, and 10A require source WSIs and upstream frozen artifacts. See [`docs/REPRODUCIBILITY.md`](docs/REPRODUCIBILITY.md).

## Data sources

Raw whole-slide images and source clinical/molecular data are not redistributed. They must be obtained from the NCI Genomic Data Commons under the TCGA Breast Invasive Carcinoma (`TCGA-BRCA`) project and, where indicated in the manuscript, from the TCGA PanCancer Atlas resources. Users are responsible for complying with the applicable source-data terms.

## Reproducibility safeguards

- Patient-level separation keeps all slides and patches from a patient in one partition.
- Preprocessing, PCA, model selection, and calibration are fitted without access to held-out test patients.
- Frozen metadata record expected sample counts, hashes, software versions, and audit gates.
- Random seeds and prespecified hyperparameter grids are embedded in the stage scripts.
- `checksums_sha256.txt` records the public release files.

## Licences and citation

Code is released under the MIT License. Original documentation and derived data in this repository are provided under CC BY 4.0; source TCGA data and UNI model resources remain governed by their original terms. See `DATA_LICENSE.md` and `CITATION.cff`.

## Contact

Correspondence: Dong-Seong Kim, Department of IT Convergence Engineering, Kumoh National Institute of Technology, Gumi, Republic of Korea; `dskim@kumoh.ac.kr`.
