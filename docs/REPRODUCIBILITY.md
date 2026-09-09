# Reproducibility guide

## Scope

The public package supports two levels of reproduction:

1. **Analysis reproduction:** rerun Stages 10B–12 from the deposited patient-level feature matrix and frozen partitions.
2. **End-to-end reconstruction:** obtain TCGA-BRCA WSIs and authorized UNI resources, then rerun Stages 9A–12 with the required upstream cohort and endpoint manifests.

## Frozen cohort

- 345,937 patches
- 347 diagnostic WSIs
- 345 patients
- 1,024-dimensional UNI representations
- ER: 337 evaluable patients
- PR: 337 evaluable patients
- HER2: 314 evaluable patients
- PIK3CA: 345 evaluable patients
- PAM50 four-class: 335 evaluable patients

The analysis seed is `20260905`. Patch selection uses seed `2026`.

## Analysis-only reproduction

From the repository root, install `requirements-analysis.txt`. The deposited `stage10a_cohort_split_freeze_v1` directory contains the frozen patient feature matrix, index, endpoint master, and partition manifests. Run the PowerShell launchers for Stages 10B, 10C, and 11, followed by Stage 12 synthesis.

Each analysis stage validates upstream hashes and structural gates before producing outputs. Existing verified checkpoints may be skipped by the resumable routines.

## WSI reconstruction

Raw WSIs are not redistributed. Download the required diagnostic slides from the NCI GDC under `TCGA-BRCA`, preserving the site directory mapping expected by the Stage 9 scripts:

```text
<WSI_ROOT>/
├── A2_remaining/
├── C8/
├── E2/
└── BH/
```

Set environment variables before running Stage 9:

```powershell
$env:TCGA_BRCA_PROJECT_ROOT = (Get-Location).Path
$env:TCGA_BRCA_WSI_ROOT = "D:\path\to\targeted_sites_v1"
$env:PYTHON_EXE = "python"
```

Stage 9A applies the frozen 0.50-micrometre-per-pixel, 256-pixel patch policy, a maximum of 1,000 patches per slide, the locked colour threshold `bg_to_red_ratio < 1.10`, and structural threshold `mean_local_std >= 0.50`, without backfilling rejected patches. Stage 9B requires a CUDA-capable PyTorch environment and authorized access to `hf-hub:MahmoodLab/uni`.

Stage 10A additionally requires the frozen Batch01 and multisite patch-level embedding artifacts plus the full-cohort endpoint master. These large upstream files are not duplicated in this software release; their hashes are recorded in the Stage 10A metadata.

## Leakage controls

- All slides and patches from a patient remain in the same partition.
- Standardization, normalization, PCA, and hyperparameter selection are fitted only on the relevant training partition.
- Binary Platt scaling and multiclass temperature scaling use development patients only.
- Outer-test and site-held-out test patients remain untouched until final prediction.
- Bootstrap resampling is performed at the patient level.

## Public-release path sanitization

The production scripts originally ran in a fixed Windows workspace. For public release, only path resolution and launcher portability were changed: project and WSI roots are resolved relative to the repository or environment variables. Statistical procedures, thresholds, seeds, cohort counts, and model specifications are unchanged.
