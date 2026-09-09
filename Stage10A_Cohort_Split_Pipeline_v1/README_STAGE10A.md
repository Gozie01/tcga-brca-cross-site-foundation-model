# Stage 10A Cohort and Split Freeze v1

This pipeline consolidates the locked Batch01 and multisite UNI artifacts,
reproduces the validated patch-mean → slide-mean → patient-mean aggregation,
attaches the frozen endpoint labels, and freezes leakage-safe modeling splits.

## Audited expectations

- 345,937 patch embeddings
- 347 WSIs
- 345 unique patients
- Zero Batch01–multisite patient or slide overlap
- Patient feature matrix: `(345, 1024)`, float32
- ER: 337 evaluable
- PR: 337 evaluable
- HER2: 314 evaluable
- PIK3CA: 345 evaluable
- PAM50 primary four-class: 335 evaluable

## Frozen split families

1. Endpoint-specific, patient-stratified internal five-fold assignments.
2. Endpoint-specific TCGA site-held-out experiments for A2, C8, E2, and BH.
3. For each site-held-out experiment, the held-out site is test-only. The remaining
   sites are divided into training and validation sets deterministically within
   site/class strata.
4. Every slide and patch belonging to a patient remains in the same partition.

## Run

Extract this ZIP directly into:

`<PROJECT_ROOT>`

Then execute:

```powershell
Set-Location "<PROJECT_ROOT>"
Set-ExecutionPolicy -Scope Process Bypass
& ".\Stage10A_Cohort_Split_Pipeline_v1\run_stage10a_cohort_split_freeze.ps1"
```

The script reads the existing embedding matrices in memory-mapped mode. It does
not alter or duplicate the patch-level matrices. The resulting patient feature
matrix is small (345 × 1024).

Successful completion must report `Status: LOCKED_COMPLETE` and every gate as
`PASS`. Outputs are written to `stage10a_cohort_split_freeze_v1`.

Do not begin classifier training until this freeze passes.
