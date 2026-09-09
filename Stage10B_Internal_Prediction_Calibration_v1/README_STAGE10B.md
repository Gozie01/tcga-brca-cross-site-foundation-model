# Stage 10B Internal Prediction and Calibration v1

This pipeline evaluates the frozen Stage 10A patient-level UNI features on all
five endpoints using leakage-safe nested cross-validation.

## Endpoints

- ER binary: 337 patients
- PR binary: 337 patients
- HER2 binary: 314 patients
- PIK3CA binary: 345 patients
- PAM50 primary four-class: 335 patients

## Protocol

- Frozen endpoint-specific outer five-fold patient assignments from Stage 10A.
- Standardization fitted only on each outer-training set.
- L2 logistic-regression linear probe with class balancing.
- `C` selected inside the outer-training set by stratified inner CV.
- Binary selection metric: average precision.
- PAM50 selection metric: negative log loss.
- Sigmoid calibration learned by stratified CV entirely inside the outer-training
  set; the outer test fold is untouched until final inference.
- Raw and calibrated predictions are both preserved.
- Metrics include AUROC, AUPRC, balanced accuracy, macro-F1, MCC, Brier score,
  log loss, and 10-bin ECE; binary endpoints also include sensitivity and
  specificity.
- Patient-level 95% bootstrap confidence intervals use 1,000 replicates.

Stage 10B is internal evaluation only. The frozen A2, C8, E2, and BH TCGA
site holdouts are reserved for Stage 10C.

## Run

Extract the package directly into:

`<PROJECT_ROOT>`

Then run:

```powershell
Set-ExecutionPolicy -Scope Process Bypass
Set-Location "<PROJECT_ROOT>"
& ".\Stage10B_Internal_Prediction_Calibration_v1\run_stage10b_internal_prediction_calibration.ps1"
```

Each endpoint/fold is checkpointed with its raw model, calibrated model,
predictions, metrics, hashes, and task fingerprint. Rerunning safely skips
verified checkpoints.

Outputs are written to `stage10b_internal_prediction_calibration_v1`.
Successful completion must report `LOCKED_COMPLETE` with every gate `PASS`.
