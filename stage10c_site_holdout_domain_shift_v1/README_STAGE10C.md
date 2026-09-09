# Stage 10C Site-Held-Out Domain Shift v1

Runs 20 within-TCGA site-held-out experiments: five endpoints × A2, C8, E2, and BH. The held-out
institution is test-only. Model selection uses training patients only;
calibration uses the frozen validation patients only. Binary endpoints use
Platt sigmoid calibration and PAM50 uses one-parameter temperature scaling.

Extract into the project root, then run:

```powershell
Set-Location "<PROJECT_ROOT>"
Set-ExecutionPolicy -Scope Process Bypass
& ".\Stage10C_Site_Holdout_Domain_Shift_v1\run_stage10c_site_holdout.ps1"
```

The run is checkpointed per endpoint/site and reports raw/calibrated metrics,
95% bootstrap confidence intervals, and AUROC/AUPRC gaps relative to Stage 10B
internal OOF performance.
