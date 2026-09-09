# Stage 11 Statistical Analysis and Ablations v1

This package performs the frozen statistical-analysis and representation-
ablation phase after successful completion of Stages 10A, 10B, and 10C.

It produces:

- paired patient-level bootstrap comparisons of raw versus calibrated internal
  predictions;
- paired internal-OOF versus site-held-out domain-shift comparisons;
- cross-site heterogeneity and worst-site summaries;
- leakage-safe representation ablations using L2-normalized UNI features and
  training-fold-fitted PCA at 32, 64, and 128 dimensions;
- paired bootstrap comparisons of every ablation against the frozen full
  1024-dimensional Stage 10B reference.

There are 100 ablation experiments: five endpoints, five outer folds, and four
ablation variants. Each experiment is checkpointed and safely resumable.

Run from PowerShell:

```powershell
Set-Location "<PROJECT_ROOT>"
Set-ExecutionPolicy -Scope Process Bypass
& ".\Stage11_Statistical_Analysis_Ablations_v1\run_stage11_statistics_ablations.ps1"
```

The final status must be `LOCKED_COMPLETE` with all audit gates marked `PASS`.
