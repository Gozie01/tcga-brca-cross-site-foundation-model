# Stage 9B Multisite UNI Pipeline v1

This package embeds the locked Stage 9A coordinates for A2, C8, E2, and BH
using the already validated Batch01 UNI configuration.

## Locked configuration

- Model: `hf-hub:MahmoodLab/uni`
- Output dimension: 1024
- CUDA inference: FP16 autocast
- Stored embeddings: float32
- Batch size: 32
- UNI input/preprocessing: resolved by `timm`, matching the validated Batch01 run
- Patch source: level-0 SVS coordinates at 0.50 MPP; intermediate RGB size 256 px
- No permanent patch-image corpus

## Expected cohort

| Site | WSIs/patients | Retained patches |
|---|---:|---:|
| A2 | 11 | 10,987 |
| C8 | 46 | 45,999 |
| E2 | 71 | 70,998 |
| BH | 126 | 125,180 |
| Total | 254 | 253,164 |

## Installation and run

Extract the ZIP directly into:

`<PROJECT_ROOT>`

Then run in PowerShell:

```powershell
Set-Location "<PROJECT_ROOT>"
Set-ExecutionPolicy -Scope Process Bypass
& ".\Stage9B_Multisite_UNI_Pipeline_v1\run_stage9b_multisite_uni.ps1"
```

The launcher processes A2, C8, E2, and BH and then performs the final audit.
Each slide receives an independently hashed checkpoint. If the run is stopped,
run the same command again; verified slide checkpoints are skipped.

## Outputs

Per-site files are written below:

`uni_multisite_embeddings_v1\<SITE>`

The final locked multisite artifact is written below:

`uni_multisite_embeddings_v1\multisite_frozen_v1`

Successful completion must print:

`Status : LOCKED_COMPLETE`

and all audit gates must show `PASS`.

Do not delete the slide checkpoints until the final multisite audit/freeze passes
and the frozen artifacts have been backed up.

## Recovery and troubleshooting

- Interruption or restart: rerun the same launcher.
- A corrupt/incomplete checkpoint is automatically recomputed.
- CUDA unavailable: the pipeline intentionally aborts rather than silently use CPU.
- CUDA out-of-memory: edit only the launcher's `--batch-size 32` value to 16.
- Missing WSI: confirm `A2_remaining`, `C8`, `E2`, and `BH` exist under
  `<WSI_ROOT>`.
- Do not edit Stage 9A frozen manifests or their metadata.
