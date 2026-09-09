# Stage 9A — Multi-site patch-coordinate extraction

This package extends the frozen Batch01 method to the 254 acquired TCGA-BRCA
slides from A2, C8, E2, and BH.

## Method preserved

- 0.50 MPP physical normalization
- 256 × 256 output patches
- non-overlapping coarse grid
- low-resolution tissue fraction ≥ 0.50
- low-resolution green fraction ≤ 0.05
- high-resolution tissue fraction ≥ 0.35
- white/background fraction ≤ 0.60
- green/cyan artifact fraction ≤ 0.01
- focus recorded but not used for rejection
- maximum 1,000 initially accepted coordinates per slide
- deterministic per-slide seed derived from global seed 2026 and filename
- color rejection when `bg_to_red_ratio ≥ 1.10`
- structural rejection when `mean_local_std < 0.50`
- final-gate rejects are not backfilled
- coordinate-only outputs; patches remain on-the-fly

The exact color metric is `((mean_G + mean_B)/2) / max(mean_R, 1e-6)`.
The exact structural metric is the mean of the 9×9 local grayscale standard
deviation computed with OpenCV box filters.

## Run order

1. Copy the complete `Stage9A_Multisite_Patch_Pipeline` folder into:
   `<PROJECT_ROOT>`
2. Open PowerShell in that audit folder.
3. Run the A2 pilot:

   ```powershell
   Set-ExecutionPolicy -Scope Process Bypass
   .\Stage9A_Multisite_Patch_Pipeline\run_stage9a_A2_pilot.ps1
   ```

4. Do not launch the remaining sites until the A2 audit ends with `Status: LOCKED`
   and every displayed gate is `PASS`.
5. After A2 approval, run:

   ```powershell
   .\Stage9A_Multisite_Patch_Pipeline\run_stage9a_remaining_sites.ps1
   ```

## Resumption

Each slide is written as an atomic CSV checkpoint plus JSON metadata and hash.
Repeating a command verifies and skips completed slides. A failed slide is
processed again. Never edit checkpoint files manually.

## Dependencies

The same Python environment used for Batch01 should be used. Required imports:
`numpy`, `Pillow`, `openslide-python`, `opencv-python`, and `pandas`.

## Output

Outputs are written below:
`multisite_patch_coordinates_v1\<SITE>`.

For each site the audit creates a combined pre-final-gate manifest, retained
manifest, rejected manifest, WSI summary, SHA-256 hashes, and locked metadata.
