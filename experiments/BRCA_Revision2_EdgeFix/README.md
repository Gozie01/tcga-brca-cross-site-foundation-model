# BRCA revision 2: edge-patch repair and guarded resume

The full extraction log completed 178/347 slides, then stopped while processing
slide 179, TCGA-BH-A0W5, because the new Hibou reader rejected a frozen patch
that extended beyond a level-zero slide boundary. The completed slide means
are preserved in the existing output folder.

## Cause and correction

The original Stage9B UNI `read_patch` directly requested the full frozen region
through OpenSlide, converted its RGBA result to RGB, then resized to 256 pixels
with Lanczos. Its reader did not reject a request crossing a right or bottom
edge. The first Hibou extractor introduced an incompatible upper-bound guard.

The repaired extractor restores the original reader behavior. Coordinates,
requested region sizes, patch counts, model checkpoint, processor, float32
inference, batch size 32, and equal-slide patient pooling remain fixed. It does
not clip, shift, exclude, replace, or resample coordinates. OpenSlide returns
the requested region; its handling of pixels outside the image is retained
through the same RGB conversion as the original UNI script.

The original source excerpt and its full-file hash are in UNI_EDGE_REFERENCE.txt.
OpenSlide Python API: https://openslide.org/api/python/

## Preserve completed work with explicit lineage

The launcher installs three side-by-side repair files. It leaves the original
extract_hibou.py untouched and first checks its known SHA256 and repair hashes.
Twelve regression tests run locally before any checkpoint migration.

The repair uses a new extraction protocol version and code hashes. Migration
only accepts the exact known old protocol, including model identity and all
other settings. Before changing extraction metadata it checks every inherited
checkpoint's vector hash, dtype, dimensions, finite/nonzero values, patient,
slide identity, patch count, coordinate hash, full WSI hash and in-bounds status.
Only entirely in-bounds inherited patches are eligible, because their pixel
path is identical before and after the fix.

After all checks pass, the old protocol and metadata are backed up under
boundary_v1_backup. The new protocol fingerprint is recorded in the checkpoint
metadata. The `.npy` vector files are retained byte-for-byte without GPU
re-extraction. boundary_checkpoint_migration.json records this lineage. An
interrupted metadata transition can be retried with the same launcher.

Each normal resume also checks the slide and vector hashes. As a result, source
images are reread during verification; initial progress may take time. The old
178 vectors are not recomputed. The next new GPU extraction is slide 179.

## Run

Save BRCA_Revision2_EdgeFix.zip in C:\Users\Gozie\Desktop\BioMedical, then:

```powershell
Set-Location "C:\Users\Gozie\Desktop\BioMedical"
Expand-Archive -LiteralPath ".\BRCA_Revision2_EdgeFix.zip" -DestinationPath "." -Force -ErrorAction Stop
powershell -NoProfile -ExecutionPolicy Bypass `
    -File ".\BRCA_Revision2_EdgeFix\Resume-BRCARevision2.ps1" `
    -ProjectRoot "C:\Users\Gozie\Desktop\BioMedical"
```

Keep the entire brca_revision2_hibou_complete folder. If interrupted, run the
same launcher again. Use the new launcher for subsequent resumes.

## New audit outputs

- slide_edge_audit.csv records slide dimensions, patch counts, edge-crossing
  counts, fully-outside counts and maximum right/bottom overflow.
- edge_patch_coordinates.csv identifies every frozen patch crossing an edge.
- feature_freeze.json includes hashes of both audit tables and migration lineage.

Final gates still require 347 slides, 345,937 patches and a finite 345 x 768
patient feature matrix. The complete output is zipped automatically. Return
brca_revision2_hibou_complete.zip for verification and paired encoder evaluation.

The repair has synthetic edge-read and checkpoint-migration regression coverage.
The real edge-slide read and actual 178-checkpoint preflight must run locally;
this package does not claim that the remaining extraction has already completed.
