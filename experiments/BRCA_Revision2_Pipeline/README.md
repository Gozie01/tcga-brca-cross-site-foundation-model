# BRCA revision 2: image audit, second encoder and external-cohort preparation

Start with the offline asset audit below. It needs only your existing Python,
checks the exact 347 required SVS filenames and current environment, and writes
a small return ZIP. It does not install packages, download weights, extract
features or change your original project files.

## First command

After downloading and extracting BRCA_Revision2_Pipeline.zip to Downloads:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
    -File "$env:USERPROFILE\Downloads\BRCA_Revision2_Pipeline\Start-BRCARevision2.ps1" `
    -Mode Audit `
    -ProjectRoot "$env:USERPROFILE\Desktop\BioMedical"
```

The project-root example comes from the earlier extraction metadata. If the
images were moved, supply their current parent folder or add another root:

```powershell
    -AdditionalWsiRoots "D:\your_actual_WSI_folder"
```

Return the generated `brca_revision2_asset_audit_YYYYMMDD_HHMMSS.zip`. A missing
or duplicate slide is listed individually. Overlapping search roots are
deduplicated. Whole-slide contents are not read or hashed during this audit.
READY_FOR_PILOT means image filenames are available; model access and working
GPU/OpenSlide dependencies still need a live pilot.

## What is already available

- All 345,937 original coordinates were recovered from the two hash-verified
  UNI embedding indices: 92,773 initial patches and 253,164 multisite patches.
- The package preserves 345 patients, 347 WSIs, original labels, patient order,
  endpoint-specific folds and the original site train/validation/test manifests.
- The selected comparator is Hibou-B. It uses a different pathology encoder
  while retaining the same physical tissue samples and hierarchical mean pooling.
- A matched evaluator is implemented for all five endpoints: 25 internal
  endpoint/fold fits and 20 site-holdout fits, with 1,000 paired patient
  bootstrap resamples for UNI-versus-Hibou contrasts.
- Completed revision-1 controls are used in manuscript text, tables and a
  replacement matched-calibration figure. No second-encoder values are invented.

## Current completion boundaries

See STATUS.json for the exact status. The first revision controls and recovered
coordinates are complete. Real Hibou extraction and the encoder comparison have
not run: the WSIs and model weights are not in this workspace. The public
CPTAC metadata API probe timed out, so independent-cohort patient mapping and
endpoint eligibility remain unverified. No external-validation result exists.

Ten tests passed: three offline inventory/metadata tests, five core pipeline
tests, and two manuscript-patch tests. The fitting adapters reproduced the
released UNI recipes for one internal fold and one site holdout for every
endpoint. Synthetic slide tests verify coordinate reading and equal-slide
patient pooling. These do not replace a live Hibou checkpoint/GPU pilot.
The PowerShell wrapper was reviewed but could not be executed here.

## Model access and environment after the audit

Hibou-B's publisher page currently requires a Hugging Face access agreement:
https://huggingface.co/histai/hibou-b

Use your own account to review the publisher conditions and obtain access.
The pipeline uses the standard Hugging Face cached credentials or a local
snapshot. No token belongs in command-line arguments, logs or return ZIPs.

Use an isolated Python 3.11+ environment. Install a compatible PyTorch and
torchvision GPU build using https://pytorch.org/get-started/locally/, then:

```powershell
python -m pip install -r .\requirements-core.txt
python -m pip install -r .\requirements-model.txt
```

The model adapter follows the official processor/AutoModel example and uses
pooler_output. Its live dependency compatibility is checked by the pilot;
the first inventory command does not require these installations. On Windows,
openslide-bin supplies the OpenSlide binaries when compatible with the Python
environment. If OpenSlide import fails, retain that error for the next step.

## Pilot, production extraction and comparison

After an audit reports READY_FOR_PILOT and dependencies/model access are ready:

```powershell
.\Start-BRCARevision2.ps1 -Mode Pilot `
    -ProjectRoot "C:\your_actual_project" `
    -AuditDirectory "C:\your_actual_project\brca_revision2_asset_audit_TIMESTAMP" `
    -Device cuda -BatchSize 32
```

The pilot reads at most 128 patches from one locked slide. It writes only a
PILOT_ONLY result, which cannot pass the complete-feature gate. Inspect the
pilot and its timing before choosing the production batch size. The original
UNI metadata recorded an RTX 4070 SUPER; the inventory checks current hardware
instead of assuming that configuration remains available.

```powershell
.\Start-BRCARevision2.ps1 -Mode Extract `
    -ProjectRoot "C:\your_actual_project" `
    -AuditDirectory "C:\your_actual_project\brca_revision2_asset_audit_TIMESTAMP" `
    -Device cuda -BatchSize 32 `
    -OutputDirectory "C:\your_actual_project\brca_revision2_hibou_complete"
```

Production reads the exact saved coordinates at level zero, resamples the
same physical patch to 256 pixels, and uses the official Hibou image processor.
Inference is float32 without autocast. The model commit, snapshot files,
processor, software, GPU, source code, coordinates and full WSI hashes are
recorded. Per-slide means use float64 accumulation and a float32 result;
patient features average slide means equally, matching the original UNI
pooling. No new patch sampling or QC threshold is introduced.

Each slide is checkpointed. To resume, use the same output directory,
environment, model snapshot/revision, precision and batch size. Changed inputs
or checkpoint contents stop the run. Full WSI hashing is performed at extraction
and resume, so a resume includes image-read overhead. No patch-level bitmap or
large embedding corpus is retained; slide means and patient features suffice
for this planned mean-pooled comparison.

```powershell
.\Start-BRCARevision2.ps1 -Mode Evaluate `
    -ProjectRoot "C:\your_actual_project" `
    -FeatureDirectory "C:\your_actual_project\brca_revision2_hibou_complete"
```

The evaluator requires all 345 x 768 real features and the original patient
index. It uses standardized balanced logistic regression, training-only
hyperparameter selection, original internal sigmoid ensembles, and original
site validation Platt/temperature calibration. Encoder feature dimensions and
official pixel normalization differ as required by the architectures; spatial
sampling and downstream protocol remain matched. Results are exploratory,
without multiplicity adjustment. Paired intervals condition on fitted models
and do not include retraining variability. Patient predictions are compared
only when identities and labels align exactly.

## External cohort: metadata before images

The candidate is non-TCGA CPTAC-BRCA, with the cBioPortal study brca_cptac_2020
as a molecular/clinical source. TCIA registry counts are not an eligible
cohort denominator. We must establish the imaging/molecular intersection and
assay provenance before choosing downloads or claiming independent validation.

```powershell
.\Start-BRCARevision2.ps1 -Mode ExternalMetadata -ProjectRoot "C:\your_actual_project"
```

This makes read-only public API requests for study, sample, attribute, clinical
and molecular-profile metadata. It does not download WSIs or assign endpoint
labels. Candidate field names are collected for review. ER/PR/HER2 must come
from compatible assay-status fields; PAM50 must be molecularly assigned;
PIK3CA wild type requires a verified mutation-profiled denominator. A documented
sample/patient/slide mapping and pretraining-overlap review are required.

After those checks, independent extraction must use a frozen external QC and
sampling protocol and the same encoder transformations. Final TCGA-trained
models and calibrators must be frozen before external predictions. External
outcomes cannot be used for tuning, recalibration or choosing an encoder.
That external extraction/evaluation stage is intentionally not represented as
ready before the cohort and schema are verified.

## Manuscript edits supported now

The manuscript_edits folder contains replacement Abstract, calibration Results,
Discussion, controlled preprocessing Results, Methods, four numerical tables
and a new calibration figure. These use completed revision-1 results only.

To apply them as a new working copy beside your CURRENT source:

```powershell
python .\code\patch_current_manuscript.py `
    --source "C:\your_manuscript\main.tex" `
    --output "C:\your_manuscript\main_revision_controls.tex"
```

The patcher requires recognizable current-source anchors. It rejects the older
September source with stale mean-fold ablation text and never overwrites the
original. Compile and visually inspect the working copy before submission.
The latest submitted source was not available here, so this package is not an
edited, submission-ready main manuscript. Encoder-comparison and external
results must be integrated only after those experiments actually finish.

## Provenance and reuse

MODEL_AND_COHORT_SOURCES.md lists primary sources. The code uses an immutable
model commit once a run starts and records hashes for all relevant snapshot
files. Unknown pretraining overlap must remain an explicit limitation.
Original code and derived data retain the accompanying licences. Raw WSIs and
UNI/Hibou weights are not redistributed.
