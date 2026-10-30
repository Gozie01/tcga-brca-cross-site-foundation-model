# Reproducing the stage-1 linkage

The pilot launcher does not require these steps. All source evidence and the locked manifests are already supplied.

`code/prepare_inventory.py` uses the cohort-builder subset, the archived CPTAC2/BR specimen records and the cBioPortal sample/clinical records. It saves the current per-image JSON and bounded remote file probes, verifies patient and image identities, and creates the download manifest. Exact and uniquely abbreviated specimen links are distinguished. The public export contains a pair of swapped node identifiers; the code resolves this only through a unique already observed record with the exact patient and image, recording both old and new node identifiers.

The original 193 MB cohort-builder file contains unrelated collections and is not duplicated in this package. Its hash and extraction rule are recorded in `SOURCE_PROVENANCE.json`; the exact CPTAC-BRCA subset is supplied. The original specimen ZIP hash and the extracted member name are also recorded. `code/refresh_sources.py` can fetch a new raw snapshot into a separate directory. It contacts public sources and does not require a login. A refreshed snapshot must be reviewed as a new version rather than replacing the frozen acquisition silently.

To rebuild the existing linkage from the supplied source subset, with Python 3.11+:

```powershell
py -3.11 .\code\prepare_inventory.py `
    --sources .\evidence\source_metadata `
    --clinical .\evidence\clinical `
    --out .\Reaudit `
    --workers 4
```

This performs network checks for the slide records and first 4,096 bytes of each candidate file. It does not download whole slides. To reuse cached checks, use a separate copy of the delivered package as `--out`, preserving its `evidence/pathdb_nodes` and `evidence/file_probes` directories.

`code/freeze_labels.py` validates the PIK3CA variant-class set against the supplied original TCGA v2 source and validates PAM50 class order against the original frozen TCGA partitions. It reconciles the WXS matrix, sequenced list and API sample flags before assigning a mutation-negative reference. To regenerate it, use the supplied `evidence/tcga_reference` directory as the TCGA pipeline and its `pik3ca_audit_v2.py` as the TCGA definition. Use the source and clinical directories above and an output directory already containing the patient crosswalk.

```powershell
py -3.11 .\code\freeze_labels.py `
    --source .\evidence\source_metadata `
    --clinical .\evidence\clinical `
    --tcga-definition .\evidence\tcga_reference\pik3ca_audit_v2.py `
    --tcga-pipeline .\evidence\tcga_reference `
    --out .\Reaudit
```

The code is standard-library Python except for `check_slides.py`, which uses OpenSlide and Pillow in the existing BRCA environment. Local software tests run with:

```powershell
py -3.11 .\code\test_download.py
```

Tests exercise interrupted-transfer resume, incorrect HTTP ranges, source changes, local corruption, partial identity mismatches, existing-file preservation, completed-part recovery and conservative identifier/section matching. The HTTP transport is simulated; these tests do not require network access or download slide images. Actual source access was separately tested using the saved bounded probes.

The entire package is protected by `PACKAGE_SHA256.json`, checked before either Audit or download mode. A receipt records each downloaded file's full SHA-256. A local receipt is not an independently published source checksum; the remote records did not supply one. Original data and pipeline outputs are preserved.
