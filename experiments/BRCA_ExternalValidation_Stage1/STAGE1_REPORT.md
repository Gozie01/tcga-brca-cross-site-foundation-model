# CPTAC external validation: stage 1

The image/specimen linkage and remote file preflight are complete. This package starts a technical download pilot; it does not compute external predictions or performance.

## Acquisition result

- Molecular study: cBioPortal `brca_cptac_2020`, 122 patients / primary-tumor samples.
- Patients with matching image identifiers: 120.
- Ready patient-level tumor-image links: 120 patients, 391 slides, 67.91 GB (63.25 GiB).
- Specimen links: {"exact": 298, "unique_abbreviated_specimen_alias": 93}.
- Molecular records without image matches: CPT000814, CPT001846.
- Technical pilot: 3 slides, 514.5 MB.

All selected file links responded successfully to a bounded HTTPS range request and returned TIFF headers. File sizes, modification times and the SHA-256 hashes of the first 4,096 bytes were recorded. Original full-file checksums are not published in the retrieved records; the downloader records full-file SHA-256 after transfer, checks source identity/size/range behaviour, and requires an OpenSlide decode check. This is an acquisition integrity audit, not independent validation of every source pixel.

## Matching and exclusions

The cBioPortal identifiers of form `XNNBRNNN` map to `NNBRNNN`; other supported identifiers remain unchanged. The mapping is one-to-one and is supported by the matching CPTAC cohort records. Age concordance is assessed at completed-year precision because one source supplies fractional years: {"AGREE": 27, "MISSING_ONE_OR_BOTH": 18, "AGREE_COMPLETED_YEARS": 77}. Sex comparison: {"AGREE": 105, "MISSING_ONE_OR_BOTH": 17}.

The public cohort-builder export contains 654 rows and 199 distinct subject strings. These inventory counts differ from the historical collection landing-page statistics. We do not use the landing-page counts as the eligible cohort. Each selected row is independently checked against its current individual PathDB image record and source filename. 2 exported node identifiers were corrected through unique exact patient-and-image matches to the current records; the originals and corrected identifiers are both preserved.

Tumor status comes from the archived ESAC specimen records preserved by the IDC conversion project. Full stripped slide labels are matched exactly when available. Abbreviated identifiers such as `[fe]-c3b8ed_D1_D1` are matched only when their first two and last six hexadecimal characters, section suffix and patient agree with exactly one image and one specimen. Abbreviated links are separately identified, allowing a strict-link sensitivity analysis. Repeated nested lists of all slides in a case are not accepted as proof of the individual specimen identity.

Images without a molecular record, without a unique specimen link, or without documented tumor-tissue status are excluded. Exclusion counts: {"NO_MATCHED_MOLECULAR_RECORD": 259, "NO_UNIQUE_SPECIMEN_LINK": 4}. `excluded_images.csv` retains all reasons. Selection uses metadata and acquisition integrity, never model predictions or endpoint performance. All qualifying slides are retained for the future hierarchical pooling protocol; no best-performing slide is selected.

Linkage is at patient level. The image has not been proven to be from the exact block used for molecular assays. The archived specimen records also do not reliably identify frozen versus FFPE preparation, so equivalence to the TCGA diagnostic preparation is not established. These limitations must be retained in the external-validation report. The current stage documents tumor annotations and decode readiness; visual quality review and image-based exclusions are still pending.

## Endpoint reference definitions

| Endpoint | Reference + image linkage, before image QC | Class counts |
| --- | --- | --- |
| ER_binary | 118 | {'0': 37, '1': 81} |
| PR_binary | 112 | {'0': 45, '1': 67} |
| HER2_binary | 95 | {'0': 82, '1': 13} |
| PAM50_primary_4class | 114 | {'2': 27, '3': 14, '0': 56, '1': 17} |
| PIK3CA_binary | 120 | {'0': 81, '1': 39} |

The counts are before visual image QC and are not final evaluated sample sizes. Binary labels are 0 = negative and 1 = positive. PAM50 uses the original TCGA model order: 0 = LumA, 1 = LumB, 2 = Basal, 3 = Her2. Normal-like and missing PAM50 are excluded. Missing and equivocal HER2 are excluded; gene amplification, proteomic status and PAM50 are not substituted for clinical HER2.

PIK3CA positivity requires a reported mutation in the original TCGA qualifying variant-class set. The cBioPortal sequenced list, sample `sequenced` flags and WXS profile matrix agree on all 122 samples. A negative reference means no qualifying reported mutation in this profiled denominator, not guaranteed base-level callability or clinically confirmed wild type. The original TCGA class set was checked directly against the supplied v2 labeling script. Source somatic filtering differs across cohorts and remains a limitation.

The endpoint definitions are frozen before external feature extraction or prediction. The future TCGA-trained models and calibrators must also be frozen before inspecting external outcomes. No external-outcome tuning, model selection or recalibration is allowed in the primary analysis.

## Pretraining documentation

For original UNI, the authors report excluding public TCGA, CPTAC and PAIP slide collections from pretraining (https://doi.org/10.1038/s41591-024-02857-3). The Hibou paper describes pretraining exclusively on a private dataset (https://arxiv.org/html/2406.05074v1). These are author-reported provenance statements, not an independent patient-level overlap audit. A CPTAC-specific training manifest for the precise Hibou checkpoint is not available here; do not claim proven zero overlap at the participant level.

## Run the pilot

Use `Start-BRCAExternalStage1.ps1 -Mode Pilot` with your existing BioMedical project and `D:\CPTAC_BRCA_WSI` as the new image directory. The script uses the existing `.venv-brca-revision` environment to decode the downloaded slides. No new Hugging Face login is needed for this acquisition stage.

| Patient | Pilot slide | Size (MB) |
| --- | --- | --- |
| 01BR001 | 0684a407-f446-486d-9160-b483cb | 253.7 |
| 01BR008 | 32269c0c-3ca2-4de1-b314-e9b07c | 90.9 |
| 01BR009 | 0ecabbfa-8095-4875-acdb-051582 | 170.0 |

The pilot selects one slide from each of the first three patient identifiers, without inspecting labels. It verifies downloading, TIFF/WSI decoding, dimensions and physical pixel spacing, and generates thumbnails for review. It does not estimate performance. Partial files can be resumed by rerunning the same command. Existing completed files are checked against their receipts rather than overwritten.

At completion, upload the report ZIP printed by the script. It includes acquisition receipts, decode checks and thumbnails, but no full slide images or credentials. If a transfer is interrupted, the script preserves partial files and writes an interruption report where possible.

`-Mode Audit` only checks package integrity and required disk space for the pilot. `-Mode DownloadAll` downloads the full locked inventory; it is provided for the later acquisition step. Start with the pilot so its decoder checks and thumbnails can be reviewed before the full transfer. The original TCGA images, features and results are not changed.

## Subsequent work

After a successful pilot: review slide quality; complete acquisition; freeze an outcome-independent patch/QC protocol consistent with the TCGA protocol; extract UNI and Hibou-B on the same external patches; finalize and checksum the TCGA-trained prediction objects; and evaluate each endpoint once with patient-level uncertainty. Same-patient links, abbreviated-specimen sensitivity, source assay differences and pretraining uncertainty must accompany the final results. No external accuracy or AUROC is available at stage 1.

Sources: https://www.cancerimagingarchive.net/collection/cptac-brca/ ; https://pathdb.cancerimagingarchive.net/ ; https://github.com/ImagingDataCommons/idc-wsi-conversion ; https://www.cbioportal.org/study?id=brca_cptac_2020 . Dataset citation: National Cancer Institute CPTAC (2020), CPTAC-BRCA Version 1, TCIA, DOI 10.7937/TCIA.CAEM-YS80, CC BY 4.0. Full source URLs, recorded hashes, public metadata snapshots and current per-image records are included in the evidence directory.
