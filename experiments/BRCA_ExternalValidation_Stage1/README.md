# Start the CPTAC technical pilot

The matching audit is complete. Read `STAGE1_REPORT.md` for the cohort, reference definitions and limits.

1. Extract this ZIP into `C:\Users\Gozie\Desktop\BioMedical`.
2. Open PowerShell in that project and run:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass `
    -File ".\BRCA_ExternalValidation_Stage1\Start-BRCAExternalStage1.ps1" `
    -Mode Pilot `
    -ProjectRoot "C:\Users\Gozie\Desktop\BioMedical" `
    -WSIRoot "D:\CPTAC_BRCA_WSI"
```

3. Upload the report ZIP printed at the end. It contains decode checks and thumbnails, not whole-slide images.

This starts only the three-slide technical pilot. It does not run external feature extraction or predictions. Rerun the same command to resume an interrupted transfer.

The existing BRCA environment must be available at `.venv-brca-revision`. The downloader uses Python's standard library; slide checks use the already installed OpenSlide and Pillow. No package installation or model login is performed.

For integrity and space checks alone, use `-Mode Audit`. Full acquisition is available with `-Mode DownloadAll` after the pilot is reviewed. Detailed source and reproducibility instructions are in `REPRODUCIBILITY.md`.
