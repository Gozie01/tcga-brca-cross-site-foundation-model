# Existing-slide stain provenance check

The GDC API matched the exact TCGA-5L-AAT0 slide filename to file b3780253-78fe-4907-9f5a-230b2bb4b24e, classified as Diagnostic Slide / Slide Image, 349684623 bytes, MD5 b03df56d20c11c51b8ad6c3e132a78de. Expanded slide metadata did not provide a stain field. Therefore stain type remains unconfirmed; unusual color is not sufficient to assign IHC or exclude the slide.

Run inside this extracted folder:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File .\Start-StainProvenance.ps1
```

The default search is D:\TCGA_BRCA_WSI. If necessary, specify the existing file:

```powershell
py -3.11 .\check_stain_provenance.py --slide "FULL PATH TO THE EXISTING SVS"
```

Upload the generated brca_stain_provenance_*.zip. It contains slide properties, available label/macro images, a thumbnail and a byte-identity check. No new slide download, patch extraction, model training, source mutation or eligibility change is performed. Requires the Python 3.11 OpenSlide/Pillow environment already used for reference collection. Syntax was checked here; execution against the Windows slide awaits the local run.

Existing frozen results remain primary. Revised QC, if later justified, is post-hoc exploratory. Tumor presence and content are not established by this technical check. Validation panels remain reserved.
