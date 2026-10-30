"""Read-only provenance collection for one already downloaded TCGA slide."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import zipfile

NAME = 'TCGA-5L-AAT0-01Z-00-DX1.5E171263-30BF-4C6B-88A1-E8EA0522A861.svs'
GDC_ID = 'b3780253-78fe-4907-9f5a-230b2bb4b24e'
EXPECTED_MD5 = 'b03df56d20c11c51b8ad6c3e132a78de'
EXPECTED_SIZE = 349684623

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--slide', type=Path)
    p.add_argument('--wsi-root', type=Path, default=Path('D:/TCGA_BRCA_WSI'))
    a = p.parse_args()
    if a.slide:
        path = a.slide
    else:
        matches = list(a.wsi_root.rglob(NAME)) if a.wsi_root.exists() else []
        if len(matches) != 1:
            p.error('Expected exactly one slide. Supply --slide "FULL PATH TO SVS". Found: ' + str(len(matches)))
        path = matches[0]
    if not path.is_file():
        p.error('Slide file does not exist: ' + str(path))
    import openslide
    out = Path.cwd() / ('brca_stain_provenance_' + datetime.datetime.now().strftime('%Y%m%d_%H%M%S'))
    out.mkdir(exist_ok=False)
    print('Checking existing slide identity; no download or model training.', flush=True)
    digest = hashlib.md5()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(8 * 1024 * 1024), b''):
            digest.update(block)
    record = {'source_path': str(path.resolve()), 'file_name': path.name,
              'file_size': path.stat().st_size, 'md5': digest.hexdigest(),
              'gdc_file_id': GDC_ID, 'expected_md5': EXPECTED_MD5,
              'expected_size': EXPECTED_SIZE,
              'stain_status': 'NOT_ESTABLISHED', 'models_or_source_changed': False}
    record['gdc_bytes_match'] = record['md5'] == EXPECTED_MD5 and record['file_size'] == EXPECTED_SIZE
    with openslide.OpenSlide(str(path)) as slide:
        record['dimensions'] = slide.dimensions
        record['properties'] = dict(slide.properties)
        record['associated_image_names'] = list(slide.associated_images)
        slide.get_thumbnail((2000, 2000)).convert('RGB').save(out / 'whole_slide_thumbnail.jpg', quality=95)
        for key in ('label', 'macro'):
            if key in slide.associated_images:
                im = slide.associated_images[key].convert('RGB')
                im.thumbnail((2400, 2400))
                im.save(out / (key + '.jpg'), quality=95)
    (out / 'provenance.json').write_text(json.dumps(record, indent=2), encoding='utf-8')
    (out / 'README.txt').write_text('Technical source identity and stain-provenance evidence only. No diagnostic review or eligibility change. Stain type is not inferred automatically. GDC diagnostic-slide classification does not by itself establish the stain.\n', encoding='utf-8')
    report = out.with_suffix('.zip')
    with zipfile.ZipFile(report, 'w', zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):
            z.write(f, f.name)
    print('GDC byte identity match:', record['gdc_bytes_match'])
    print('UPLOAD REPORT:', report.resolve())

if __name__ == '__main__':
    main()
