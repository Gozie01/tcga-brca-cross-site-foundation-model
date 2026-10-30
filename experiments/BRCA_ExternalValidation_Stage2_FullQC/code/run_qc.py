"""Outcome-blind CPTAC full-cohort QC using the supplied, unmodified TCGA patch policy."""
import argparse, csv, hashlib, json, sys, zipfile
from pathlib import Path
from datetime import datetime, timezone

def sha(p):
    return hashlib.file_digest(Path(p).open('rb'), 'sha256').hexdigest()

def read_csv(p):
    with Path(p).open(encoding='utf-8-sig', newline='') as f:
        return list(csv.DictReader(f))

def select_pilot(rows):
    groups = {}
    for r in rows:
        groups.setdefault(float(r['mpp_x']), []).append(r)
    chosen = []
    for mpp, group in sorted(groups.items()):
        patients = set()
        for r in sorted(group, key=lambda x: x['file_name']):
            if r['patient_id'] not in patients:
                chosen.append(r); patients.add(r['patient_id'])
            if len(patients) == 2:
                break
    return chosen

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--project-root', type=Path, required=True)
    ap.add_argument('--wsi-root', type=Path, required=True)
    a = ap.parse_args()
    import locked_tcga_patch as policy
    import numpy as np
    from PIL import Image, ImageDraw
    source = a.project_root / 'BRCA_ExternalValidation_Stage1/data/download_manifest.csv'
    rows = read_csv(source)
    if len(rows) != 391 or len({r['file_name'] for r in rows}) != 391:
        raise RuntimeError('Expected 391 unique Stage1 slides')
    # Check the previously frozen acquisition manifest before using it.
    protocol = json.loads((source.parents[1] / 'acquisition_protocol.json').read_text())
    expected = protocol['locks']['data/download_manifest.csv']
    if sha(source) != expected:
        raise RuntimeError('Stage1 download manifest checksum mismatch')
    selected = sorted(rows, key=lambda r: r["file_name"])
    if len(selected) != 391 or len({r['patient_id'] for r in selected}) != 120:
        raise RuntimeError('Expected 391 slides from 120 patients')
    index = {}
    for p in a.wsi_root.rglob('*.svs'):
        if p.name in index:
            raise RuntimeError('Duplicate WSI basename: ' + p.name)
        index[p.name] = p
    for r in selected:
        p = index.get(r['file_name'])
        if p is None or p.stat().st_size != int(r['bytes']):
            raise RuntimeError('Missing or size-mismatched slide: ' + r['file_name'])
    stamp = datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
    out = a.project_root / 'BRCA_ExternalValidation_Stage2_FullQC' / 'full_qc_v1'
    out.mkdir(parents=True, exist_ok=True)
    lock = dict(manifest_sha256=sha(source), policy_sha256=sha(Path(policy.__file__)), runner_sha256=sha(__file__), wsi_root=str(a.wsi_root.resolve()))
    lock_path = out/'run_lock.json'
    if lock_path.exists() and json.loads(lock_path.read_text()) != lock:
        raise RuntimeError('Run inputs/code changed; do not mix checkpoints. Use a separate output folder.')
    policy.atomic_json(lock_path, lock)
    policy.atomic_json(out/'run_status.json', {'status':'RUNNING_REVIEW_PENDING'})
    (out / 'selection.json').write_text(json.dumps(selected, indent=2))
    # Never parse CPTAC patient identifiers as TCGA identifiers.
    patient_map = {r['file_name']: r['patient_id'] for r in selected}
    policy.patient_id = lambda name: patient_map[name]
    policy.tss = lambda name: 'CPTAC'
    records = []
    retained_total = 0
    retained_patients = set()
    combined_tmp = out/'retained_coordinates.csv.tmp'
    combined_file = combined_tmp.open('w',encoding='utf-8-sig',newline='')
    combined_writer = csv.DictWriter(combined_file,fieldnames=policy.FIELDS)
    combined_writer.writeheader()
    for number, r in enumerate(selected, 1):
        p = index[r['file_name']]
        with policy.openslide.OpenSlide(str(p)) as slide:
            mx = float(slide.properties[policy.openslide.PROPERTY_NAME_MPP_X])
            my = float(slide.properties[policy.openslide.PROPERTY_NAME_MPP_Y])
            if abs(mx-float(r['mpp_x'])) > 1e-6 or abs(my-float(r['mpp_y'])) > 1e-6:
                raise RuntimeError('MPP mismatch: '+p.name)
            if abs(mx-my) > 1e-6:
                raise RuntimeError('Anisotropic MPP needs separate protocol review')
        meta = policy.process_slide('CPTAC', number, len(selected), p, out)
        if meta['status'] != 'COMPLETE':
            raise RuntimeError(meta.get('error', 'Patch extraction failed'))
        patches = read_csv(meta['checkpoint_csv'])
        kept = [x for x in patches if x['retain_after_new_locked_qc'].lower() == 'true']
        retained_total += len(kept)
        if kept: retained_patients.add(r['patient_id'])
        if any(x['patient_id'] != r['patient_id'] or x['file_name'] != r['file_name'] for x in patches):
            raise RuntimeError('Checkpoint identity mismatch')
        if len({(x['x_level0'],x['y_level0']) for x in kept}) != len(kept):
            raise RuntimeError('Duplicate retained coordinates')
        for x in kept:
            if not (float(x['tissue_fraction_hr']) >= 0.35-1e-6 and float(x['white_fraction']) <= 0.60+1e-6 and float(x['artifact_fraction']) <= 0.01+1e-6 and float(x['bg_to_red_ratio']) < 1.10 and float(x['mean_local_std']) >= 0.50):
                raise RuntimeError('Retained patch violates locked QC')
        combined_writer.writerows(kept)
        combined_file.flush()
        canvas = Image.new('RGB', (1024, 1120), 'white')
        draw = ImageDraw.Draw(canvas)
        edge_count = 0
        with policy.openslide.OpenSlide(str(p)) as slide:
            width, height = slide.dimensions
            for x in kept:
                xx, yy, size = (int(x[k]) for k in ('x_level0','y_level0','read_size_level0'))
                edge_count += int(xx < 0 or yy < 0 or xx+size > width or yy+size > height)
            # Evenly spaced through retained coordinate order, independent of labels/scores.
            preview = np.linspace(0, len(kept)-1, min(16,len(kept)), dtype=int) if kept else []
            for i, j in enumerate(preview):
                x = kept[int(j)]; size = int(x['read_size_level0'])
                patch = slide.read_region((int(x['x_level0']),int(x['y_level0'])),0,(size,size)).convert('RGB')
                patch = patch.resize((256,256), Image.Resampling.LANCZOS)
                left, top = (i%4)*256, (i//4)*280
                canvas.paste(patch,(left,top))
                draw.text((left+3,top+258), 'patch '+x['patch_id'], fill='black')
        canvas.save(out / ('preview_'+str(number)+'.jpg'), quality=92)
        meta.update(specimen_link_method=r['specimen_link_method'], zero_retained_flag=not bool(kept), source_mpp_x=mx, source_mpp_y=my, retained_edge_crossing_patches=edge_count,
                    visual_tumor_qc='PENDING_HUMAN_REVIEW')
        records.append(meta)
        print(f'FULL QC {number}/391: retained={len(kept)}, edge-crossing={edge_count}', flush=True)
    combined_file.close()
    combined_tmp.replace(out/'retained_coordinates.csv')
    summary = dict(status='FULL_TECHNICAL_QC_COMPLETE_REVIEW_PENDING', slides=records,
        manifest_sha256=sha(source), locked_policy_sha256=sha(Path(policy.__file__)),
        retained_patches=retained_total, retained_patients=len(retained_patients), zero_retained_slides=sum(x['zero_retained_flag'] for x in records), coordinates_sha256=sha(out/'retained_coordinates.csv'), target_mpp=policy.TARGET_MPP, output_size=policy.PATCH_OUT,
        max_initial_accepted_per_slide=policy.MAX_PATCHES_PER_SLIDE, seed=policy.RANDOM_SEED,
        selection='All 391 frozen Stage1 slides sorted by filename; no label/model-score selection',
        edge_handling='Preserved Stage9A/9B full read_region request -> RGB -> 256 Lanczos; crossing audited, no clipping/exclusion',
        human_visual_qc_complete=False, tumor_content_confirmed=False,
        feature_extraction_run=False, external_predictions_computed=False,
        python=sys.version, numpy=np.__version__)
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    policy.atomic_json(out/'run_status.json', {'status':summary['status']})
    report = out.parent / ('brca_stage2_full_qc_'+stamp+'.zip')
    with zipfile.ZipFile(report,'w',zipfile.ZIP_DEFLATED) as z:
        for p in sorted(out.rglob('*')):
            if p.is_file():z.write(p,p.relative_to(out))
    print('REPORT ZIP: '+str(report),flush=True)

if __name__ == '__main__':
    main()
