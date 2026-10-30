"""Check actual SVS decoding and pixel metadata without molecular prediction."""
import argparse
import csv
import hashlib
import json
from pathlib import Path


def check(manifest, wsi_root, out):
    import openslide
    from PIL import ImageStat
    rows=list(csv.DictReader(manifest.open(encoding='utf-8')))
    records=[]
    failures=[]
    thumbdir=out.parent/'qc_thumbnails';thumbdir.mkdir(exist_ok=True)
    for r in rows:
        path=wsi_root/r['patient_id']/r['file_name']
        item={'patient_id':r['patient_id'],'file_name':r['file_name']}
        try:
            if path.stat().st_size!=int(r['bytes']):raise ValueError('File size mismatch')
            with openslide.OpenSlide(str(path)) as slide:
                width,height=slide.dimensions
                if width!=int(r['width']) or height!=int(r['height']):raise ValueError('Decoded dimensions differ from source record')
                mpp=[float(slide.properties.get(k,0)) for k in ['openslide.mpp-x','openslide.mpp-y']]
                if min(mpp)<=0:raise ValueError('Missing or invalid slide pixel spacing')
                for got,wanted in zip(mpp,[r['mpp_x'],r['mpp_y']]):
                    if wanted and abs(got-float(wanted))>0.001:raise ValueError('Pixel spacing differs from remote metadata')
                thumb=slide.get_thumbnail((1200,1200)).convert('RGB')
                thumb.save(thumbdir/(r['pathdb_node_id']+'.jpg'),quality=88)
                # Decode tiles at several positions; blank regions are allowed, decode failures are not.
                tile_hashes=[]
                for fx,fy in [(.25,.25),(.5,.5),(.75,.75)]:
                    tile=slide.read_region((min(int(width*fx),width-256),min(int(height*fy),height-256)),0,(256,256)).convert('RGB')
                    tile_hashes.append(hashlib.sha256(tile.tobytes()).hexdigest())
                item.update({'status':'READABLE','width':width,'height':height,'mpp_x':mpp[0],'mpp_y':mpp[1],
                             'vendor':slide.properties.get('openslide.vendor',''),'level_count':slide.level_count,
                             'decoded_tile_hashes':tile_hashes,'thumbnail_RGB_stddev':ImageStat.Stat(thumb).stddev,
                             'thumbnail':str((thumbdir/(r['pathdb_node_id']+'.jpg')).relative_to(out.parent)),
                             'visual_tumor_QC':'PENDING_HUMAN_REVIEW','feature_extraction_run':False})
        except Exception as e:
            item.update({'status':'FAILED','error':str(e)});failures.append(item)
        records.append(item);print(item['status'],r['file_name'],flush=True)
    out.write_text(json.dumps({'status':'DECODE_CHECKS_PASSED' if not failures else 'DECODE_CHECKS_FAILED',
        'slides':len(rows),'failures':len(failures),'records':records,'external_predictions_computed':False},indent=2))
    if failures:raise ValueError('Some slides failed decoding; see report. No feature extraction was run.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--wsi-root',type=Path,required=True);p.add_argument('--out',type=Path,required=True)
    a=p.parse_args();check(a.manifest,a.wsi_root,a.out)
