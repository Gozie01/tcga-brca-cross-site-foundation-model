"""Read exact frozen WSI coordinates; descriptive diagnostics only, no QC exclusions."""
from pathlib import Path
import argparse,csv,json,hashlib,zipfile,datetime,sys,platform,os
from PIL import Image,ImageDraw
import numpy as np
R=Path(__file__).resolve().parent

def source_patch(slide,row):
 return slide.read_region((int(row['x_level0']),int(row['y_level0'])),0,(int(row['read_size_level0']),int(row['read_size_level0']))).convert('RGB').resize((int(row['output_size']),int(row['output_size'])),Image.Resampling.LANCZOS)

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--all-flagged',action='store_true');ap.add_argument('--wsi-root',type=Path);ap.add_argument('--qc-zip',type=Path);ap.add_argument('--out',type=Path);args=ap.parse_args()
 try:import openslide
 except ImportError:raise SystemExit('OpenSlide is unavailable in this Python environment. Run this script using the same environment that successfully extracted the BRCA features. No slides have been changed.')
 stamp=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');out=args.out or R/('trace_'+stamp);out.mkdir(parents=True,exist_ok=False)
 rows=list(csv.DictReader((R/'patch_trace_manifest.csv').open(newline='',encoding='utf-8-sig')))
 if not args.all_flagged:rows=[x for x in rows if x['priority']=='FIRST_TECHNICAL_INSPECTION']
 candidates=[]
 if args.qc_zip:candidates=[args.qc_zip]
 else:
  candidates+=list(R.glob('brca_stage2_full_qc_20261021_041052_426646*.zip'))
  base=Path.home()/'Desktop'/'BioMedical';folder=base/'BRCA_ExternalValidation_Stage2_FullQC'
  if folder.exists():candidates+=list(folder.rglob('brca_stage2_full_qc_20261021_041052_426646*.zip'))
  downloads=Path.home()/'Downloads'
  if downloads.exists():candidates+=list(downloads.glob('brca_stage2_full_qc_20261021_041052_426646*.zip'))
 if not candidates:raise SystemExit('QC report ZIP not found. Copy brca_stage2_full_qc_20261021_041052_426646.zip beside this script, or supply --qc-zip "FULL PATH".')
 qc=candidates[0];expected=json.loads((R/'trace_plan.json').read_text());actual=hashlib.sha256(qc.read_bytes()).hexdigest()
 if actual!=expected['qc_zip_sha256']:raise SystemExit('QC ZIP hash differs from the traced evidence. No WSI read started; use the exact original report ZIP.')
 z=zipfile.ZipFile(qc);results=[];errors=[];sizes={}
 try:
  for rid in sorted({x['review_id'] for x in rows}):
   group=[x for x in rows if x['review_id']==rid];first=group[0];src=Path(first['source_path'])
   if args.wsi_root:src=args.wsi_root/first['patient_id']/first['file_name']
   print(rid,'reading',first['file_name'],flush=True)
   try:
    if not src.is_file():raise FileNotFoundError(str(src))
    if src.stat().st_size!=int(first['expected_source_bytes']):raise ValueError('WSI file size differs from the original QC record')
    prior=Image.open(z.open(first['preview_name'])).convert('RGB');assert prior.size==(1024,1120)
    panel=Image.new('RGB',(2048,1120),'white');draw=ImageDraw.Draw(panel)
    with openslide.OpenSlide(str(src)) as slide:
     sizes[rid]={'dimensions':slide.dimensions,'source_bytes':src.stat().st_size,'mtime_ns':src.stat().st_mtime_ns,'path':str(src)}
     for row in group:
      j=int(row['preview_slot']);cx=(j%4)*256;cy=(j//4)*280;ref=prior.crop((cx,cy,cx+256,cy+256));patch=source_patch(slide,row)
      a=np.asarray(ref,dtype=float);b=np.asarray(patch,dtype=float);gray=b.mean(2);sat=b.max(2)-b.min(2)
      # Gradient is descriptive only; no validated focus threshold is claimed.
      gradient=float((np.abs(np.diff(gray,axis=0)).mean()+np.abs(np.diff(gray,axis=1)).mean())/2)
      results.append({**row,'preview_vs_reread_mean_abs_rgb':float(np.abs(a-b).mean()),'preview_vs_reread_p95_abs_rgb':float(np.percentile(np.abs(a-b),95)),'reread_mean_rgb':float(b.mean()),'reread_mean_channel_range':float(sat.mean()),'reread_gray_std':float(gray.std()),'reread_mean_absolute_gradient':gradient,'source_inside_bounds':int(row['x_level0'])>=0 and int(row['y_level0'])>=0 and int(row['x_level0'])+int(row['read_size_level0'])<=slide.dimensions[0] and int(row['y_level0'])+int(row['read_size_level0'])<=slide.dimensions[1]})
      x=(j%4)*512;y=(j//4)*280;panel.paste(ref,(x,y));panel.paste(patch,(x+256,y));draw.text((x+4,y+258),f'{rid} patch {row["patch_id"]} | preview / WSI reread',fill='black')
    panel.save(out/(rid+'_comparison.png'))
   except Exception as e:errors.append({'review_id':rid,'error':str(e)})
 finally:z.close()
 if results:
  with (out/'patch_comparison_metrics.csv').open('w',newline='',encoding='utf-8-sig') as f:
   w=csv.DictWriter(f,fieldnames=list(results[0]));w.writeheader();w.writerows(results)
 summary={'status':'COMPLETE' if not errors else 'INCOMPLETE_WITH_ERRORS','slides_requested':len({x['review_id'] for x in rows}),'patches_requested':len(rows),'patches_reread':len(results),'errors':errors,'source_files':sizes,'qc_zip_sha256':actual,'python':platform.python_version(),'openslide_python':openslide.__version__,'numpy':np.__version__,'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),'source_WSI_hash_recomputed':False,'features_or_models_changed':False,'diagnosis_or_eligibility_assigned':False,'interpretation':'Preview/WSI pixel differences include contact-sheet JPEG compression. Inspect paired images; descriptive metrics are not validated exclusion rules. This task does not establish tumor content.'}
 (out/'summary.json').write_text(json.dumps(summary,indent=2));archive=out.with_suffix('.zip')
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as zz:
  for p in out.iterdir():
   if p.is_file():zz.write(p,p.name)
 print('UPLOAD THIS REPORT:',archive,flush=True)
 if errors:raise SystemExit(1)
if __name__=='__main__':main()
