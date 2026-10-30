"""Collect outcome-free QC reference patches from frozen TCGA development patients.
Does not filter patches, extract model features, fit models or read outcome labels.
"""
from pathlib import Path
import argparse,hashlib,json,datetime,zipfile,platform
import numpy as np,pandas as pd
from PIL import Image,ImageDraw
from qc_measurements import measure
R=Path(__file__).resolve().parent
EXPECTED_COORD='eb5591a68242ff6c8a97759cdf51e9bc70df050eda7c1cd0115f31be044c5e66'
EXPECTED_SPLIT='3e4e2a23d8bbd1dc2989d466a77ce8e79c75ea24dc7a0a9f57df1a26395222b1'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def find_verified(name,explicit,expected):
 if explicit:
  p=explicit.resolve()
  if not p.is_file() or sha(p)!=expected:raise ValueError('Specified input missing or differs from frozen hash: '+str(p))
  return p
 base=Path.home()/'Desktop'/'BioMedical';roots=[R,base/'BRCA_Revision2_Pipeline',base/'TCGA_BRCA_Mutation_Audit',base/'TCGA_BRCA_Mutations'];seen=set()
 for root in roots:
  if not root.exists():continue
  candidates=[root/name,root/'data'/name]
  if root!=R:candidates+=list(root.rglob(name))
  for p in candidates:
   if p in seen:continue
   seen.add(p)
   if p.is_file() and sha(p)==expected:return p
 raise FileNotFoundError('Cannot find frozen '+name+'. Supply its full path with --coordinates or --splits. No WSI reads have started.')

def main():
 ap=argparse.ArgumentParser();ap.add_argument('--coordinates',type=Path);ap.add_argument('--splits',type=Path);ap.add_argument('--wsi-root',type=Path,default=Path(r'D:\TCGA_BRCA_WSI'));args=ap.parse_args()
 try:import openslide
 except ImportError:raise SystemExit('Use the same Python environment as the successful BRCA trace; OpenSlide is missing here.')
 coord=find_verified('patch_coordinates.csv.gz',args.coordinates,EXPECTED_COORD);split=find_verified('stage10a_internal_5fold_manifest_v1.csv',args.splits,EXPECTED_SPLIT)
 # Deliberately load only IDs and partition numbers, never molecular labels.
 folds=pd.read_csv(split,usecols=['patient_id','internal_fold']);held=set(folds.loc[folds.internal_fold.eq(0),'patient_id']);eligible=sorted(set(folds.patient_id)-held)
 if len(eligible)<40:raise ValueError(f'Only {len(eligible)} patients outside all endpoint fold-0 tests; 40 required. Do not silently change this protocol.')
 rng=np.random.default_rng(20261022);rng.shuffle(eligible);chosen=eligible[:40]
 table=pd.read_csv(coord);required=['patient_id','file_name','x_level0','y_level0','read_size_level0']
 if not set(required)<=set(table):raise ValueError('Coordinate schema missing required fields: '+str(set(required)-set(table)))
 table=table[table.patient_id.isin(chosen)].copy();assert not set(table.patient_id)&held
 selected=[]
 for i,pid in enumerate(chosen):
  g=table[table.patient_id.eq(pid)].sort_values(['file_name','x_level0','y_level0']).reset_index(drop=True)
  if len(g)<8:raise ValueError('Insufficient retained patches: '+pid)
  g=g.iloc[rng.choice(len(g),8,replace=False)].copy();g['qc_case_id']=f'QC{i+1:03d}';g['qc_partition']='development' if i<20 else 'validation';g['panel_slot']=range(8);selected.append(g)
 selected=pd.concat(selected,ignore_index=True)
 # Resolve every WSI before starting reads. Existing absolute paths take priority.
 pathfield=next((k for k in ['source_path','path','slide_path','wsi_path'] if k in selected),None);paths={};unresolved=[]
 for fname,g in selected.groupby('file_name'):
  p=Path(str(g.iloc[0][pathfield])) if pathfield else None
  if p and p.is_file():paths[fname]=p
  else:unresolved.append(fname)
 if unresolved:
  inventory={}
  if args.wsi_root.exists():
   for p in args.wsi_root.rglob('*.svs'):
    if p.name in unresolved:inventory.setdefault(p.name,[]).append(p)
  for name in unresolved:
   matches=inventory.get(name,[])
   if len(matches)!=1:raise FileNotFoundError('Missing or ambiguous WSI '+name+'. Supply --wsi-root with the original TCGA slide root. No WSI reads have started.')
   paths[name]=matches[0]
 now=datetime.datetime.now().strftime('%Y%m%d_%H%M%S');out=R/('tcga_qc_reference_'+now);out.mkdir(exist_ok=False)
 plan={'status':'FROZEN_BEFORE_WSI_READS','seed':20261022,'coordinate_sha256':sha(coord),'split_sha256':sha(split),'case_count':40,'patch_count':320,'patient_partition':'20 development / 20 untouched technical validation; patient-disjoint; all chosen patients excluded from every endpoint fold-0 test in frozen manifest','development_ids':chosen[:20],'validation_ids':chosen[20:],'molecular_label_columns_loaded':False,'external_data_used_for_reference_sampling':False,'thresholds_selected':False,'script_sha256':sha(Path(__file__)),'measurement_code_sha256':sha(R/'qc_measurements.py')}
 (out/'collection_protocol.json').write_text(json.dumps(plan,indent=2));selected.to_csv(out/'selected_coordinates_PRIVATE.csv',index=False)
 metrics=[];images={};errors=[]
 for fname,g in selected.groupby('file_name',sort=True):
  print('READ',fname,flush=True)
  try:
   with openslide.OpenSlide(str(paths[fname])) as slide:
    for _,row in g.iterrows():
     size=int(row.read_size_level0);x=int(row.x_level0);y=int(row.y_level0)
     if x<0 or y<0 or x+size>slide.dimensions[0] or y+size>slide.dimensions[1]:raise ValueError('Edge-crossing sample; stop and review frozen boundary behavior rather than clipping.')
     patch=slide.read_region((x,y),0,(size,size)).convert('RGB').resize((256,256),Image.Resampling.LANCZOS);key=(row.qc_case_id,int(row.panel_slot));images[key]=patch
     metrics.append({'qc_case_id':row.qc_case_id,'qc_partition':row.qc_partition,'panel_slot':int(row.panel_slot),'patient_id':row.patient_id,'file_name':fname,'x_level0':x,'y_level0':y,**measure(np.asarray(patch))})
  except Exception as e:errors.append({'file_name':fname,'error':str(e)})
 for case,g in selected.groupby('qc_case_id',sort=True):
  canvas=Image.new('RGB',(1024,560),'white');draw=ImageDraw.Draw(canvas);complete=True
  for _,row in g.iterrows():
   j=int(row.panel_slot);im=images.get((case,j))
   if im is None:complete=False;continue
   left=j%4*256;top=j//4*280;canvas.paste(im,(left,top));draw.text((left+3,top+258),f'{case} field {j+1}',fill='black')
  folder=out/(g.qc_partition.iloc[0]+'_panels');folder.mkdir(exist_ok=True)
  if complete:canvas.save(folder/(case+'.png'))
 pd.DataFrame(metrics).to_csv(out/'qc_reference_measurements.csv',index=False)
 (out/'summary.json').write_text(json.dumps({'status':'COMPLETE' if not errors else 'INCOMPLETE','patients_requested':40,'patches_requested':320,'patches_read':len(metrics),'errors':errors,'dataset_or_models_changed':False,'diagnostic_review_performed':False,'software':{'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'openslide':openslide.__version__}},indent=2))
 archive=out.with_suffix('.zip')
 with zipfile.ZipFile(archive,'w',compression=zipfile.ZIP_DEFLATED) as z:
  for p in out.rglob('*'):
   if p.is_file():z.write(p,p.relative_to(out))
 print('UPLOAD REPORT:',archive,flush=True)
 if errors:raise SystemExit(1)
if __name__=='__main__':main()
