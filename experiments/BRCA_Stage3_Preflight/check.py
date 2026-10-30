"""Offline environment and exact boundary-patch audit; no encoder inference."""
import argparse, csv, hashlib, importlib.metadata, io, json, os, sys, zipfile
from pathlib import Path
from datetime import datetime

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(4*1024*1024),b''): h.update(b)
    return h.hexdigest()

def main():
    p=argparse.ArgumentParser(); p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--report',type=Path); a=p.parse_args()
    report=a.report or a.project_root/'BRCA_ExternalValidation_Stage2_FullQC'/'brca_stage2_full_qc_20261021_041052_426646.zip'
    if not report.is_file(): raise FileNotFoundError('QC report missing; pass -Report with its full ZIP path.')
    out=Path(__file__).resolve().parent/('report_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f'));out.mkdir()
    result={'python':sys.version,'packages':{},'errors':[], 'encoder_inference_run':False,'network_access_used':False}
    for name in ['torch','torchvision','timm','transformers','huggingface-hub','openslide-python','Pillow','numpy','pandas']:
        try: result['packages'][name]=importlib.metadata.version(name)
        except importlib.metadata.PackageNotFoundError: result['packages'][name]=None
    try:
        import torch
        result['cuda_available']=torch.cuda.is_available();result['torch_cuda']=torch.version.cuda
        result['gpus']=[{'name':torch.cuda.get_device_name(i),'memory_bytes':torch.cuda.get_device_properties(i).total_memory} for i in range(torch.cuda.device_count())]
    except Exception as e: result['errors'].append('torch: '+str(e))
    roots={Path.home()/'.cache'/'huggingface'/'hub'}
    if os.environ.get('HF_HOME'): roots.add(Path(os.environ['HF_HOME'])/'hub')
    for key in ['HF_HUB_CACHE','HUGGINGFACE_HUB_CACHE']:
        if os.environ.get(key): roots.add(Path(os.environ[key]))
    result['cached_models']=[]
    for root in sorted(roots):
        for repo in ['models--MahmoodLab--uni','models--histai--hibou-b']:
            folder=root/repo/'snapshots'
            if not folder.exists(): continue
            for snapshot in sorted(folder.iterdir()):
                if not snapshot.is_dir(): continue
                item={'repository':repo,'revision':snapshot.name,'path':str(snapshot),'files':[]}
                for f in sorted(snapshot.rglob('*')):
                    if f.is_file() and f.suffix in ['.json','.safetensors','.bin','.py']:
                        print('Hashing cached file:',f.name,flush=True)
                        item['files'].append({'name':str(f.relative_to(snapshot)),'bytes':f.stat().st_size,'sha256':digest(f)})
                result['cached_models'].append(item)
    try:
        import openslide
        from PIL import Image, ImageDraw
        with zipfile.ZipFile(report) as z:
            summary=json.loads(z.read('summary.json'));result['qc_coordinates_sha256']=summary['coordinates_sha256']
            result['boundary_patches']=[]
            for slide in summary['slides']:
                if not slide['retained_edge_crossing_patches']:continue
                path=Path(slide['source_path']);checkpoint='slide_checkpoints/'+slide['checkpoint_csv'].replace('\\','/').split('/')[-1]
                data=z.read(checkpoint)
                if hashlib.sha256(data).hexdigest()!=slide['checkpoint_sha256']:raise ValueError('Changed checkpoint '+checkpoint)
                rows=list(csv.DictReader(io.StringIO(data.decode('utf-8-sig'))))
                with openslide.OpenSlide(str(path)) as wsi:
                    width,height=wsi.dimensions
                    for row in rows:
                        if row['retain_after_new_locked_qc'].lower() not in ['true','1','1.0']:continue
                        x,y,size=[int(row[k]) for k in ['x_level0','y_level0','read_size_level0']]
                        if x>=0 and y>=0 and x+size<=width and y+size<=height:continue
                        im=wsi.read_region((x,y),0,(size,size)).convert('RGB').resize((256,256),Image.Resampling.LANCZOS)
                        name=str(slide['site_slide_index'])+'_'+str(row['patch_id'])+'.png';im.save(out/name)
                        result['boundary_patches'].append({'file_name':slide['file_name'],'patch_id':row['patch_id'],'x':x,'y':y,'read_size':size,'slide_dimensions':[width,height],'preview':name})
        if len(result['boundary_patches'])!=5:result['errors'].append('Expected 5 boundary patches; recovered '+str(len(result['boundary_patches'])))
    except Exception as e:result['errors'].append('boundary audit: '+type(e).__name__+': '+str(e))
    result['status']='PREFLIGHT_REVIEW_REQUIRED'
    (out/'report.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    destination=out.with_suffix('.zip')
    with zipfile.ZipFile(destination,'w',zipfile.ZIP_DEFLATED) as z:
        for f in sorted(out.iterdir()):z.write(f,f.name)
    print('UPLOAD REPORT ZIP:',destination,flush=True)
    if result['errors']:print('Recorded issues:',json.dumps(result['errors']))

if __name__=='__main__':main()
