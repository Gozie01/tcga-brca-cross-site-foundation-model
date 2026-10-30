"""Collect bounded project artifacts for provenance review, without loading models."""
import argparse,csv,hashlib,importlib.metadata,json,os,sys,zipfile
from pathlib import Path
from datetime import datetime

MAX_FILE=20*1024**2
MAX_TOTAL=200*1024**2
SKIP={'.git','node_modules','__pycache__','site-packages','features_v1','runs'}

def wanted(path):
    name=path.name.lower();ext=path.suffix.lower();parts=[p.lower() for p in path.parts]
    if any(x.startswith('.venv') or x in SKIP for x in parts):return False
    if any(x in name for x in ['token','credential','secret','password','auth.json']):return False
    if ext in ['.py','.ps1','.joblib','.pkl','.pickle']:return True
    if ext=='.npy':return 'patient' in name and 'patch' not in name
    if ext in ['.json','.yaml','.yml']:
        return any(x in name for x in ['metadata','freeze','lock','summary','calibrat','cohort','label','clinical','patient','specimen','link','mapping','audit','model','sample','molecular','mutation','pam50','subtype'])
    if ext in ['.csv','.tsv','.txt','.md']:
        return any(x in name for x in ['patient','cohort','split','label','clinical','specimen','link','mapping','endpoint','prediction','metric','readme','requirements','manifest','sample','molecular','mutation','pam50','subtype']) and not any(x in name for x in ['patch','coordinate','embedding_index'])
    return False

def roots(base):
    prefixes=('tcga_brca','brca_revision','brca_externalvalidation_stage1')
    return sorted(p for p in base.iterdir() if p.is_dir() and p.name.lower().startswith(prefixes))

def walk(root):
    for directory,folders,files in os.walk(root,followlinks=False):
        folders[:]=sorted(f for f in folders if f.lower() not in SKIP and not f.lower().startswith('.venv') and not Path(directory,f).is_symlink())
        for name in sorted(files):
            p=Path(directory,name)
            if not p.is_symlink():yield p

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--project-root',type=Path,required=True);a=ap.parse_args()
    base=a.project_root.resolve()
    if not base.is_dir():raise FileNotFoundError(base)
    out=Path(__file__).resolve().parent/('brca_asset_audit_'+datetime.now().strftime('%Y%m%d_%H%M%S_%f')+'.zip')
    report={'status':'ASSET_COLLECTION_FOR_REVIEW','roots':[],'files':[],'skipped':[],'errors':[],'python':sys.version,'packages':{},'models_loaded':False,'training_run':False,'predictions_computed':False,'collection_byte_limit':MAX_TOTAL}
    for package in ['scikit-learn','joblib','numpy','pandas','torch','timm','transformers']:
        try:report['packages'][package]=importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:report['packages'][package]=None
    total=0
    with zipfile.ZipFile(out,'w',zipfile.ZIP_DEFLATED) as z:
        for root in roots(base):
            report['roots'].append(str(root));print('Reviewing project folder:',root,flush=True)
            for path in walk(root):
                if not wanted(path):continue
                relative=path.relative_to(base).as_posix()
                try:
                    size=path.stat().st_size
                    if size>MAX_FILE or total+size>MAX_TOTAL:
                        report['skipped'].append({'path':relative,'bytes':size,'reason':'size_limit'});continue
                    data=path.read_bytes()
                    if len(data)>MAX_FILE or total+len(data)>MAX_TOTAL:raise ValueError('File size changed beyond collection limit')
                    z.writestr('assets/'+relative,data);total+=len(data)
                    report['files'].append({'path':relative,'bytes':len(data),'sha256':hashlib.sha256(data).hexdigest()})
                except Exception as e:report['errors'].append({'path':relative,'error':type(e).__name__+': '+str(e)})
        report['collected_bytes']=total;report['file_count']=len(report['files'])
        report['model_file_count']=sum(Path(x['path']).suffix.lower() in ['.joblib','.pkl','.pickle'] for x in report['files'])
        if not report['roots']:report['errors'].append({'error':'No matching BRCA project folders found under ProjectRoot'})
        if report['errors'] or report['skipped']:report['status']='ASSET_COLLECTION_WITH_GAPS'
        z.writestr('asset_inventory.json',json.dumps(report,indent=2))
    print('Collected files:',report['file_count'],'Model files:',report['model_file_count'],'Skipped:',len(report['skipped']),flush=True)
    print('UPLOAD REPORT ZIP:',out,flush=True)

if __name__=='__main__':main()
