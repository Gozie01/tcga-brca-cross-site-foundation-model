"""Checkpoint primitives: all identities and arrays validated before reuse."""
import csv, hashlib, io, json, os
from pathlib import Path
import numpy as np

def sha(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for b in iter(lambda:f.read(4194304),b''):h.update(b)
    return h.hexdigest()

def identity(obj):
    return hashlib.sha256(json.dumps(obj,sort_keys=True,default=str).encode()).hexdigest()

def atomic_json(path,obj):
    path=Path(path);tmp=path.with_suffix('.json.tmp')
    tmp.write_text(json.dumps(obj,indent=2,default=str),encoding='utf-8');os.replace(tmp,path)

def atomic_npy(path,array):
    path=Path(path);tmp=path.with_suffix('.npy.tmp')
    with tmp.open('wb') as f:np.save(f,array,allow_pickle=False)
    os.replace(tmp,path)

def read_rows(z,slide):
    key='slide_checkpoints/'+slide['checkpoint_csv'].replace('\\','/').split('/')[-1]
    data=z.read(key)
    if hashlib.sha256(data).hexdigest()!=slide['checkpoint_sha256']:raise ValueError('Coordinate checksum changed: '+key)
    rows=[r for r in csv.DictReader(io.StringIO(data.decode('utf-8-sig'))) if r['retain_after_new_locked_qc'].lower() in ['true','1','1.0']]
    if len(rows)!=slide['final_retained'] or not rows:raise ValueError('Retained count mismatch')
    coords=[(int(r['x_level0']),int(r['y_level0'])) for r in rows]
    if len(set(coords))!=len(rows):raise ValueError('Duplicate retained coordinates')
    if any(r['file_name']!=slide['file_name'] or r['patient_id']!=slide['patient_id'] for r in rows):raise ValueError('Slide identity mismatch')
    return rows

def valid_array(a,rows,dim):
    return a.shape==(rows,dim) and a.dtype==np.float32 and bool(np.isfinite(a).all()) and bool(np.all(np.linalg.norm(a,axis=1)>0))

def checkpoint(folder,token,expected,dim):
    array=folder/(token+'.npy');meta=folder/(token+'.json')
    if not meta.exists():return None
    m=json.loads(meta.read_text(encoding='utf-8'))
    if m.get('status')!='COMPLETE' or m.get('identity')!=expected:raise ValueError('Checkpoint identity mismatch: '+token)
    if not array.is_file() or sha(array)!=m['feature_sha256']:raise ValueError('Checkpoint feature hash mismatch: '+token)
    a=np.load(array,allow_pickle=False,mmap_mode='r')
    if not valid_array(a,m['rows'],dim):raise ValueError('Checkpoint values invalid: '+token)
    return m

def commit(folder,token,a,expected,details):
    atomic_npy(folder/(token+'.npy'),a)
    m={**details,'status':'COMPLETE','identity':expected,'rows':len(a),'feature_sha256':sha(folder/(token+'.npy'))}
    atomic_json(folder/(token+'.json'),m)
    return m

def pool(folder,slides,dim):
    means=[];index=[]
    for s in slides:
        a=np.load(folder/(s['token']+'.npy'),allow_pickle=False,mmap_mode='r')
        means.append(np.asarray(a,dtype=np.float64).mean(axis=0).astype(np.float32));index.append({'patient_id':s['patient_id'],'file_name':s['file_name'],'patches':s['rows']})
    slide_features=np.stack(means)
    patients=sorted(set(r['patient_id'] for r in index))
    features=np.stack([slide_features[[i for i,r in enumerate(index) if r['patient_id']==p]].astype(np.float64).mean(axis=0).astype(np.float32) for p in patients])
    assert features.shape==(len(patients),dim)
    return slide_features,index,features,patients
