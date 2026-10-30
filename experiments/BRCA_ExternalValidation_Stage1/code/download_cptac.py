"""Download the locked CPTAC SVS inventory with restart and integrity checks."""
from __future__ import annotations
import argparse
import csv
from datetime import datetime, timezone
import hashlib
import http.client
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
import urllib.request
import zipfile

ROOT=Path(__file__).resolve().parents[1]
MAGICS=(b'II*\x00',b'MM\x00*',b'II+\x00',b'MM\x00+')
AGENT='BRCA-CPTAC-Research-Download/1.0'


def atomic_json(path, value):
    path=Path(path); tmp=path.with_name(path.name+'.tmp')
    tmp.write_text(json.dumps(value,indent=2),encoding='utf-8');tmp.replace(path)


def digest(path):
    h=hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda:f.read(8*1024*1024),b''):h.update(block)
    return h.hexdigest()


def verify_package(root=ROOT):
    lock=json.loads((root/'PACKAGE_SHA256.json').read_text())
    for name, expected in lock.items():
        if digest(root/name)!=expected:raise ValueError('Package file changed: '+name)
    return len(lock)


def read_manifest(path):
    with Path(path).open(encoding='utf-8',newline='') as f:rows=list(csv.DictReader(f))
    seen=set()
    for r in rows:
        url=urllib.parse.urlsplit(r['url'])
        if url.scheme!='https' or url.hostname!='pathdb.cancerimagingarchive.net':raise ValueError('Unexpected download URL')
        if not re.fullmatch(r'\d{2}BR\d{3}',r['patient_id']):raise ValueError('Unexpected patient identifier')
        if Path(r['file_name']).name!=r['file_name'] or '/' in r['file_name'] or '\\' in r['file_name']:raise ValueError('Unsafe file name')
        if r['file_name']!=r['patient_id']+'-'+r['slide_id']+'.svs':raise ValueError('Patient/slide filename mismatch')
        if r['file_name'] in seen:raise ValueError('Duplicate file name')
        if int(r['bytes'])<4096 or int(r['prefix_bytes'])!=4096:raise ValueError('Invalid remote size')
        seen.add(r['file_name'])
    if not rows:raise ValueError('Empty manifest')
    return rows


def remote_prefix(row, opener=urllib.request.urlopen):
    req=urllib.request.Request(row['url'],headers={'Range':'bytes=0-4095','User-Agent':AGENT})
    with opener(req,timeout=45) as r:
        prefix=r.read(4096)
        if r.status!=206 or r.headers.get('Content-Range')!=f"bytes 0-4095/{row['bytes']}":raise ValueError('Remote size/range changed')
        if not prefix.startswith(MAGICS) or hashlib.sha256(prefix).hexdigest()!=row['prefix_sha256']:raise ValueError('Remote TIFF prefix changed')
        if row.get('last_modified') and r.headers.get('Last-Modified')!=row['last_modified']:raise ValueError('Remote modification time changed')
        if row.get('etag') and r.headers.get('ETag')!=row['etag']:raise ValueError('Remote ETag changed')


def local_prefix(path, row):
    with path.open('rb') as f:prefix=f.read(4096)
    if len(prefix)!=4096 or hashlib.sha256(prefix).hexdigest()!=row['prefix_sha256']:raise ValueError('Local TIFF prefix mismatch: '+str(path))


def download_one(row, destination, opener=urllib.request.urlopen, progress=print):
    expected=int(row['bytes']); destination=Path(destination)
    destination.mkdir(parents=True,exist_ok=True)
    final=destination/row['file_name']; part=final.with_name(final.name+'.part')
    state=final.with_name(final.name+'.part.json'); receipt=final.with_name(final.name+'.receipt.json')
    identity={k:row.get(k,'') for k in ['url','bytes','prefix_sha256','last_modified','etag']}
    remote_prefix(row,opener)
    if final.exists():
        if not receipt.exists():raise ValueError('Existing file lacks this downloader\'s receipt; no overwrite: '+str(final))
        saved=json.loads(receipt.read_text())
        if saved['remote_identity']!=identity or final.stat().st_size!=expected or digest(final)!=saved['sha256']:raise ValueError('Existing file/receipt mismatch')
        local_prefix(final,row);progress('VERIFIED EXISTING '+row['file_name']);return saved
    if part.exists():
        if not state.exists() or json.loads(state.read_text())!=identity:raise ValueError('Partial file identity does not match locked source')
        if part.stat().st_size>expected:raise ValueError('Partial file exceeds expected size')
        if part.stat().st_size>=4096:local_prefix(part,row)
    else:
        atomic_json(state,identity)
    for attempt in range(4):
        offset=part.stat().st_size if part.exists() else 0
        if offset==expected:break
        request=urllib.request.Request(row['url'],headers={'Range':f'bytes={offset}-','User-Agent':AGENT})
        try:
            with opener(request,timeout=90) as response:
                expected_range=f'bytes {offset}-{expected-1}/{expected}'
                if response.status!=206 or response.headers.get('Content-Range')!=expected_range:raise ValueError('Server returned an incorrect resume range')
                if row.get('last_modified') and response.headers.get('Last-Modified')!=row['last_modified']:raise ValueError('Remote object changed during transfer')
                if row.get('etag') and response.headers.get('ETag')!=row['etag']:raise ValueError('Remote ETag changed during transfer')
                count=offset;last_update=time.monotonic()
                with part.open('ab') as f:
                    while True:
                        chunk=response.read(min(2*1024*1024,expected-count+1))
                        if not chunk:break
                        if count+len(chunk)>expected:raise ValueError('Server sent excess bytes')
                        f.write(chunk);count+=len(chunk)
                        if time.monotonic()-last_update>=5:
                            progress(f"{row['file_name']}: {count/1e6:.1f}/{expected/1e6:.1f} MB");last_update=time.monotonic()
                if count!=expected:raise OSError('Transfer ended before the expected byte count')
            break
        except ValueError:raise
        except (OSError, TimeoutError, http.client.IncompleteRead):
            if attempt==3:raise
            progress('Transfer interrupted; retrying from saved bytes.');time.sleep(min(2**attempt,8))
    if part.stat().st_size!=expected:raise ValueError('Incomplete file')
    local_prefix(part,row)
    content_hash=digest(part)
    record={'file_name':row['file_name'],'patient_id':row['patient_id'],'bytes':expected,'sha256':content_hash,
            'remote_identity':identity,'completed_utc':datetime.now(timezone.utc).isoformat(),
            'source_full_file_checksum_available':False,
            'integrity_basis':'Locked remote size/modification/prefix; correct HTTP ranges; local full-file SHA-256 receipt'}
    # Write the receipt first so interruption after the atomic rename is recoverable.
    atomic_json(receipt,record);part.replace(final)
    state.unlink(missing_ok=True)
    progress('COMPLETE '+row['file_name']);return record


def run(args):
    print(f'PACKAGE VERIFIED: {verify_package()} files',flush=True)
    manifest=ROOT/'data'/('download_manifest.csv' if args.mode=='DownloadAll' else 'pilot_manifest.csv')
    rows=read_manifest(manifest)
    args.wsi_root.mkdir(parents=True,exist_ok=True)
    existing=sum(min((args.wsi_root/r['patient_id']/(r['file_name']+'.part')).stat().st_size,int(r['bytes']))
                 if (args.wsi_root/r['patient_id']/(r['file_name']+'.part')).exists() else int(r['bytes'])
                 if (args.wsi_root/r['patient_id']/r['file_name']).exists() else 0 for r in rows)
    remaining=sum(int(r['bytes']) for r in rows)-existing
    free=shutil.disk_usage(args.wsi_root).free
    print(f'{len(rows)} slides; {sum(int(r["bytes"]) for r in rows)/1e9:.2f} GB total; {remaining/1e9:.2f} GB remaining; {free/1e9:.1f} GB free',flush=True)
    if args.mode=='Audit':return
    if free < remaining+2*1024**3:raise OSError('Insufficient disk space; allow 2 GiB above remaining slide bytes')
    run_dir=args.project_root/'brca_external_stage1'/datetime.now(timezone.utc).strftime('%Y%m%d_%H%M%S_%f')
    run_dir.mkdir(parents=True,exist_ok=False)
    shutil.copyfile(manifest,run_dir/'selected_manifest.csv')
    shutil.copyfile(ROOT/'linkage_summary.json',run_dir/'linkage_summary.json')
    try:
        results=[]
        for i,row in enumerate(rows,1):
            print(f'DOWNLOAD {i}/{len(rows)} {row["patient_id"]}',flush=True)
            results.append(download_one(row,args.wsi_root/row['patient_id'],progress=lambda x:print(x,flush=True)))
            atomic_json(run_dir/'download_receipts.json',results)
        if args.slide_python:
            subprocess.run([str(args.slide_python),str(ROOT/'code/check_slides.py'),'--manifest',str(manifest),
                '--wsi-root',str(args.wsi_root),'--out',str(run_dir/'slide_preflight.json')],check=True)
        else:
            atomic_json(run_dir/'slide_preflight.json',{'status':'NOT_RUN_NO_SLIDE_PYTHON'})
        atomic_json(run_dir/'run_status.json',{'status':'DOWNLOAD_AND_PREFLIGHT_COMPLETE' if args.slide_python else 'DOWNLOAD_COMPLETE_PREFLIGHT_PENDING',
            'mode':args.mode,'slides':len(rows),'wsi_root':str(args.wsi_root),'image_feature_extraction_run':False,
            'external_predictions_computed':False})
    except BaseException as e:
        atomic_json(run_dir/'run_status.json',{'status':'INTERRUPTED' if isinstance(e,KeyboardInterrupt) else 'FAILED',
            'mode':args.mode,'error':str(e),'error_type':type(e).__name__,
            'image_feature_extraction_run':False,'external_predictions_computed':False})
        raise
    finally:
        archive=args.project_root/('brca_external_stage1_'+run_dir.name+'.zip')
        with zipfile.ZipFile(archive,'w',zipfile.ZIP_DEFLATED) as z:
            for p in sorted(run_dir.rglob('*')):
                if p.is_file():z.write(p,p.relative_to(run_dir))
        print('Upload this report ZIP: '+str(archive),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--mode',choices=['Audit','Pilot','DownloadAll'],default='Audit')
    p.add_argument('--wsi-root',type=Path,required=True);p.add_argument('--project-root',type=Path,required=True)
    p.add_argument('--slide-python',type=Path)
    try:run(p.parse_args())
    except KeyboardInterrupt:print('\nCancelled; partial downloads are retained for the next run.',flush=True);sys.exit(130)
