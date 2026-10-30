"""Build an outcome-independent CPTAC image/specimen linkage and transfer inventory."""
from __future__ import annotations
import argparse
from collections import Counter, defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import time
import urllib.parse
import urllib.request

USER_AGENT = 'BRCA-CPTAC-Research-Audit/1.0'
TIFF_MAGICS = (b'II*\x00', b'MM\x00*', b'II+\x00', b'MM\x00+')


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def write_csv(path, rows):
    if not rows:
        Path(path).unlink(missing_ok=True)
        return
    with Path(path).open('w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def canonical_patient(identifier):
    if re.fullmatch(r'X\d{2}BR\d{3}', identifier): return identifier[1:]
    if re.fullmatch(r'\d{2}BR\d{3}|CPT\d{6}', identifier): return identifier
    raise ValueError('Unrecognized external identifier: '+identifier)


def specimen_matches(slide, specimen):
    # Only a documented suffix is removed. Section suffixes remain meaningful.
    stem = re.sub(r'_D1_D1(?=_\d+$|$)', '', slide)
    for key in ['stripped_label', 'specimen_id']:
        value = str(specimen.get(key, ''))
        value = re.sub(r'_D1_D1(?=_\d+$|$)', '', value)
        if stem == value: return 'exact'
    alias = re.fullmatch(r'\[([0-9a-f]{2})\]-([0-9a-f]{6})_D1_D1(_\d+)?', str(specimen.get('specimen_id','')), re.I)
    if alias:
        match = re.fullmatch(r'([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{6})(_\d+)?', stem, re.I)
        if match and match[1].startswith(alias[1]) and match[1].endswith('-'+alias[2]) and (match[2] or '') == (alias[3] or ''):
            return 'unique_abbreviated_specimen_alias'
    return None


def fetch_cached(url, path):
    path = Path(path)
    if path.exists(): return json.loads(path.read_text())
    for attempt in range(3):
        try:
            req = urllib.request.Request(url, headers={'User-Agent':USER_AGENT, 'Accept':'application/json'})
            with urllib.request.urlopen(req, timeout=35) as response: raw = response.read()
            result = json.loads(raw)
            path.write_bytes(raw)
            return result
        except Exception:
            if attempt == 2: raise
            time.sleep(1+attempt)


def node_value(node, key, default=''):
    values = node.get(key, [])
    return values[0].get('value', default) if values else default


def probe_file(url, path):
    path = Path(path)
    if path.exists(): return json.loads(path.read_text())
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'Range':'bytes=0-4095','User-Agent':USER_AGENT})
            with urllib.request.urlopen(request,timeout=35) as response:
                header = response.read(4096)
                content_range = response.headers.get('Content-Range','')
                match = re.fullmatch(r'bytes 0-4095/(\d+)',content_range)
                if response.status != 206 or not match: raise ValueError('Server did not honour bounded range: '+url)
                if len(header)!=4096 or not header.startswith(TIFF_MAGICS): raise ValueError('Missing TIFF header: '+url)
                result = {'url':url, 'status':response.status,'bytes':int(match[1]),
                    'last_modified':response.headers.get('Last-Modified',''),
                    'etag':response.headers.get('ETag',''), 'prefix_bytes':4096,
                    'prefix_sha256':hashlib.sha256(header).hexdigest(), 'tiff_magic_hex':header[:4].hex()}
            path.write_text(json.dumps(result,indent=2)); return result
        except Exception:
            if attempt == 2: raise
            time.sleep(1+attempt)


def prepare(sources, clinical, out, workers=4):
    out.mkdir(parents=True,exist_ok=True)
    data = out/'data'; data.mkdir(exist_ok=True)
    evidence = out/'evidence'; evidence.mkdir(exist_ok=True)
    nodes = evidence/'pathdb_nodes'; nodes.mkdir(exist_ok=True)
    probes = evidence/'file_probes'; probes.mkdir(exist_ok=True)
    inventory = list(csv.DictReader((sources/'pathdb_cptac_brca_inventory.csv').open(encoding='utf-8')))
    cases = {str(c['case_id']):c for c in json.loads((sources/'archived_BR.json').read_text())}
    samples = json.loads((clinical/'samples.json').read_text())
    canonical = {canonical_patient(r['patientId']):r['patientId'] for r in samples}
    assert len(canonical)==len(samples)==122
    by_patient = defaultdict(list)
    for r in inventory:
        assert r['collection']=='CPTAC-BRCA' and '\ufffd' not in ''.join(r.values())
        by_patient[r['patient_id']].append(r)
    candidates,excluded = [],[]
    for row in inventory:
        pid = row['patient_id']
        if pid not in canonical:
            excluded.append({**row,'reason':'NO_MATCHED_MOLECULAR_RECORD'});continue
        matches = [(sp,kind) for sp in cases.get(pid,{}).get('specimens',[]) if (kind:=specimen_matches(row['slide_id'],sp))]
        if len(matches)!=1:
            excluded.append({**row,'reason':'NO_UNIQUE_SPECIMEN_LINK'});continue
        sp, kind = matches[0]
        # A shortened alias must map back uniquely within this patient's entire image inventory.
        if kind!='exact' and sum(bool(specimen_matches(r['slide_id'],sp)) for r in by_patient[pid])!=1:
            excluded.append({**row,'reason':'AMBIGUOUS_ABBREVIATED_ALIAS'});continue
        if sp.get('specimen_type')!='tumor_tissue':
            excluded.append({**row,'reason':'NOT_DOCUMENTED_TUMOR_TISSUE'});continue
        candidates.append({**row,'cbio_patient_id':canonical[pid], 'specimen_id':sp['specimen_id'],
                           'specimen_type':sp['specimen_type'],'specimen_link_method':kind})
    write_csv(data/'excluded_images.csv',excluded)
    write_csv(data/'candidate_linkage_before_remote_checks.csv',candidates)
    demographics = json.loads((clinical/'patient_clinical_data.json').read_text())
    attr = {(r['patientId'],r['clinicalAttributeId']):r['value'] for r in demographics}
    crosswalk = []
    for pid, cbio in sorted(canonical.items()):
        case = cases.get(pid,{})
        row = {'cbio_patient_id':cbio,'tcia_patient_id':pid,
               'id_rule':'remove_one_leading_X_from_XNNBRNNN' if cbio.startswith('X') else 'identity',
               'pathdb_slides':len(by_patient.get(pid,[])),
               'linked_tumor_slides':sum(r['patient_id']==pid for r in candidates),
               'ESAC_case_present':bool(case)}
        for source, target in [('age','AGE'),('sex','SEX')]:
            a,b=str(case.get(source,'')).strip(),attr.get((cbio,target),'').strip()
            row[source+'_esac']=a;row[source+'_cbio']=b
            agreement='MISSING_ONE_OR_BOTH' if not a or not b or a.lower() in ['not available','na'] or b.lower()=='na' else 'AGREE' if a.lower()==b.lower() else 'DISAGREE'
            if source=='age' and agreement=='DISAGREE':
                try:
                    if math.floor(float(a))==math.floor(float(b)):agreement='AGREE_COMPLETED_YEARS'
                except ValueError:pass
            row[source+'_agreement']=agreement
        crosswalk.append(row)
    write_csv(data/'patient_crosswalk.csv',crosswalk)

    def inspect(row):
        nid = row['camic_id']; url = f'https://pathdb.cancerimagingarchive.net/node/{nid}?_format=json'
        node = fetch_cached(url,nodes/(nid+'.json'))
        original_nid=nid;resolution='EXPORTED_NODE_CONFIRMED'
        if node_value(node,'clinicaltrialsubjectid')!=row['patient_id'] or node_value(node,'imageid')!=row['slide_id']:
            # The public cohort export has two swapped node identifiers. Resolve
            # only by a unique, already observed node with this exact patient/image.
            alternatives=[]
            for other in candidates:
                other_path=nodes/(other['camic_id']+'.json')
                if other['patient_id']!=row['patient_id'] or not other_path.exists():continue
                observed=json.loads(other_path.read_text())
                if node_value(observed,'clinicaltrialsubjectid')==row['patient_id'] and node_value(observed,'imageid')==row['slide_id']:
                    alternatives.append((other['camic_id'],observed))
            if len(alternatives)!=1:raise ValueError('Unresolved inventory/node identity conflict: '+nid)
            nid,node=alternatives[0];resolution='UNIQUE_EXACT_PATIENT_AND_IMAGE_NODE_CORRECTION'
        if node_value(node,'studyid')!='CPTAC_BRCA':raise ValueError('Unexpected slide study')
        file_url = node['field_wsiimage'][0]['url']
        parsed = urllib.parse.urlsplit(file_url)
        if parsed.hostname!='pathdb.cancerimagingarchive.net':raise ValueError('Unexpected file host')
        file_url = urllib.parse.urlunsplit(('https',parsed.netloc,parsed.path,parsed.query,''))
        expected = row['patient_id']+'-'+row['slide_id']+'.svs'
        if urllib.parse.unquote(Path(parsed.path).name)!=expected:raise ValueError('File/patient identity conflict: '+nid)
        probe = probe_file(file_url,probes/(nid+'.json'))
        return {'patient_id':row['patient_id'],'cbio_patient_id':row['cbio_patient_id'],
                'slide_id':row['slide_id'],'pathdb_node_id':nid,'file_name':expected,'url':file_url,
                'exported_node_id':original_nid,'node_resolution':resolution,
                'bytes':probe['bytes'],'last_modified':probe['last_modified'],'etag':probe['etag'],
                'prefix_bytes':4096,'prefix_sha256':probe['prefix_sha256'],
                'width':node_value(node,'imagedvolumewidth'),'height':node_value(node,'imagedvolumeheight'),
                'mpp_x':node_value(node,'referencepixelphysicalvaluex'),'mpp_y':node_value(node,'referencepixelphysicalvaluey'),
                'specimen_id':row['specimen_id'],'specimen_type':row['specimen_type'],
                'specimen_link_method':row['specimen_link_method'],
                'visual_qc_status':'PENDING','same_assayed_block_verified':'NO_PATIENT_LEVEL_LINK_ONLY'}
    ready, failures = [],[]
    with ThreadPoolExecutor(max_workers=workers) as executor:
        futures = {executor.submit(inspect,row):row for row in candidates}
        for i,future in enumerate(as_completed(futures),1):
            row = futures[future]
            try:ready.append(future.result())
            except Exception as e:
                failures.append({'patient_id':row['patient_id'],'slide_id':row['slide_id'],'error':str(e)})
                print(f"UNRESOLVED {row['camic_id']}: {e}",flush=True)
            if i%20==0 or i==len(candidates):print(f'REMOTE INVENTORY {i}/{len(candidates)} checked; {len(failures)} unresolved',flush=True)
    ready.sort(key=lambda r:(r['patient_id'],r['slide_id']))
    write_csv(data/'download_manifest.csv',ready)
    write_csv(data/'remote_check_failures.csv',failures)
    # Technical pilot only; determined from identifiers, with no outcome-based selection.
    pilot=[]
    for pid in sorted({r['patient_id'] for r in ready})[:3]:pilot.append(next(r for r in ready if r['patient_id']==pid))
    write_csv(data/'pilot_manifest.csv',pilot)
    summary = {'status':'LINKAGE_AND_REMOTE_PREFLIGHT_COMPLETE' if not failures else 'REMOTE_PREFLIGHT_INCOMPLETE',
        'molecular_patients':len(samples),'pathdb_inventory_rows':len(inventory),
        'pathdb_inventory_patient_ids':len(by_patient),'patients_with_any_pathdb_image':sum(bool(by_patient.get(pid)) for pid in canonical),
        'candidate_linked_tumor_slides':len(candidates),'ready_slides':len(ready),
        'ready_patients':len({r['patient_id'] for r in ready}),'ready_bytes':sum(r['bytes'] for r in ready),
        'linkage_counts':dict(Counter(r['specimen_link_method'] for r in ready)),
        'excluded_counts':dict(Counter(r['reason'] for r in excluded)), 'remote_failures':len(failures),
        'patients_without_any_image':[pid for pid in canonical if not by_patient.get(pid)],
        'pilot_slides':len(pilot),'pilot_bytes':sum(r['bytes'] for r in pilot),
        'demographic_agreement':{x:dict(Counter(r[x+'_agreement'] for r in crosswalk)) for x in ['age','sex']},
        'same_assayed_block_proven':False,'visual_QC_complete':False,'external_features_extracted':False,
        'external_metrics_computed':False,'pretraining_overlap':'NOT_YET_ESTABLISHED',
        'selection':'All uniquely linked tumor-tissue slides with matching patient records and passing remote checks; no selection by outcome or model score'}
    (out/'linkage_summary.json').write_text(json.dumps(summary,indent=2))
    print(json.dumps(summary,indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--sources',type=Path,required=True);p.add_argument('--clinical',type=Path,required=True)
    p.add_argument('--out',type=Path,required=True);p.add_argument('--workers',type=int,default=4)
    a=p.parse_args();prepare(a.sources,a.clinical,a.out,a.workers)
