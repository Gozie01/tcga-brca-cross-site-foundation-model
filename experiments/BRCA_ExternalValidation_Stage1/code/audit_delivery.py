"""Independently reconcile delivered transfer rows with current nodes and reference data."""
import csv
from collections import Counter
import hashlib
import json
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]


if __name__=='__main__':
    rows=list(csv.DictReader((ROOT/'data/download_manifest.csv').open()))
    assert len(rows)==391 and len({r['patient_id'] for r in rows})==120
    assert len({r['url'] for r in rows})==len({r['file_name'] for r in rows})==391
    corrected=[]
    for r in rows:
        node=json.loads((ROOT/'evidence/pathdb_nodes'/(r['pathdb_node_id']+'.json')).read_text())
        probe=json.loads((ROOT/'evidence/file_probes'/(r['pathdb_node_id']+'.json')).read_text())
        assert node['clinicaltrialsubjectid'][0]['value']==r['patient_id']
        assert node['imageid'][0]['value']==r['slide_id']
        assert node['studyid'][0]['value']=='CPTAC_BRCA'
        assert probe['url']==r['url'] and str(probe['bytes'])==r['bytes'] and probe['prefix_sha256']==r['prefix_sha256']
        assert probe['status']==206 and probe['tiff_magic_hex'] in ['49492a00','4d4d002a','49492b00','4d4d002b']
        assert node['field_wsiimage'][0]['url'].removeprefix('http://')==r['url'].removeprefix('https://')
        if r['exported_node_id']!=r['pathdb_node_id']:corrected.append((r['exported_node_id'],r['pathdb_node_id']))
    assert sorted(corrected)==[('211933','211934'),('211934','211933')]
    summary=json.loads((ROOT/'linkage_summary.json').read_text())
    assert summary['ready_bytes']==sum(int(r['bytes']) for r in rows) and summary['remote_failures']==0
    assert Counter(r['specimen_link_method'] for r in rows)=={'exact':298,'unique_abbreviated_specimen_alias':93}
    pilot=list(csv.DictReader((ROOT/'data/pilot_manifest.csv').open()))
    assert len(pilot)==3 and sum(int(r['bytes']) for r in pilot)==514539339
    endpoints=list(csv.DictReader((ROOT/'data/endpoint_reference_freeze.csv').open()))
    counts=Counter(r['endpoint'] for r in endpoints if r['eligibility_before_visual_QC']=='True')
    assert counts=={'ER_binary':118,'PR_binary':112,'HER2_binary':95,'PAM50_primary_4class':114,'PIK3CA_binary':120}
    protocol=json.loads((ROOT/'endpoint_protocol.json').read_text())
    assert protocol['PAM50_class_order']=={'LumA':0,'LumB':1,'Basal':2,'Her2':3}
    assert protocol['files_sha256']['endpoint_reference_freeze.csv']==hashlib.sha256((ROOT/'data/endpoint_reference_freeze.csv').read_bytes()).hexdigest()
    # This prefix screen can identify potential duplicates; it does not prove complete-file uniqueness.
    prefixes=Counter(r['prefix_sha256'] for r in rows)
    result={'status':'VERIFIED','image_records_checked':len(rows),'patients':120,'endpoint_rows':len(endpoints),
        'metadata_node_corrections':corrected,'identical_remote_prefix_groups':sum(v>1 for v in prefixes.values()),
        'source_full_file_checksums_available':False,'image_decoding':'PENDING_LOCAL_PILOT',
        'external_features_or_metrics_computed':False}
    (ROOT/'delivery_audit.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
