"""Freeze explicit public-reference endpoint mappings before external prediction."""
import argparse
import ast
from collections import Counter
import csv
import hashlib
import json
from pathlib import Path

QUALIFYING_CLASSES={'Missense_Mutation','In_Frame_Del','In_Frame_Ins','Nonsense_Mutation',
    'Frame_Shift_Del','Frame_Shift_Ins','Splice_Site','Translation_Start_Site','Nonstop_Mutation'}
PAM50={'LumA':0,'LumB':1,'Basal':2,'Her2':3}
FIELDS={'ER_binary':'ER_UPDATED_CLINICAL_STATUS','PR_binary':'PR_CLINICAL_STATUS','HER2_binary':'ERBB2_UPDATED_CLINICAL_STATUS'}


def freeze(clinical,source,out,tcga_definition,tcga_pipeline):
    # Read the original source as data; do not execute a legacy script.
    tree=ast.parse(tcga_definition.read_text(encoding='utf-8-sig'))
    prior=next(ast.literal_eval(n.value) for n in tree.body if isinstance(n,ast.Assign)
               and any(isinstance(t,ast.Name) and t.id=='QUALIFYING_CLASSES' for t in n.targets))
    assert prior==QUALIFYING_CLASSES,'Variant-class definition differs from TCGA v2'
    with (tcga_pipeline/'data/stage10a_patient_endpoint_master_v1.csv').open(encoding='utf-8-sig',newline='') as f:
        tcga_names={r['patient_id']:r['pam50_primary_4class_label'] for r in csv.DictReader(f)}
    with (tcga_pipeline/'data/stage10a_internal_5fold_manifest_v1.csv').open(encoding='utf-8-sig',newline='') as f:
        tcga_rows=[r for r in csv.DictReader(f) if r['endpoint']=='PAM50_primary_4class']
    observed={}
    for r in tcga_rows:
        name=tcga_names[r['patient_id']].removeprefix('BRCA_');label=int(r['label'])
        if name in observed:assert observed[name]==label
        observed[name]=label
    assert observed==PAM50,'External PAM50 class order differs from the TCGA models'
    samples=json.loads((clinical/'samples.json').read_text())
    mutations=json.loads((clinical/'pik3ca_mutations.json').read_text())
    seq=next(x for x in json.loads((clinical/'sample_lists.json').read_text()) if x['sampleListId']=='brca_cptac_2020_sequenced')
    matrix=list(csv.DictReader((source/'cbio_gene_matrix_resolved.txt').open(),delimiter='\t'))
    seq_ids=set(seq['sampleIds']); matrix_ids={r['SAMPLE_ID'] for r in matrix if r['mutations']=='WXS'}
    sample_ids={r['sampleId'] for r in samples}
    assert seq_ids==matrix_ids==sample_ids and len(sample_ids)==122
    assert all(r.get('sequenced') is True and r['sampleType']=='Primary Solid Tumor' for r in samples)
    assert all(r['sampleId'] in sample_ids and r['entrezGeneId']==5290 and r['gene']['hugoGeneSymbol']=='PIK3CA' for r in mutations)
    qualifying=Counter(r['sampleId'] for r in mutations if r['mutationType'] in QUALIFYING_CLASSES)
    attrs={}
    for r in json.loads((clinical/'patient_clinical_data.json').read_text())+json.loads((clinical/'sample_clinical_data.json').read_text()):
        k=(r['patientId'],r['clinicalAttributeId'])
        if k in attrs:assert attrs[k]==r['value'],'Conflicting clinical values'
        attrs[k]=r['value']
    with (out/'data/patient_crosswalk.csv').open() as f:cross={r['cbio_patient_id']:r for r in csv.DictReader(f)}
    endpoint_rows=[]
    for s in samples:
        c=cross[s['patientId']];has_image=int(c['linked_tumor_slides'])>0
        for endpoint,field in FIELDS.items():
            raw=attrs.get((s['patientId'],field),'')
            label={'positive':1,'negative':0}.get(raw.lower())
            endpoint_rows.append({'patient_id':c['tcia_patient_id'],'cbio_patient_id':s['patientId'],'endpoint':endpoint,
                'source_field':field,'raw_reference':raw,'label':'' if label is None else label,
                'reference_evaluable':label is not None,'has_linked_tumor_image':has_image,
                'eligibility_before_visual_QC':has_image and label is not None})
        raw=attrs.get((s['patientId'],'PAM50'),'');label=PAM50.get(raw)
        endpoint_rows.append({'patient_id':c['tcia_patient_id'],'cbio_patient_id':s['patientId'],'endpoint':'PAM50_primary_4class',
            'source_field':'PAM50','raw_reference':raw,'label':'' if label is None else label,
            'reference_evaluable':label is not None,'has_linked_tumor_image':has_image,
            'eligibility_before_visual_QC':has_image and label is not None})
        positive=bool(qualifying[s['sampleId']])
        endpoint_rows.append({'patient_id':c['tcia_patient_id'],'cbio_patient_id':s['patientId'],'endpoint':'PIK3CA_binary',
            'source_field':'brca_cptac_2020_mutations + sequenced sample list + WXS matrix',
            'raw_reference':'QUALIFYING_REPORTED_MUTATION' if positive else 'NO_QUALIFYING_REPORTED_MUTATION_IN_WXS_PROFILE',
            'label':int(positive),'reference_evaluable':True,'has_linked_tumor_image':has_image,
            'eligibility_before_visual_QC':has_image})
    with (out/'data/endpoint_reference_freeze.csv').open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(endpoint_rows[0]));w.writeheader();w.writerows(endpoint_rows)
    summary={}
    for endpoint in list(FIELDS)+['PAM50_primary_4class','PIK3CA_binary']:
        selected=[r for r in endpoint_rows if r['endpoint']==endpoint and r['eligibility_before_visual_QC']]
        summary[endpoint]={'eligible_before_visual_QC':len(selected),'class_counts':dict(Counter(str(r['label']) for r in selected))}
    protocol={'status':'REFERENCE_DEFINITIONS_FROZEN_BEFORE_EXTERNAL_PREDICTION','summary':summary,
        'binary_class_order':{'0':'Negative or no qualifying reported PIK3CA mutation','1':'Positive or qualifying reported PIK3CA mutation'},
        'PAM50_class_order':PAM50,'PIK3CA_qualifying_classes':sorted(QUALIFYING_CLASSES),
        'PAM50_class_order_checked_against_original_TCGA_frozen_partitions':True,
        'PIK3CA_profiled_denominator':122,'PIK3CA_variants':len(mutations),'PIK3CA_qualifying_variants':sum(qualifying.values()),
        'clinical_HER2_proxy_substitution':False,'equivocal_or_missing_HER2_excluded':True,
        'external_recalibration_allowed':False,'external_outcome_model_selection_allowed':False,
        'image_QC_pending':True,'same_assayed_block_linkage':'Not demonstrated; analysis is linked at patient level',
        'PIK3CA_limitations':'Negatives mean no qualifying reported mutation in the WXS-profiled reference, not guaranteed base-level callability or clinical wild type. Upstream somatic filtering differs across source cohorts.',
        'files_sha256':{'endpoint_reference_freeze.csv':hashlib.sha256((out/'data/endpoint_reference_freeze.csv').read_bytes()).hexdigest(),
                        'TCGA_pik3ca_audit_v2.py':hashlib.sha256(tcga_definition.read_bytes()).hexdigest(),
                        'WXS_profile_matrix':hashlib.sha256((source/'cbio_gene_matrix_resolved.txt').read_bytes()).hexdigest()}}
    (out/'endpoint_protocol.json').write_text(json.dumps(protocol,indent=2))
    print(json.dumps(protocol['summary'],indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser()
    for key in ['clinical','source','out','tcga-definition','tcga-pipeline']:p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();freeze(a.clinical,a.source,a.out,a.tcga_definition,a.tcga_pipeline)
