#!/usr/bin/env python3
"""Acquire source feature tables and preserve probe identity, quality flags and batches."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import json, hashlib, urllib.request, gzip, io, time
import pandas as pd

P=Path(__file__).resolve().parents[1]
meta=pd.read_csv(P/'02_audit/GSE96102_sample_inventory.tsv',sep='\t')
expected_raw=pd.read_csv(P/'raw_data_sources.tsv',sep='\t').set_index('sample')
ann=pd.read_csv(P/'01_sources/GPL4134_annotation.tsv',sep='\t',dtype=str).fillna('')
lookup=json.loads((P/'01_sources/Ensembl_mouse_lookup.json').read_text())
hom=pd.read_csv(P/'02_audit/Ensembl_orthology.tsv',sep='\t')
mods=[m for m in json.loads((P.parent/'03_protocol/modules_v1.json').read_text()) if m['family']=='primary_transport']
genes={g for m in mods for g in m['genes']}
hom=hom[hom.human_symbol.isin(genes)].copy()
hom['mouse_symbol']=hom.mouse_ensembl.map(lambda g:lookup[g]['display_name'])
hom.to_csv(P/'02_audit/primary_orthology.tsv',sep='\t',index=False)
unique=ann.groupby('SPOT_ID').GENE_SYMBOL.agg(lambda x:sorted(set(x)-{''}))
valid_names={k:v[0] for k,v in unique.items() if len(v)==1}
target_names={k:v for k,v in valid_names.items() if v in set(hom.mouse_symbol)}
pd.DataFrame([{'probe':k,'mouse_symbol':v} for k,v in target_names.items()]).to_csv(P/'02_audit/target_probe_map.tsv',sep='\t',index=False)

def acquire(row):
    dest=P/'01_sources/raw'/row.filename
    dest.parent.mkdir(exist_ok=True)
    original=P/'01_sources'/row.filename
    if original.exists() and not dest.exists(): original.rename(dest)
    if not dest.exists():
        for attempt in range(3):
            try:
                with urllib.request.urlopen(row.source_url,timeout=50) as r: body=r.read()
                assert body[:2]==b'\x1f\x8b'
                dest.write_bytes(body);break
            except Exception:
                if attempt==2: raise
                time.sleep(2)
    payload=dest.read_bytes();digest=hashlib.sha256(payload).hexdigest()
    expected=expected_raw.loc[row['sample']]
    if len(payload)!=expected['bytes'] or digest!=expected.sha256:
        raise ValueError('Source bytes or checksum changed: '+dest.name)
    with gzip.open(dest,'rt') as f:
        lines=[]; params={}
        for l in f:
            lines.append(l)
            if l.startswith('FEPARAMS\t'): keys=l.rstrip('\n').split('\t')
            elif l.startswith('DATA\t') and not params:params=dict(zip(keys,l.rstrip('\n').split('\t')))
            if l.startswith('FEATURES\t'):break
        dat=pd.read_csv(io.StringIO(l+f.read()),sep='\t',low_memory=False)
    assert dat.FeatureNum.is_unique
    tar=dat[dat.ProbeName.isin(target_names)].copy()
    tar['mouse_symbol']=tar.ProbeName.map(target_names)
    tar.insert(0,'sample',row['sample'])
    keep=['sample','FeatureNum','ProbeName','mouse_symbol','ControlType','GeneName','SystematicName',
          'LogRatio','gSurrogateUsed','rSurrogateUsed','gIsFound','rIsFound',
          'gIsWellAboveBG','rIsWellAboveBG','gIsSaturated','rIsSaturated',
          'gIsFeatNonUnifOL','rIsFeatNonUnifOL','gIsFeatPopnOL','rIsFeatPopnOL','IsManualFlag']
    optional_missing=[k for k in ['GeneName','SystematicName','gSurrogateUsed','rSurrogateUsed','gIsFound','rIsFound'] if k not in tar]
    for k in optional_missing: tar[k]=float('nan')
    for k in keep: assert k in tar,k
    tar=tar[keep]
    stat={'sample':row['sample'],'file':dest.name,'source_url':row.source_url,
          'bytes':len(payload),'sha256':digest,
          'scan_date':params.get('Scan_Date'),'array_barcode':params.get('FeatureExtractor_Barcode'),
          'design':params.get('FeatureExtractor_DesignFileName'),'extraction_version':params.get('FeatureExtractor_Version'),
          'optional_annotation_fields_absent':optional_missing,
          'features':len(dat),'noncontrol_features':int(dat.ControlType.eq(0).sum()),
          'target_features':len(tar),'target_genes':tar.mouse_symbol.nunique(),
          'r_above_background_fraction':float(dat.loc[dat.ControlType.eq(0),'rIsWellAboveBG'].mean()),
          'g_above_background_fraction':float(dat.loc[dat.ControlType.eq(0),'gIsWellAboveBG'].mean())}
    return tar,stat

results=[]; receipts=[]
with ThreadPoolExecutor(max_workers=3) as pool:
    fs={pool.submit(acquire,r):r['sample'] for _,r in meta.iterrows()}
    for f in as_completed(fs):
        tar,rec=f.result();results.append(tar);receipts.append(rec)
        if len(receipts)%5==0: print('Verified',len(receipts),'of',len(meta),'raw feature tables',flush=True)
table=pd.concat(results).sort_values(['sample','FeatureNum'])
table.to_csv(P/'01_sources/target_feature_values.tsv',sep='\t',index=False)
receipts.sort(key=lambda r:r['sample'])
(P/'09_qa/raw_acquisition_receipt.json').write_text(json.dumps(receipts,indent=2))
batches=meta.merge(pd.DataFrame(receipts),on='sample',validate='one_to_one')
batches.to_csv(P/'02_audit/sample_batches.tsv',sep='\t',index=False)
print('Total bytes',sum(r['bytes'] for r in receipts))
print(batches.groupby(['duration_days','recovery_days','condition','scan_date']).size().to_string())
print('All target symbols',sorted(table.mouse_symbol.unique()))
