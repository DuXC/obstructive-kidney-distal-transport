#!/usr/bin/env python3
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import gzip,io,json,hashlib
import pandas as pd
import numpy as np
P=Path(__file__).resolve().parents[1]
freeze=json.loads((P/'03_protocol/freeze_receipt.json').read_text())
assert hashlib.sha256((P/freeze['file']).read_bytes()).hexdigest()==freeze['sha256']
meta=pd.read_csv(P/'02_audit/GSE96102_sample_inventory.tsv',sep='\t')
meta=meta[meta.condition.ne('Normal')].sort_values('sample')

def extract(row):
    with gzip.open(P/'01_sources/raw'/row.filename,'rt') as f:
        for l in f:
            if l.startswith('FEATURES\t'):break
        x=pd.read_csv(io.StringIO(l+f.read()),sep='\t',low_memory=False)
    x=x[x.ControlType.eq(0)].copy()
    x['feature_key']=x.FeatureNum.astype(str)+'|'+x.ProbeName
    assert x.feature_key.is_unique
    return row['sample'],x.set_index('feature_key')[['FeatureNum','ProbeName','rMedianSignal','rBGMedianSignal']]
with ThreadPoolExecutor(max_workers=3) as pool: tables=dict(pool.map(extract,[r for _,r in meta.iterrows()]))
common=set.intersection(*(set(x.index) for x in tables.values()))
keys=sorted(common,key=lambda x:int(x.split('|')[0]))
fore=pd.DataFrame({s:d.reindex(keys).rMedianSignal for s,d in tables.items()})
back=pd.DataFrame({s:d.reindex(keys).rBGMedianSignal for s,d in tables.items()})
valid=np.isfinite(fore).all(axis=1)&np.isfinite(back).all(axis=1)
fore=fore.loc[valid];back=back.loc[valid]
for name,d in [('cy5_foreground',fore),('cy5_background',back)]:d.to_csv(P/'01_sources'/(name+'.tsv.gz'),sep='\t',index_label='feature_key',compression={'method':'gzip','mtime':0})
fmap=tables[next(iter(tables))].loc[fore.index,['FeatureNum','ProbeName']]
fmap.to_csv(P/'02_audit/retained_feature_map.tsv',sep='\t',index_label='feature_key')
qc={'samples':len(meta),'common_noncontrol_features':len(common),'finite_retained_features':len(fore),
    'sample_ids':meta['sample'].tolist(),'protocol_sha256':freeze['sha256'],
    'foreground_range':[float(fore.min().min()),float(fore.max().max())],
    'background_range':[float(back.min().min()),float(back.max().max())]}
(P/'09_qa/cy5_input_receipt.json').write_text(json.dumps(qc,indent=2))
print(json.dumps(qc,indent=2))
