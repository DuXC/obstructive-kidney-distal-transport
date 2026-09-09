from pathlib import Path
import pandas as pd,json,gzip,re
P=Path(__file__).resolve().parents[1]
sn=pd.read_csv(next((P/'01_sources').glob('GSE183277*Metadata_0328*')),sep='\t',index_col=0,low_memory=False)
sc=pd.read_csv(next((P/'01_sources').glob('GSE183276*Metadata_0328*')),sep='\t',index_col=0,low_memory=False)
for d in [sn,sc]:d['patient']=d['patient'].astype(str).str.strip()
rows=[]
for label,d in [('sn',sn),('sc',sc)]:
 print(label,'cells',len(d),'unique_barcodes',d.index.nunique(),'donors',d.patient.nunique())
 fields=['patient','condition.l1','condition.l2','condition.l3','region.l1','region.l2','percent.cortex','tissue_type','sex']
 samp=d[fields].drop_duplicates().sort_values(['condition.l1','patient'])
 samp.to_csv(P/'05_results'/f'{label}_metadata_strata.tsv',sep='\t',index=False)
 print(samp.to_string(index=False) if label=='sn' else samp.groupby(['condition.l1','condition.l3','tissue_type'],dropna=False).patient.nunique().to_string())
 print('all_lineage counts',d['subclass.l1'].value_counts().to_dict())
 print('state names',d['state.l2'].value_counts().to_dict())
 print('eligible cortex-containing nonstone counts')
 sel=d.copy()
 if label=='sn':sel=d[(pd.to_numeric(d['percent.cortex'],errors='coerce')>0)&(d['condition.l3']!='Stone')]
 for ct in ['TAL','DCT']:
  n=sel[sel['subclass.l1']==ct].groupby(['condition.l1','patient']).size()
  print(ct, n[n>=20].groupby(level=0).size().to_dict())
 print('reference patient list',sorted(d[d['condition.l1']=='Ref'].patient.unique()))
overlap=sorted(set(sn.patient)&set(sc.patient))
print('overlap sn/sc',len(overlap),overlap)
print('sc nonoverlap donors',sc[~sc.patient.isin(overlap)].groupby('condition.l1').patient.nunique().to_dict())
(P/'05_results/metadata_overlap.json').write_text(json.dumps({'overlap_count':len(overlap),'overlap_patients':overlap},indent=2))

