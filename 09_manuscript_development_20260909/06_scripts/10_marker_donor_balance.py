#!/usr/bin/env python3
"""Pool technical libraries within donor before equal-donor marker summaries."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
DEV=Path(__file__).resolve().parents[1];ROOT=DEV.parent
m=pd.read_csv(ROOT/'05_results/marker_audit_by_donor.tsv.gz',sep='\t')
d=m[m.level.eq('celltype')].copy()
d['detected_sum']=d.fraction_detected*d.n_nuclei
d['expression_sum']=d.mean_log1p_cp10k*d.n_nuclei
donors=d.groupby(['donor','group','annotation','gene'],as_index=False).agg(n_nuclei=('n_nuclei','sum'),detected_sum=('detected_sum','sum'),expression_sum=('expression_sum','sum'),technical_records=('n_nuclei','size'))
donors['fraction_detected']=donors.detected_sum/donors.n_nuclei
donors['mean_log1p_cp10k']=donors.expression_sum/donors.n_nuclei
assert not donors.duplicated(['donor','annotation','gene']).any()
assert np.allclose(donors.detected_sum,donors.detected_sum.round(),atol=1e-6)
assert donors.fraction_detected.between(0,1+1e-12).all()
qc=pd.read_csv(ROOT/'05_results/donor_celltype_qc.tsv',sep='\t').set_index(['donor','celltype'])
assert all(row.n_nuclei==qc.loc[(row.donor,row.annotation),'n_nuclei'] for row in donors.itertuples())
balanced=donors.groupby(['group','annotation','gene'],as_index=False).agg(mean_log1p_cp10k=('mean_log1p_cp10k','mean'),mean_fraction_detected=('fraction_detected','mean'),n_donors=('donor','nunique'))
donors.to_csv(DEV/'02_analysis/markers_pooled_within_donor.tsv',sep='\t',index=False)
balanced.to_csv(DEV/'02_analysis/marker_donor_balanced.tsv',sep='\t',index=False)
(DEV/'07_qa/marker_balance_receipt.json').write_text(json.dumps({'status':'PASS','pooled_donor_gene_records':len(donors),'equal_donor_summary_records':len(balanced),'method':'nucleus-weighted pooling across technical libraries within each donor, followed by equal-donor means','matches_donor_nucleus_counts':True,'primary_model_inputs_changed':False},indent=2))
print('PASS: pooled technical libraries within donor before equal-donor summaries')
