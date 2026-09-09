#!/usr/bin/env python3
"""Donor-level composition/covariation and independent checks of primary statistics."""
from pathlib import Path
from itertools import combinations
import json,hashlib
import numpy as np,pandas as pd
from scipy import stats
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'05_results'
def bh(p):
    p=np.asarray(p,float);out=np.full(len(p),np.nan);k=np.flatnonzero(np.isfinite(p));o=k[np.argsort(p[k])]
    if len(o):out[o]=np.minimum(1,np.minimum.accumulate((p[o]*len(o)/np.arange(1,len(o)+1))[::-1])[::-1])
    return out
def compare(y,g):
    a=y[g=='UUO'];b=y[g=='Control'];tt=stats.ttest_ind(a,b,equal_var=False);ci=tt.confidence_interval()
    obs=a.mean()-b.mean();den=len(y)-len(a)
    vals=np.array([(y[list(c)].sum()/len(a)-(y.sum()-y[list(c)].sum())/den) for c in combinations(range(len(y)),len(a))])
    return dict(effect=obs,ci_low=ci.low,ci_high=ci.high,p_welch=tt.pvalue,p_permutation=float((np.abs(vals)>=abs(obs)-1e-12).mean()))
cmp=pd.read_csv(R/'module_comparisons_primary_run.tsv',sep='\t');sc=pd.read_csv(R/'module_scores_all_scenarios.tsv',sep='\t');main=sc[sc.scenario.eq('main')]
checks=[]
for _,r in cmp[cmp.family.eq('primary_transport')].iterrows():
    ss=main[(main.celltype==r.celltype)&(main.module==r.module)]
    z=pd.read_csv(ROOT/f'02_data/derived/logcpm_{r.celltype}.tsv.gz',sep='\t',index_col=0)
    mm=next(m for m in json.loads((ROOT/'03_protocol/modules_v1.json').read_text()) if m['module']==r.module)
    independent=z.loc[mm['genes'],ss.donor+'__'+r.celltype].mean(axis=0).to_numpy()
    assert np.allclose(independent,ss.score_mean_log2cpm,atol=1e-9)
    rr=compare(independent,ss.group.to_numpy())
    for a,b in [('effect','effect'),('ci_low','ci_low'),('ci_high','ci_high'),('p_welch','p_value'),('p_permutation','permutation_p')]:assert abs(rr[a]-r[b])<1e-8,(r.module,a,rr[a],r[b])
    checks.append(dict(module=r.module,status='PASS_independent_score_Welch_CI_exact_permutation'))
sel=cmp.family.eq('primary_transport');assert np.allclose(bh(cmp.loc[sel,'p_value']),cmp.loc[sel,'fdr'])

co=pd.read_csv(R/'captured_composition.tsv',sep='\t');cr=[]
for den in ['fraction_all_captured','fraction_captured_epithelia']:
    for ct in ['TAL','DCT','PC','IC-A','IC-B']:
        d=co[co.celltype.eq(ct)].sort_values('donor');rr=compare(d[den].to_numpy()*100,d.group.to_numpy())
        cr.append(dict(celltype=ct,denominator=den,n_control=sum(d.group=='Control'),n_uuo=sum(d.group=='UUO'),**rr))
cr=pd.DataFrame(cr);cr['fdr_permutation']=cr.groupby('denominator').p_permutation.transform(bh)
cr.to_csv(R/'captured_composition_comparisons.tsv',sep='\t',index=False)

ar=[]
for _,r in cmp[cmp.family.str.contains('transport')].iterrows():
    ss=main[(main.celltype==r.celltype)&(main.module==r.module)][['donor','group','score_mean_log2cpm']].rename(columns={'score_mean_log2cpm':'transport'})
    for ctxt in ['epithelial_injury','chemokine_adhesion']:
        d=main[(main.celltype==r.celltype)&(main.module==ctxt)][['donor','score_mean_log2cpm']].merge(ss,on='donor',validate='one_to_one')
        n=len(d);g=(d.group=='UUO').astype(float).to_numpy()
        if n<8 or min(sum(g==0),sum(g==1))<3:continue
        x=d.transport.to_numpy();y=d.score_mean_log2cpm.to_numpy();design=np.column_stack([np.ones(n),g])
        rx=x-design@np.linalg.lstsq(design,x,rcond=None)[0];ry=y-design@np.linalg.lstsq(design,y,rcond=None)[0]
        raw=stats.pearsonr(x,y);adj=np.corrcoef(rx,ry)[0,1];df=n-3;t=adj*np.sqrt(df/max(1-adj**2,1e-15));pv=2*stats.t.sf(abs(t),df)
        z=np.arctanh(np.clip(adj,-.999999,.999999));lo,hi=np.tanh(z+np.array([-1,1])*stats.norm.ppf(.975)/np.sqrt(n-4))
        ar.append(dict(celltype=r.celltype,module=r.module,context=ctxt,n=n,n_control=int(sum(g==0)),n_uuo=int(sum(g==1)),pooled_r=raw.statistic,pooled_p=raw.pvalue,partial_r=adj,partial_ci_low=lo,partial_ci_high=hi,partial_p=pv,ci_method='approximate Fisher conditional on group'))
ar=pd.DataFrame(ar);ar['partial_fdr']=bh(ar.partial_p);ar.to_csv(R/'transport_context_associations.tsv',sep='\t',index=False)

alls=pd.read_csv(R/'module_comparisons_all_scenarios.tsv',sep='\t');sens=[]
for _,r in cmp[cmp.family.str.contains('transport')].iterrows():
    d=alls[(alls.module==r.module)&(alls.celltype==r.celltype)&alls.scenario.str.startswith('LODO')]
    sens.append(dict(celltype=r.celltype,module=r.module,main_effect=r.effect,n_lodo=len(d),n_same_direction=int((np.sign(d.effect)==np.sign(r.effect)).sum()),lodo_min=d.effect.min(),lodo_max=d.effect.max(),n_lodo_nominal_p05=int((d.p_value<.05).sum()),interpretation='direction stability; significance is not required in every reduced sample'))
pd.DataFrame(sens).to_csv(R/'module_lodo_summary.tsv',sep='\t',index=False)

# Equal-donor marker summaries are produced by the subsequent donor-marker step.
rawq=pd.read_csv(ROOT/'02_data/derived/nucleus_qc.tsv.gz',sep='\t')
qcsummary=rawq.groupby(['donor','group','celltype'],as_index=False).agg(n_nuclei=('barcode','size'),median_umis=('n_counts','median'),median_genes=('n_genes','median'),median_mito_percent=('percent_mito','median'))
qcsummary.to_csv(R/'donor_celltype_rna_qc.tsv',sep='\t',index=False)
freeze=json.loads((ROOT/'03_protocol/freeze_receipt.json').read_text());assert all(hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h for p,h in freeze['files'].items())
summary=dict(status='PASS',independent_primary_statistic_checks=checks,primary_BH_verified=True,frozen_files_unchanged=True,composition_tests=len(cr),partial_correlations=len(ar),partial_correlations_fdr05=int((ar.partial_fdr<.05).sum()))
(ROOT/'08_logs/statistical_verification.json').write_text(json.dumps(summary,indent=2));print(json.dumps(summary,indent=2))
print(ar[['celltype','module','context','pooled_r','partial_r','partial_fdr']].to_string(index=False))
