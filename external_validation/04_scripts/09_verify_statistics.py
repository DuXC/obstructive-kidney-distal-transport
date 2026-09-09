from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from scipy import stats
P=Path(__file__).resolve().parents[1]
def bh(p):
 p=np.asarray(p,float);out=np.full(len(p),np.nan);ok=np.isfinite(p);order=np.argsort(p[ok]);v=p[ok][order]
 q=np.minimum.accumulate((v*10/np.arange(1,len(v)+1))[::-1])[::-1].clip(0,1);back=np.empty_like(q);back[order]=q;out[ok]=back;return out
audits=[]
for assay in ['sn','sc']:
 d=pd.read_csv(P/f'05_results/{assay}_program_comparisons.tsv',sep='\t');s=pd.read_csv(P/f'05_results/{assay}_donor_scores.tsv',sep='\t')
 primary=d[d.scenario=='primary'].copy();assert len(primary)==10
 for _,r in primary.iterrows():
  z=s[(s.scenario=='primary')&(s.module==r.module)&(s.contrast==r.contrast)]
  if r.status!='TESTED':continue
  assert not z.patient.duplicated().any()
  a=z[z.group==r.contrast].score.to_numpy();b=z[z.group=='Ref'].score.to_numpy()
  assert len(a)==r.n_disease and len(b)==r.n_reference
  se=np.sqrt(np.var(a,ddof=1)/len(a)+np.var(b,ddof=1)/len(b))
  df=(np.var(a,ddof=1)/len(a)+np.var(b,ddof=1)/len(b))**2/((np.var(a,ddof=1)/len(a))**2/(len(a)-1)+(np.var(b,ddof=1)/len(b))**2/(len(b)-1))
  effect=a.mean()-b.mean();ci=np.asarray([effect-stats.t.ppf(.975,df)*se,effect+stats.t.ppf(.975,df)*se])
  tt=stats.ttest_ind(a,b,equal_var=False)
  assert np.allclose([effect,*ci,tt.pvalue],[r.effect,r.ci_low,r.ci_high,r.p_value],atol=1e-9,rtol=1e-8)
 assert np.allclose(bh(primary.p_value),primary.fdr_10,equal_nan=True)
 coverage=pd.read_csv(P/f'05_results/{assay}_member_coverage.tsv',sep='\t');assert coverage.represented.all()
 lodo=d[d.scenario.str.startswith('LODO_')]
 summary=lodo.groupby(['contrast','module']).agg(n_omissions=('effect','size'),minimum=('effect','min'),maximum=('effect','max'),n_negative=('effect',lambda x:(x<0).sum()),n_positive=('effect',lambda x:(x>0).sum())).reset_index()
 summary.to_csv(P/f'05_results/{assay}_lodo_summary.tsv',sep='\t',index=False)
 audits.append({'assay':assay,'primary_tests':10,'primary_tested':int((primary.status=='TESTED').sum()),'fdr_lt_05':int((primary.fdr_10<.05).sum()),'welch_ci_bh_independently_recomputed':True,'all_fixed_members_represented':True})
 print(assay,primary[['contrast','module','n_reference','n_disease','effect','ci_low','ci_high','fdr_10','status']].to_string(index=False))
receipt=json.loads((P/'03_protocol/freeze_receipt.json').read_text())
for name,value in receipt['files'].items():assert hashlib.sha256((P/'03_protocol'/name).read_bytes()).hexdigest()==value
(P/'08_qa/independent_statistics_check.json').write_text(json.dumps({'status':'PASS','frozen_plan_and_memberships_unchanged':True,'assays':audits},indent=2))

