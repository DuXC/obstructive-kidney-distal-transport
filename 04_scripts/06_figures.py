#!/usr/bin/env python3
"""Generate four source-linked, editable quantitative research figures."""
from pathlib import Path
import json
import numpy as np,pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm
from matplotlib.lines import Line2D
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'05_results';O=ROOT/'06_figures'
plt.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7,'axes.titlesize':8,'axes.labelsize':7,'xtick.labelsize':6.5,'ytick.labelsize':6.5,'axes.spines.top':False,'axes.spines.right':False,'axes.linewidth':.6,'lines.linewidth':1,'svg.fonttype':'none','pdf.fonttype':42,'legend.frameon':False,'savefig.facecolor':'white'})
COL={'Control':'#356F85','UUO':'#C16B3A'};DON=[f'Control{x}' for x in range(1,8)]+[f'UUO{x}' for x in range(1,6)]
CT=['TAL','DCT','PC','IC-A','IC-B'];PRIMARY=['TAL_salt','TAL_CaMg','DCT_NaCl','DCT_Mg','DCT_Ca']
SHORT={'TAL_salt':'TAL salt','TAL_CaMg':'TAL Ca/Mg','DCT_NaCl':'DCT NaCl','DCT_Mg':'DCT Mg','DCT_Ca':'DCT Ca','PC_water':'PC water','PC_NaK':'PC Na/K','ICA_acid':'IC-A acid','ICB_base':'IC-B base'}
cmp=pd.read_csv(R/'module_comparisons_primary_run.tsv',sep='\t');scores=pd.read_csv(R/'module_scores_all_scenarios.tsv',sep='\t');main=scores[scores.scenario.eq('main')];allcmp=pd.read_csv(R/'module_comparisons_all_scenarios.tsv',sep='\t')
def panel(ax,letter,title):
    ax.set_title(title,loc='left',pad=9,fontweight='normal');ax.text(-.14,1.075,letter,transform=ax.transAxes,fontweight='bold',fontsize=10,va='bottom')
def legend(ax,loc='best'):
    ax.legend(handles=[Line2D([],[],marker='o',linestyle='',color=COL[g],label=('Unobstructed reference' if g=='Control' else 'UTUC + obstruction'),markersize=4) for g in COL],loc=loc,fontsize=6)
def save(fig,name):
    for ext in ['svg','pdf','png']:fig.savefig(O/f'{name}.{ext}',dpi=301 if ext=='png' else 300)
    plt.close(fig)
def donorplot(ax,frame,categories,value,labelmap=None):
    for i,cat in enumerate(categories):
        for gi,g in enumerate(['Control','UUO']):
            d=frame[(frame.category==cat)&(frame.group==g)].sort_values('donor');y=d[value].to_numpy()
            off=(-.17,.17)[gi];j=np.linspace(-.07,.07,len(y))
            ax.scatter(i+off+j,y,s=15,color=COL[g],edgecolor='white',linewidth=.25,zorder=3)
            if len(y):ax.plot([i+off-.10,i+off+.10],[y.mean()]*2,color=COL[g],lw=1.6,zorder=4)
    ax.set_xticks(range(len(categories)),[(labelmap or {}).get(x,x) for x in categories]);ax.set_xlim(-.55,len(categories)-.45);ax.grid(axis='y',alpha=.15)
def effectforest(ax,frame,order,showq=True):
    for i,mo in enumerate(order):
        rr=frame[frame.module.eq(mo)].iloc[0];tested=pd.notna(rr.ci_low)
        if tested:
            ax.plot([rr.ci_low,rr.ci_high],[i,i],color='#536878',lw=1.2)
            ax.scatter(rr.effect,i,s=23,facecolor='#356F85' if rr.fdr<.05 else 'white',edgecolor='#536878',zorder=3)
        else:ax.scatter(rr.effect,i,marker='x',color='#888888',s=20)
        if showq:ax.text(1.015,i,f'{rr.fdr:.3f}' if tested else 'n.a.',transform=ax.get_yaxis_transform(),va='center',fontsize=6)
    ax.axvline(0,color='#aaaaaa',lw=.7,ls=':');ax.set_yticks(range(len(order)),[SHORT[x] for x in order]);ax.set_ylim(len(order)-.5,-.6);ax.set_xlabel('Difference in mean log2 CPM')
    if showq:ax.text(1.015,1.01,'FDR',transform=ax.transAxes,fontsize=6)

# Figure 1: sampling and author-lineage audit.
fig=plt.figure(figsize=(183/25.4,162/25.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1.2,1],width_ratios=[1.1,1])
ax=fig.add_subplot(gs[0,0]);lib=pd.read_csv(R/'donor_library_celltype_counts.tsv',sep='\t').groupby(['donor','Library']).n_nuclei.sum().unstack(fill_value=0).reindex(DON);frac=lib.div(lib.sum(axis=1),axis=0)*100
im=ax.imshow(frac,cmap='Blues',vmin=0,vmax=100,aspect='auto')
for i in range(12):
    for j in range(6):
        if frac.iloc[i,j]>0:ax.text(j,i,f'{frac.iloc[i,j]:.0f}',ha='center',va='center',fontsize=6,color='white' if frac.iloc[i,j]>60 else '#222222')
ax.set_xticks(range(6),[f'L{x}' for x in range(1,7)]);ax.set_yticks(range(12),DON);panel(ax,'a','Donor contribution by library (%)');fig.colorbar(im,ax=ax,fraction=.035,pad=.02)
ax=fig.add_subplot(gs[0,1]);cnt=pd.read_csv(R/'donor_celltype_qc.tsv',sep='\t').pivot(index='donor',columns='celltype',values='n_nuclei').reindex(index=DON,columns=CT)
im=ax.imshow(np.log10(cnt),cmap='Blues',vmin=0,vmax=3.3,aspect='auto')
for i in range(12):
    for j in range(5):ax.text(j,i,str(int(cnt.iloc[i,j])),ha='center',va='center',fontsize=6,color='white' if cnt.iloc[i,j]>260 else '#222222',fontweight='bold' if cnt.iloc[i,j]<20 else 'normal')
ax.set_xticks(range(5),CT);ax.set_yticks(range(12),DON);panel(ax,'b','Nuclei per donor and lineage');ax.text(0,-.16,'Bold counts: below the 20-nucleus threshold',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,0]);mk=pd.read_csv(R/'marker_audit_donor_balanced.tsv',sep='\t');mk=mk[mk.group.eq('Control')];cts=['TAL','DCT','CNT','PC','IC-A','IC-B','PT','Myeloid Cell'];markers=['UMOD','SLC12A1','SLC12A3','PVALB','CALB1','AQP2','AQP3','SLC4A1','ATP6V1B1','SLC26A4','LRP2','LST1']
for i,ct in enumerate(cts):
    for j,g in enumerate(markers):
        row=mk[(mk.annotation==ct)&(mk.gene==g)].iloc[0]
        im=ax.scatter(j,i,s=5+row.mean_fraction_detected*60,c=row.mean_log1p_cp10k,cmap='Blues',vmin=0,vmax=4,edgecolor='#8193A1',linewidth=.2)
ax.set_xticks(range(len(markers)),markers,rotation=65,ha='right');ax.set_yticks(range(len(cts)),cts);ax.set_ylim(len(cts)-.5,-.5);ax.set_xlim(-.7,len(markers)-.3);panel(ax,'c','Classical markers in reference donors');fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label='Mean log1p CP10K')
ax.text(0,-.46,'Size: donor-mean detection fraction (0-100%)',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,1]);qc=pd.read_csv(R/'donor_celltype_rna_qc.tsv',sep='\t');qc['category']=qc.celltype;donorplot(ax,qc,CT,'median_umis');ax.set_ylabel('Median RNA UMIs per nucleus');panel(ax,'d','RNA quality across donors');legend(ax,loc='upper left')
save(fig,'Figure1_donor_QC')

# Figure 2: capture/state composition and all primary fixed modules.
fig=plt.figure(figsize=(183/25.4,170/25.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.2],width_ratios=[1,1.15])
ax=fig.add_subplot(gs[0,0]);co=pd.read_csv(R/'captured_composition.tsv',sep='\t');co['category']=co.celltype;co['percent']=co.fraction_captured_epithelia*100;donorplot(ax,co,['TAL','DCT'],'percent');ax.set_ylabel('% of captured epithelial nuclei');panel(ax,'a','Relative capture composition');legend(ax)
ax=fig.add_subplot(gs[0,1]);st=pd.read_csv(R/'author_state_composition.tsv',sep='\t');st['injured']=st.state.str.contains('Injured|Inflammatory');st=st[st.celltype.isin(['TAL','DCT'])]
ss=st[st.injured].groupby(['donor','group','celltype']).fraction_within_celltype.sum().rename('fraction').reset_index();base=st[['donor','group','celltype']].drop_duplicates();ss=base.merge(ss,how='left').fillna({'fraction':0});ss['percent']=ss.fraction*100;ss['category']=ss.celltype;donorplot(ax,ss,['TAL','DCT'],'percent');ax.set_ylim(-3,103);ax.set_ylabel('% of nuclei within lineage');panel(ax,'b','Author-labelled injury states');ax.text(0,-.22,'Injured + inflammatory states; descriptive annotation',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,0]);effectforest(ax,cmp,PRIMARY);panel(ax,'c','Prespecified transport effects');ax.text(0,-.25,'95% Welch intervals; FDR across five tests',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,1]);d=main[main.module.isin(PRIMARY)].copy();d['category']=d.module;donorplot(ax,d,PRIMARY,'score_mean_log2cpm',SHORT);ax.tick_params(axis='x',labelrotation=35);ax.set_ylabel('Module mean log2 CPM');panel(ax,'d','Individual donor expression');ax.text(0,-.31,'TAL: 7 reference / 5 UUO; DCT: 6 reference / 5 UUO',transform=ax.transAxes,fontsize=6)
save(fig,'Figure2_primary_transport')

# Figure 3: all fixed members in the two primary lineages, collecting-duct exploration.
fig=plt.figure(figsize=(183/25.4,183/25.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1.5,1])
gene=pd.read_csv(R/'fixed_transport_gene_results.tsv',sep='\t');mods=json.loads((ROOT/'03_protocol/modules_v1.json').read_text())
for col,ct in enumerate(['TAL','DCT']):
    ax=fig.add_subplot(gs[0,col]);genes=list(dict.fromkeys(g for mm in mods if mm['celltype']==ct and mm['family']=='primary_transport' for g in mm['genes']))
    for i,g in enumerate(genes):
        z=gene[(gene.celltype==ct)&(gene.gene==g)]
        if len(z):
            r=z.iloc[0];ax.plot([r['CI.L'],r['CI.R']],[i,i],color='#647889');ax.scatter(r.logFC,i,s=19,facecolor='#356F85' if r.FDR_family<.05 else 'white',edgecolor='#536878',zorder=3)
        else:ax.text(.03,i,'Not tested: low counts',transform=ax.get_yaxis_transform(),fontsize=5.5,va='center',color='#777777')
    ax.set_yticks(range(len(genes)),genes);ax.set_ylim(len(genes)-.5,-.6);ax.axvline(0,color='#aaa',ls=':',lw=.7);ax.set_xlabel('Gene log2 fold change (UUO - reference)');panel(ax,'ab'[col],f'{ct}: prespecified member genes')
    ax.text(0,-.17,'Filled: FDR <0.05 across TAL + DCT gene tests',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,0]);extra=['PC_water','PC_NaK','ICA_acid','ICB_base'];effectforest(ax,cmp,extra);panel(ax,'c','Collecting-duct program estimates');ax.text(0,-.27,'PC 7/4; IC-A 6/5; IC-B 5/3 (reference / UUO)\nIC-B: descriptive point; no group test',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[1,1]);d=main[main.module.isin(extra)].copy();d['category']=d.module;donorplot(ax,d,extra,'score_mean_log2cpm',SHORT);ax.tick_params(axis='x',labelrotation=35);ax.set_ylabel('Module mean log2 CPM');panel(ax,'d','Collecting-duct donor expression');legend(ax,loc='lower right')
save(fig,'Figure3_genes_collecting_duct')

# Figure 4: stability and two prespecified primary-lineage salt/injury comparisons.
fig=plt.figure(figsize=(183/25.4,170/25.4),layout='constrained');gs=fig.add_gridspec(2,2,height_ratios=[1,1.2],width_ratios=[1,1.05])
ax=fig.add_subplot(gs[0,0]);ld=pd.read_csv(R/'module_lodo_summary.tsv',sep='\t').set_index('module')
for i,m in enumerate(PRIMARY):
    z=ld.loc[m];ax.plot([z.lodo_min,z.lodo_max],[i,i],lw=3,color='#a1b5c2');ax.scatter(z.main_effect,i,s=25,color='#356F85',zorder=4);ax.text(1.02,i,f'{int(z.n_same_direction)}/{int(z.n_lodo)}',transform=ax.get_yaxis_transform(),va='center',fontsize=6)
ax.set_yticks(range(5),[SHORT[x] for x in PRIMARY]);ax.set_ylim(4.6,-.6);ax.axvline(0,color='#aaa',ls=':',lw=.7);ax.set_xlabel('Difference in mean log2 CPM');panel(ax,'a','Leave-one-donor-out estimates');ax.text(0,-.24,'Bar: estimate range; dot: full-sample estimate',transform=ax.transAxes,fontsize=6)
ax.text(1.02,1.01,'Same sign',transform=ax.transAxes,fontsize=6)
ax=fig.add_subplot(gs[0,1]);scenarios=['main','min10','min50','sex_adjusted','sex_ageband_adjusted','author_healthy_states'];labs=['Main','Min 10','Min 50','+ sex','+ sex/age','Healthy']
df=allcmp[allcmp.module.isin(PRIMARY)].pivot_table(index='module',columns='scenario',values='effect').reindex(index=PRIMARY,columns=scenarios)
im=ax.imshow(df.to_numpy(),cmap='RdBu_r',vmin=-2,vmax=2,aspect='auto')
for i in range(5):
    for j in range(6):ax.text(j,i,f'{df.iloc[i,j]:.2f}' if pd.notna(df.iloc[i,j]) else 'n.a.',ha='center',va='center',fontsize=5.5,color='#222222')
ax.set_yticks(range(5),[SHORT[x] for x in PRIMARY]);ax.set_xticks(range(6),labs,rotation=40,ha='right');panel(ax,'b','Prespecified sensitivity estimates');fig.colorbar(im,ax=ax,fraction=.035,pad=.02,label='Effect');ax.text(0,-.46,'Healthy-state DCT: only 1 eligible UUO donor',transform=ax.transAxes,fontsize=6)
assoc=pd.read_csv(R/'transport_context_associations.tsv',sep='\t')
for col,(ct,mo) in enumerate([('TAL','TAL_salt'),('DCT','DCT_NaCl')]):
    ax=fig.add_subplot(gs[1,col]);d=main[(main.celltype==ct)&main.module.isin([mo,'epithelial_injury'])].pivot(index=['donor','group'],columns='module',values='score_mean_log2cpm').reset_index()
    for g in COL:
        s=d[d.group==g];ax.scatter(s.epithelial_injury,s[mo],color=COL[g],s=22,edgecolor='white',linewidth=.35,zorder=3)
        for _,rr in s.iterrows():
            offset=(-10,-8) if ct=='TAL' and rr.donor=='Control7' else (3,3)
            ax.annotate(rr.donor.replace('Control','C').replace('UUO','U'),(rr.epithelial_injury,rr[mo]),xytext=offset,textcoords='offset points',fontsize=5,color=COL[g])
    rr=assoc[(assoc.module==mo)&assoc.context.eq('epithelial_injury')].iloc[0]
    ax.set_xlabel('Epithelial injury score (mean log2 CPM)');ax.set_ylabel(SHORT[mo]+' score (mean log2 CPM)');panel(ax,'cd'[col],f'{ct} transport and epithelial injury')
    ax.text(0,-.26,f'Pooled r = {rr.pooled_r:.2f}; group-adjusted r = {rr.partial_r:.2f}\nAdjusted 95% CI {rr.partial_ci_low:.2f} to {rr.partial_ci_high:.2f}; FDR {rr.partial_fdr:.3f}',transform=ax.transAxes,fontsize=6)
save(fig,'Figure4_stability_covariation')
print('Exported four figures in SVG, PDF and PNG')
