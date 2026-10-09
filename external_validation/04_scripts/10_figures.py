from pathlib import Path
import numpy as np,pandas as pd,matplotlib as mpl,matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1]
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':8,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,'legend.frameon':False})
mods=['TAL_salt','TAL_CaMg','DCT_NaCl','DCT_Mg','DCT_Ca']
labels=['TAL salt','TAL Ca/Mg','DCT NaCl','DCT Mg','DCT Ca']
short=['TAL\nsalt','TAL\nCa/Mg','DCT\nNaCl','DCT\nMg','DCT\nCa']
blue='#356F85';orange='#C16B3A';gray='#8C8C8C'
def save(fig,name):
 for ext in ['svg','pdf']:fig.savefig(P/f'06_figures/{name}.{ext}',facecolor='white')
 fig.savefig(P/f'06_figures/{name}.png',dpi=300,facecolor='white')
 fig.savefig(P/f'06_figures/{name}.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)
# Main external comparison: both assays retain their prespecified roles.
from matplotlib.lines import Line2D
frames={a:pd.read_csv(P/f'05_results/{a}_program_comparisons.tsv',sep='\t').query("scenario == 'primary'") for a in ['sn','sc']}
score_frames={a:pd.read_csv(P/f'05_results/{a}_donor_scores.tsv',sep='\t').query("scenario == 'primary'") for a in ['sn','sc']}
all_estimates=pd.concat(frames.values(),ignore_index=True)
all_estimates.to_csv(P/'06_figures/Fig4_source_estimates.tsv',sep='\t',index=False)
lo=min(all_estimates.ci_low.min(),all_estimates.effect.min());hi=max(all_estimates.ci_high.max(),all_estimates.effect.max())
xlo=np.floor((lo-.15)*2)/2;xhi=np.ceil((hi+.15)*2)/2
fig=plt.figure(figsize=(174/25.4,156/25.4))
for row,assay in enumerate(['sn','sc']):
 top=.945-row*.46
 fig.text(.02,top,'Primary single-nucleus assay' if assay=='sn' else 'Secondary single-cell assay',fontsize=9,fontweight='bold',color='#26363D')
 for col,condition in enumerate(['AKI','CKD']):
  left=.13+col*.50;bottom=.60-row*.46
  ax=fig.add_axes([left,bottom,.255,.235]);z=frames[assay].query('contrast == @condition').set_index('module').loc[mods]
  for k,(_,v) in enumerate(z.iterrows()):
   if k%2==0:ax.axhspan(k-.44,k+.44,color='#F2F5F6',zorder=0)
   if np.isfinite(v.ci_low):
    ax.plot([v.ci_low,v.ci_high],[k,k],color=blue,lw=1.2)
    ax.plot(v.effect,k,'o',ms=4.3,mec=blue,mfc=blue if v.fdr_10<.05 else 'white',zorder=3)
   else:ax.plot(v.effect,k,marker='D',ms=4,mec=gray,mfc='white',zorder=3)
   tr=blended_transform_factory(ax.transAxes,ax.transData)
   ax.text(1.08,k,f'{int(v.n_reference)}/{int(v.n_disease)}',transform=tr,ha='center',va='center',fontsize=7.5)
   ax.text(1.35,k,f'{v.fdr_10:.3f}' if np.isfinite(v.fdr_10) else 'NA',transform=tr,ha='center',va='center',fontsize=7.5)
  ax.axvline(0,ls=':',color='#8F999D',lw=.8)
  ax.set(xlim=(xlo,xhi),ylim=(4.6,-.65),yticks=range(5),yticklabels=labels,xlabel='Mean log2 CPM difference')
  ax.tick_params(axis='y',length=0)
  ax.text(1.08,1.04,'n R/D',transform=ax.transAxes,ha='center',fontsize=7.5)
  ax.text(1.35,1.04,'FDR',transform=ax.transAxes,ha='center',fontsize=7.5)
  fig.text(left,bottom+.275,f'{condition} versus reference',fontsize=8.5)
  fig.text(left-.105,bottom+.276,'abcd'[row*2+col],fontsize=10,fontweight='bold')
fig.text(.02,.039,'Filled circles: FDR < 0.05     Open diamonds: descriptive estimates',fontsize=7.5)
save(fig,'Fig4')
# Every eligible external donor remains visible in the supplement.
fig,axes=plt.subplots(2,2,figsize=(174/25.4,166/25.4))
fig.subplots_adjust(left=.12,right=.98,top=.94,bottom=.12,wspace=.31,hspace=.44)
for row,assay in enumerate(['sn','sc']):
 for col,condition in enumerate(['AKI','CKD']):
  ax=axes[row,col];ss=score_frames[assay].query('contrast == @condition')
  for k,mod in enumerate(mods):
   for grp,off,color in [('Ref',-.14,blue),(condition,.14,orange)]:
    vals=ss[(ss.module==mod)&(ss.group==grp)].sort_values('patient').score.to_numpy()
    if not len(vals):continue
    ax.scatter(k+off+np.linspace(-.065,.065,len(vals)),vals,s=13,color=color,marker='o' if grp=='Ref' else '^',alpha=.9,edgecolors='white',linewidth=.25)
    ax.plot([k+off-.09,k+off+.09],[vals.mean()]*2,color=color,lw=1.6)
  ax.set(xticks=range(5),xticklabels=short,xlim=(-.55,4.55))
  if col==0:ax.set_ylabel('Donor mean log2 CPM score')
  ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
  ax.set_title(('Primary snRNA: ' if assay=='sn' else 'Secondary scRNA: ')+condition,loc='left',fontsize=9,pad=13)
  ax.text(-.23,1.055,'abcd'[row*2+col],transform=ax.transAxes,fontsize=10,fontweight='bold')
fig.legend([Line2D([],[],marker='o',ls='',color=blue),Line2D([],[],marker='^',ls='',color=orange)],['Reference','AKI or CKD'],loc='lower center',bbox_to_anchor=(.5,.025),ncol=2,fontsize=8)
pd.concat(score_frames.values(),ignore_index=True).to_csv(P/'06_figures/FigS5_source_scores.tsv',sep='\t',index=False)
save(fig,'FigS5')
d=pd.read_csv(P/'05_results/sn_program_comparisons.tsv',sep='\t')
scens=['primary','min10','min50','cortex_any','reference_state','exclude_COVID','sex_adjusted']
slabs=['Primary','10 nuclei','50 nuclei','Any cortex','Reference\nstate','Exclude\nCOVID','Sex adjusted']
fig=plt.figure(figsize=(183/25.4,130/25.4))
maxabs=1.5
for i,condition in enumerate(['AKI','CKD']):
 ax=fig.add_axes([.17,.58-i*.43,.69,.28])
 z=d[(d.contrast==condition)&d.scenario.isin(scens)].pivot(index='module',columns='scenario',values='effect').loc[mods,scens]
 im=ax.imshow(z.to_numpy(),vmin=-maxabs,vmax=maxabs,cmap='RdBu_r',aspect='auto')
 for row in range(5):
  for col in range(7):
   v=z.iloc[row,col]
   ax.text(col,row,f'{v:+.2f}' if np.isfinite(v) else 'NA',ha='center',va='center',color='white' if abs(v)>1.0 else '#222222',fontsize=8)
 ax.set(yticks=range(5),yticklabels=labels,xticks=range(7),xticklabels=slabs)
 ax.tick_params(length=0);ax.add_patch(Rectangle((-.49,-.49),.98,4.98,fill=False,edgecolor='#222222',lw=1.2))
 ax.set_title(f'{condition} versus reference',loc='left',fontsize=10,pad=10)
 fig.text(.035,.92-i*.43,'ab'[i],fontweight='bold',fontsize=11)
 for sp in ax.spines.values():sp.set_visible(False)
cax=fig.add_axes([.90,.20,.018,.58])
fig.colorbar(im,cax=cax,label='Mean log2 CPM difference')
d[d.scenario.isin(scens)].to_csv(P/'06_figures/FigS6_source_estimates.tsv',sep='\t',index=False)
save(fig,'FigS6')
print('Created Fig4, FigS5, FigS6 in SVG/PDF/PNG/TIFF')

