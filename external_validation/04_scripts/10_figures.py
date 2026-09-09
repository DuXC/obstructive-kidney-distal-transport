from pathlib import Path
import numpy as np,pandas as pd,matplotlib as mpl,matplotlib.pyplot as plt
from matplotlib.transforms import blended_transform_factory
from matplotlib.patches import Rectangle
P=Path(__file__).resolve().parents[1]
mpl.rcParams.update({'font.family':'sans-serif','font.sans-serif':['Arial','DejaVu Sans'],'font.size':7.5,'svg.fonttype':'none','pdf.fonttype':42,'axes.spines.right':False,'axes.spines.top':False,'axes.linewidth':.7,'legend.frameon':False})
mods=['TAL_salt','TAL_CaMg','DCT_NaCl','DCT_Mg','DCT_Ca']
labels=['TAL salt','TAL Ca/Mg','DCT NaCl','DCT Mg','DCT Ca']
short=['TAL\nsalt','TAL\nCa/Mg','DCT\nNaCl','DCT\nMg','DCT\nCa']
blue='#466F84';orange='#AC502A';gray='#8C8C8C'
def save(fig,name):
 for ext in ['svg','pdf']:fig.savefig(P/f'06_figures/{name}.{ext}',facecolor='white')
 fig.savefig(P/f'06_figures/{name}.png',dpi=300,facecolor='white')
 fig.savefig(P/f'06_figures/{name}.tiff',dpi=600,facecolor='white',pil_kwargs={'compression':'tiff_lzw'})
 plt.close(fig)
for assay,name in [('sn','Fig4'),('sc','FigS5')]:
 d=pd.read_csv(P/f'05_results/{assay}_program_comparisons.tsv',sep='\t');s=pd.read_csv(P/f'05_results/{assay}_donor_scores.tsv',sep='\t')
 d=d[d.scenario=='primary'].copy();s=s[s.scenario=='primary'].copy()
 d.to_csv(P/f'06_figures/{name}_source_estimates.tsv',sep='\t',index=False)
 s.to_csv(P/f'06_figures/{name}_source_scores.tsv',sep='\t',index=False)
 fig=plt.figure(figsize=(183/25.4,175/25.4))
 all_min=min(d.ci_low.min(),d.effect.min());all_max=max(d.ci_high.max(),.2)
 xlo=np.floor((all_min-.2)*2)/2;xhi=np.ceil((all_max+.2)*2)/2
 for col,condition in enumerate(['AKI','CKD']):
  left=.125+.5*col
  ax=fig.add_axes([left,.65,.255,.25])
  z=d[d.contrast==condition].set_index('module').loc[mods]
  for k,(_,row) in enumerate(z.iterrows()):
   y=4-k
   if k%2==0:ax.axhspan(y-.45,y+.45,color='#F1F4F6',zorder=0)
   if np.isfinite(row.ci_low):
    ax.plot([row.ci_low,row.ci_high],[y,y],color=blue,lw=1.3,zorder=2)
    ax.plot(row.effect,y,'o',ms=4.5,mec=blue,mfc=blue if row.fdr_10<.05 else 'white',zorder=3)
   else:
    ax.plot(row.effect,y,marker='D',ms=4,mec=gray,mfc='white',zorder=3)
   tr=blended_transform_factory(ax.transAxes,ax.transData)
   ax.text(1.06,y,f'{int(row.n_reference)}/{int(row.n_disease)}',transform=tr,ha='center',va='center',fontsize=7)
   ax.text(1.36,y,f'{row.fdr_10:.3f}' if np.isfinite(row.fdr_10) else 'NA',transform=tr,ha='center',va='center',fontsize=7)
  ax.axvline(0,ls=':',color='#999999',lw=.8)
  ax.set(xlim=(xlo,xhi),ylim=(-.6,4.8),yticks=range(5),yticklabels=labels[::-1],xlabel='Mean log2 CPM difference')
  ax.tick_params(axis='y',length=0)
  ax.text(1.06,1.04,'n R/D',transform=ax.transAxes,ha='center',fontsize=7)
  ax.text(1.36,1.04,'FDR',transform=ax.transAxes,ha='center',fontsize=7)
  fig.text(left,.949,f'{condition} versus reference',fontsize=9)
  fig.text(left-.1,.95,'ab'[col],fontsize=10,fontweight='bold')
  ax=fig.add_axes([.105+.5*col,.16,.35,.32])
  zz=s[s.contrast==condition]
  for k,mod in enumerate(mods):
   for grp,off,color in [('Ref',-.14,blue),(condition,.14,orange)]:
    vals=zz[(zz.module==mod)&(zz.group==grp)].sort_values('patient').score.to_numpy()
    if not len(vals):continue
    jitter=np.linspace(-.065,.065,len(vals))
    ax.scatter(k+off+jitter,vals,s=12,color=color,marker='o' if grp=='Ref' else '^',edgecolors='none',alpha=.85)
    ax.plot([k+off-.09,k+off+.09],[vals.mean()]*2,color=color,lw=1.5)
  ax.set(xticks=range(5),xticklabels=short,xlim=(-.55,4.55))
  if col==0:ax.set_ylabel('Donor mean log2 CPM score')
  ax.grid(axis='y',alpha=.14,zorder=0);ax.set_axisbelow(True)
  fig.text(.105+.5*col,.525,f'{condition}: individual donor scores',fontsize=9)
  fig.text(.025+.5*col,.525,'cd'[col],fontsize=10,fontweight='bold')
 from matplotlib.lines import Line2D
 handles=[Line2D([],[],marker='o',ls='',color=blue,label='Reference'),Line2D([],[],marker='^',ls='',color=orange,label='AKI or CKD')]
 fig.legend(handles=handles,loc='lower center',bbox_to_anchor=(.5,.045),ncol=2,fontsize=8)
 save(fig,name)
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

