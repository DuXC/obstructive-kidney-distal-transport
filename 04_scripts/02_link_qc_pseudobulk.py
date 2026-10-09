#!/usr/bin/env python3
"""Validate sparse integer RNA matrices and aggregate author-retained nuclei by donor."""
from pathlib import Path
import hashlib,json,time,resource
import numpy as np
import pandas as pd
import h5py
from scipy import sparse

ROOT=Path(__file__).resolve().parents[1]
RAW=ROOT/'02_data/raw'; OUT=ROOT/'02_data/derived'; RES=ROOT/'05_results'
meta=pd.read_csv(RAW/'GSE254185_cell_meta_data.tsv.gz',sep='\t').rename(columns={'Unnamed: 0':'barcode','Sample':'donor','Annotation.Lvl1':'celltype','Annotation.Lvl2':'state','Sex':'sex'})
assert len(meta)==46957 and meta.barcode.is_unique and meta.donor.nunique()==12
meta['group']=np.where(meta.donor.str.startswith('UUO'),'UUO','Control')
assert (meta.barcode.str.extract(r'^l(\d)_')[0].values==meta.Library.str.extract(r'Library(\d)')[0].values).all()
keys={}
for level in ['celltype','state']:
    keys[level]=pd.DataFrame(sorted(set(zip(meta.donor,meta[level]))),columns=['donor',level])
    keys[level]['pb_id']=keys[level].donor+'__'+keys[level][level]
    keys[level]=keys[level].merge(meta.groupby(['donor',level]).size().rename('n_nuclei').reset_index(),on=['donor',level])
    keys[level]=keys[level].merge(meta[['donor','group','sex']].drop_duplicates(),on='donor')
markers='SLC12A1 UMOD KCNJ1 CLCNKA CLCNKB SLC12A3 PVALB TRPM6 CALB1 TRPV5 SLC8A1 AQP2 AQP3 AQP4 AVPR2 SCNN1A SCNN1B SCNN1G ATP6V1B1 ATP6V0A4 ATP6V1G3 SLC4A1 SLC26A4 FOXI1 RHCG RHBG CA2 CLDN16 CLDN19 KCNJ10 KCNJ16 SLC34A1 LRP2 CUBN PTPRC LST1 CD3D MS4A1 PECAM1 VWF COL1A1 COL3A1 PROM1 DCDC2 SPP1 ITGB6 ITGB8 CCL2 CXCL1 ICAM1 VCAM1 HAVCR1'.split()
receipt=[];qc=[];markrows=[];pbulk={};pb_lib=[];ref_genes=None;direct_umis=0
start=time.monotonic()
for lib in range(1,7):
    path=RAW/f'GSE254185_Library{lib}_gex_matrix.h5'
    with h5py.File(path,'r') as f:
        assert list(f)==['unknown'],list(f)
        g=f['unknown'];genes=g['gene_names'][:].astype(str);ids=g['genes'][:].astype(str)
        assert np.array_equal(genes,ids),'Unexpected IDs/symbols: inspect before aggregation'
        assert len(set(genes))==len(genes),'Duplicate gene symbols require documented resolution'
        if ref_genes is None:
            ref_genes=genes
            for level in keys:pbulk[level]=np.zeros((len(genes),len(keys[level])),dtype=np.int64)
        assert np.array_equal(genes,ref_genes),'Library feature mismatch'
        barcodes=g['barcodes'][:].astype(str)
        assert len(set(barcodes))==len(barcodes)
        d=g['data'][:];ind=g['indices'][:];ptr=g['indptr'][:];shape=tuple(g['shape'][:])
        assert np.isfinite(d).all() and (d>=0).all() and (d==np.floor(d)).all()
        assert ptr[0]==0 and ptr[-1]==len(d) and (np.diff(ptr)>=0).all()
        assert ind.min()>=0 and ind.max()<shape[0] and shape==(len(genes),len(barcodes))
        x=sparse.csc_matrix((d.astype(np.int64),ind,ptr),shape=shape);del d,ind,ptr
    m=meta[meta.Library.eq(f'Library{lib}')].copy()
    lookup=pd.Index('l'+str(lib)+'_'+barcodes)
    idx=lookup.get_indexer(m.barcode)
    assert (idx>=0).all(),m.loc[idx<0,'barcode'].tolist()[:10]
    assert len(set(idx))==len(m)
    y=x[:,idx];del x
    umi=np.asarray(y.sum(axis=0)).ravel();ng=np.diff(y.indptr);mt=np.asarray(y[np.char.startswith(genes,'MT-'),:].sum(axis=0)).ravel()/np.maximum(umi,1)*100
    assert (umi>0).all()
    m['n_counts']=umi;m['n_genes']=ng;m['percent_mito']=mt
    qc.append(m);direct_umis+=int(umi.sum())
    receipt.append(dict(library=f'Library{lib}',h5_group='unknown',shape=list(map(int,shape)),stored_dtype='float64',all_values_nonnegative_integers=True,retained_nuclei=len(m),h5_barcodes_not_in_author_metadata=len(barcodes)-len(m),missing_metadata_barcodes=0,retained_umis=int(umi.sum()),gene_symbols_unique=True))
    for level,k in keys.items():
        ids=pd.Index(k.pb_id).get_indexer(m.donor+'__'+m[level]);assert (ids>=0).all()
        one=sparse.csr_matrix((np.ones(len(m),dtype=np.int64),(np.arange(len(m)),ids)),shape=(len(m),len(k)))
        add=np.asarray((y@one).toarray());pbulk[level]+=add
        if level=='celltype':
            cols=np.flatnonzero(add.sum(axis=0));pb_lib.append((lib,cols,add[:,cols]))
    present=[a for a in markers if a in set(genes)]
    z=y[pd.Index(genes).get_indexer(present),:].toarray()
    norm=np.log1p(z/umi[None,:]*1e4)
    for level in ['celltype','state']:
        for (donor,ct),mi in m.groupby(['donor',level]).indices.items():
            for j,gene in enumerate(present):
                markrows.append(dict(donor=donor,group=m.iloc[mi[0]].group,level=level,annotation=ct,gene=gene,n_nuclei=len(mi),fraction_detected=float((z[j,mi]>0).mean()),mean_log1p_cp10k=float(norm[j,mi].mean())))
    del y,z,norm
    print(json.dumps(receipt[-1]),flush=True)

q=pd.concat(qc,ignore_index=True);assert q.barcode.is_unique and len(q)==len(meta)
q.to_csv(OUT/'nucleus_qc.tsv.gz',sep='\t',index=False)
pd.DataFrame({'gene':ref_genes}).to_csv(OUT/'genes.tsv',sep='\t',index=False)
for level,counts in pbulk.items():
    assert int(counts.sum())==direct_umis
    np.savez_compressed(OUT/f'pseudobulk_{level}.npz',counts=counts,genes=ref_genes,samples=keys[level].pb_id.values.astype(str))
    pd.DataFrame(counts,index=ref_genes,columns=keys[level].pb_id).to_csv(OUT/f'pseudobulk_{level}.tsv.gz',sep='\t',index_label='gene')
    k=keys[level];k['total_umis']=counts.sum(axis=0);k['detected_genes']=(counts>0).sum(axis=0)
    k.to_csv(OUT/f'pseudobulk_{level}_samples.tsv',sep='\t',index=False)
    if level=='celltype':k.to_csv(RES/'donor_celltype_qc.tsv',sep='\t',index=False)
for lib,cols,c in pb_lib:
    np.savez_compressed(OUT/f'pseudobulk_library{lib}.npz',counts=c,genes=ref_genes,samples=keys['celltype'].pb_id.values[cols].astype(str))
pd.DataFrame(markrows).to_csv(RES/'marker_audit_by_donor.tsv.gz',sep='\t',index=False)
q.groupby(['donor','Library','group','sex','celltype']).size().rename('n_nuclei').reset_index().to_csv(RES/'donor_library_celltype_counts.tsv',sep='\t',index=False)
donor=q.groupby('donor').agg(group=('group','first'),sex=('sex','first'),n_nuclei=('barcode','size'),n_libraries=('Library','nunique'),median_umis=('n_counts','median'),median_genes=('n_genes','median'),median_percent_mito=('percent_mito','median')).reset_index()
donor['tumour_context']=np.where(donor.group.eq('UUO'),'upper urinary tract urothelial carcinoma','unobstructed renal tumour; individual histology unlinked')
donor.to_csv(RES/'donor_qc.tsv',sep='\t',index=False)
composition=q.groupby(['donor','group','celltype']).size().rename('n_nuclei').reset_index()
epithelia=['ATL','CNT','DCT','DTL','IC-A','IC-B','PC','PEC','PT','Podocyte','TAL']
total=q.groupby('donor').size();epi=q[q.celltype.isin(epithelia)].groupby('donor').size()
composition['fraction_all_captured']=composition.n_nuclei/composition.donor.map(total)
composition['fraction_captured_epithelia']=np.where(composition.celltype.isin(epithelia),composition.n_nuclei/composition.donor.map(epi),np.nan)
composition.to_csv(RES/'captured_composition.tsv',sep='\t',index=False)
state=q.groupby(['donor','group','celltype','state']).size().rename('n_nuclei').reset_index()
den=q.groupby(['donor','celltype']).size()
state['fraction_within_celltype']=[n/den.loc[(d,c)] for d,c,n in zip(state.donor,state.celltype,state.n_nuclei)]
state.to_csv(RES/'author_state_composition.tsv',sep='\t',index=False)
summary=dict(status='PASS',n_nuclei=len(q),n_donors=q.donor.nunique(),n_genes=len(ref_genes),n_libraries=6,retained_total_umis=direct_umis,all_metadata_mapped_once=True,raw_count_sum_equals_both_pseudobulk_sums=True,libraries=receipt,wall_seconds=time.monotonic()-start,peak_rss_bytes=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss)
(ROOT/'08_logs/linkage_qc_receipt.json').write_text(json.dumps(summary,indent=2))
print(json.dumps({k:v for k,v in summary.items() if k!='libraries'},indent=2))
