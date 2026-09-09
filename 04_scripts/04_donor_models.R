#!/usr/bin/env Rscript
args <- commandArgs(FALSE)
script <- sub('^--file=', '', args[grepl('^--file=',args)])
root <- dirname(dirname(normalizePath(script)))
if(nzchar(Sys.getenv('KIDNEY_R_LIB'))) .libPaths(c(Sys.getenv('KIDNEY_R_LIB'),.libPaths()))
library(edgeR);library(limma);library(jsonlite)
options(stringsAsFactors=FALSE)
writeLines(capture.output(sessionInfo()),file.path(root,'08_logs/R_sessionInfo.txt'))
stopifnot(file.exists(file.path(root,'03_protocol/freeze_receipt.json')))
out <- file.path(root,'05_results');dat <- file.path(root,'02_data/derived')
readcounts <- function(name){d<-read.delim(gzfile(file.path(dat,paste0('pseudobulk_',name,'.tsv.gz'))),check.names=FALSE);x<-as.matrix(d[,-1]);rownames(x)<-d[[1]];storage.mode(x)<-'double';x}
counts <- readcounts('celltype')
meta <- read.delim(file.path(dat,'pseudobulk_celltype_samples.tsv'),check.names=FALSE)
clinical <- read.delim(file.path(out,'donor_clinical_source_table.tsv'),check.names=FALSE)
meta$age_mid <- clinical$age_band_midpoint_approx[match(meta$donor,clinical$donor)]
mods <- fromJSON(file.path(root,'03_protocol/modules_v1.json'),simplifyVector=FALSE)
mods <- lapply(mods,function(x){x$genes<-unlist(x$genes,use.names=FALSE);x})
transport_genes <- unique(unlist(lapply(Filter(function(x) x$family!='context',mods),function(x)x$genes)))
stopifnot(identical(colnames(counts),meta$pb_id))
scores <- list(); comparisons <- list(); genes_main <- list(); genes_sens <- list(); qclog <- list(); norm_main <- list(); coverage <- list()

exact_p <- function(y,group){n1<-sum(group=='UUO');obs<-abs(mean(y[group=='UUO'])-mean(y[group=='Control']));cmb<-combn(length(y),n1);vals<-colMeans(matrix(y[cmb],nrow=n1));delta<-vals-(sum(y)-vals*n1)/(length(y)-n1);mean(abs(delta)>=obs-1e-12)}
run <- function(ct, scenario, nmin=20, omit=NULL, adjustment='none', custom=NULL, min_group=4){
  if(is.null(custom)){
    m<-meta[meta$celltype==ct,,drop=FALSE];x<-counts[,m$pb_id,drop=FALSE]
  }else{m<-custom$m;x<-custom$x}
  k<-m$n_nuclei>=nmin & !(m$donor %in% omit);m<-m[k,,drop=FALSE];x<-x[,k,drop=FALSE]
  m$group<-factor(m$group,levels=c('Control','UUO'));m$sex<-factor(m$sex)
  n0<-sum(m$group=='Control');n1<-sum(m$group=='UUO');testable<-min(n0,n1)>=min_group
  if(min(n0,n1)<2){qclog[[length(qclog)+1]]<<-data.frame(celltype=ct,scenario=scenario,n_control=n0,n_uuo=n1,status='INSUFFICIENT_COVERAGE',n_genes_tested=0);return(invisible(NULL))}
  form<-switch(adjustment,none=~group,sex=~group+sex,sex_age=~group+sex+age_mid)
  design<-model.matrix(form,m)
  if(qr(design)$rank<ncol(design)||nrow(design)-ncol(design)<3)testable<-FALSE
  yall<-DGEList(x);keep<-filterByExpr(yall,design=model.matrix(~group,m))
  y<-yall[keep,,keep.lib.sizes=TRUE];y<-calcNormFactors(y,method='TMM')
  yall$samples$norm.factors<-y$samples$norm.factors
  logcpm<-cpm(yall,log=TRUE,prior.count=2)
  if(scenario=='main'){
    write.table(data.frame(pb_id=m$pb_id,donor=m$donor,celltype=ct,group=m$group,n_nuclei=m$n_nuclei,total_umis=yall$samples$lib.size,tmm_factor=y$samples$norm.factors),file.path(out,paste0('normalization_',ct,'.tsv')),sep='\t',row.names=FALSE,quote=FALSE)
    norm_main[[ct]]<<-list(logcpm=logcpm,meta=m)
  }
  status<-if(testable)'TESTED' else 'DESCRIPTIVE_ONLY'
  qclog[[length(qclog)+1]]<<-data.frame(celltype=ct,scenario=scenario,n_control=n0,n_uuo=n1,status=status,n_genes_tested=if(testable)sum(keep) else 0)
  for(mm in mods){
    if(!(mm$celltype %in% c(ct,'ALL_TARGET')))next
    members<-intersect(mm$genes,rownames(logcpm));valid<-length(members)>=max(2,ceiling(length(mm$genes)*.8))
    if(!valid)next
    s<-colMeans(logcpm[members,,drop=FALSE])
    scores[[length(scores)+1]]<<-data.frame(scenario=scenario,celltype=ct,module=mm$module,label=mm$label,family=mm$family,donor=m$donor,group=m$group,sex=m$sex,n_nuclei=m$n_nuclei,score_mean_log2cpm=s)
    eff<-mean(s[m$group=='UUO'])-mean(s[m$group=='Control']);lo<-hi<-pv<-perm<-NA_real_
    if(testable){
      if(adjustment=='none'){
        tt<-t.test(s[m$group=='UUO'],s[m$group=='Control']);lo<-tt$conf.int[1];hi<-tt$conf.int[2];pv<-tt$p.value
        perm<-exact_p(s,m$group)
      }else{
        model<-lm.fit(design,s);model<-lm(s~design-1);ss<-summary(model)$coefficients;ii<-which(colnames(design)=='groupUUO')
        eff<-ss[ii,1];ci<-confint(model)[ii,];lo<-ci[1];hi<-ci[2];pv<-ss[ii,4]
      }
    }
    comparisons[[length(comparisons)+1]]<<-data.frame(scenario=scenario,celltype=ct,module=mm$module,label=mm$label,family=mm$family,n_control=n0,n_uuo=n1,n_genes=length(members),effect=eff,ci_low=lo,ci_high=hi,p_value=pv,permutation_p=perm,adjustment=adjustment,status=status)
    if(scenario=='main')for(g in members)coverage[[length(coverage)+1]]<<-data.frame(celltype=ct,module=mm$module,gene=g,n_donors=length(s),n_donors_nonzero=sum(x[g,]>0),n_donors_cpm1=sum(cpm(yall)[g,]>=1),passes_gene_filter=keep[match(g,rownames(x))])
  }
  if(testable){
    v<-voom(y,design,plot=FALSE);fit<-eBayes(lmFit(v,design),robust=TRUE)
    tt<-topTable(fit,coef='groupUUO',number=Inf,sort.by='none',confint=.95)
    tt$gene<-rownames(tt);tt$celltype<-ct;tt$scenario<-scenario;tt$n_control<-n0;tt$n_uuo<-n1
    tt$family<-if(ct %in% c('TAL','DCT'))'primary_gene' else 'exploratory_gene'
    if(scenario=='main')genes_main[[ct]]<<-tt else genes_sens[[length(genes_sens)+1]]<<-tt[tt$gene %in% transport_genes,,drop=FALSE]
  }
  cat(ct,scenario,status,n0,n1,'\n');flush.console()
}

for(ct in c('TAL','DCT','PC','IC-A','IC-B')){
  run(ct,'main')
  run(ct,'min10',nmin=10);run(ct,'min50',nmin=50)
  eligible<-meta[meta$celltype==ct&meta$n_nuclei>=20,,drop=FALSE]
  if(min(table(factor(eligible$group,levels=c('Control','UUO'))))>=4){
    for(d in eligible$donor)run(ct,paste0('LODO_',d),omit=d,min_group=3)
    run(ct,'sex_adjusted',adjustment='sex')
    run(ct,'sex_ageband_adjusted',adjustment='sex_age')
  }
}

# Prespecified composition probes, using separately aggregated original counts.
stcounts<-readcounts('state');stmeta<-read.delim(file.path(dat,'pseudobulk_state_samples.tsv'),check.names=FALSE)
custom_sum<-function(ct,states){
  sm<-stmeta[stmeta$state %in% states,,drop=FALSE];donors<-sort(unique(meta$donor))
  xx<-sapply(donors,function(d){ids<-sm$pb_id[sm$donor==d];if(length(ids))rowSums(stcounts[,ids,drop=FALSE]) else rep(0,nrow(stcounts))})
  mm<-meta[match(paste0(donors,'__',ct),meta$pb_id),,drop=FALSE]
  mm$n_nuclei<-sapply(donors,function(d)sum(sm$n_nuclei[sm$donor==d]));colnames(xx)<-mm$pb_id;rownames(xx)<-rownames(stcounts)
  list(x=xx,m=mm)
}
run('TAL','author_healthy_states',custom=custom_sum('TAL',c('cTAL1','cTAL2','mTAL')))
run('DCT','author_healthy_states',custom=custom_sum('DCT',c('DCT1','DCT2')))
run('TAL','macula_densa_excluded',custom=custom_sum('TAL',c('cTAL1','cTAL2','mTAL','TAL Injured','TAL Inflammatory')))

cmp<-do.call(rbind,comparisons);cmp$fdr<-NA_real_
for(s in unique(cmp$scenario))for(f in unique(cmp$family)){
  k<-cmp$scenario==s&cmp$family==f;cmp$fdr[k]<-p.adjust(cmp$p_value[k],method='BH')
}
write.table(cmp,file.path(out,'module_comparisons_all_scenarios.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(cmp[cmp$scenario=='main',],file.path(out,'module_comparisons_primary_run.tsv'),sep='\t',row.names=FALSE,quote=FALSE,na='NA')
write.table(do.call(rbind,scores),file.path(out,'module_scores_all_scenarios.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,coverage),file.path(out,'module_gene_coverage.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
gm<-do.call(rbind,genes_main);gm$FDR_family<-NA_real_
for(f in unique(gm$family)){k<-gm$family==f;gm$FDR_family[k]<-p.adjust(gm$P.Value[k],method='BH')}
write.table(gm,gzfile(file.path(out,'pseudobulk_gene_results.tsv.gz')),sep='\t',row.names=FALSE,quote=FALSE)
write.table(gm[gm$gene %in% transport_genes,],file.path(out,'fixed_transport_gene_results.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,genes_sens),file.path(out,'fixed_transport_gene_sensitivity.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
write.table(do.call(rbind,qclog),file.path(out,'model_eligibility_log.tsv'),sep='\t',row.names=FALSE,quote=FALSE)
saveRDS(norm_main,file.path(dat,'normalized_primary_logcpm.rds'))
for(ct in names(norm_main)){
  z<-norm_main[[ct]];write.table(data.frame(gene=rownames(z$logcpm),z$logcpm,check.names=FALSE),gzfile(file.path(dat,paste0('logcpm_',ct,'.tsv.gz'))),sep='\t',row.names=FALSE,quote=FALSE)
}
summary<-list(completed_at=as.character(Sys.time()),method='edgeR TMM + limma-voom robust eBayes; Welch module comparisons',n_gene_tests=nrow(gm),primary_transport=cmp[cmp$scenario=='main'&cmp$family=='primary_transport',],status='MODELS_COMPLETED')
write_json(summary,file.path(root,'08_logs/model_run_receipt.json'),pretty=TRUE,auto_unbox=TRUE,na='null')
print(summary$primary_transport)
