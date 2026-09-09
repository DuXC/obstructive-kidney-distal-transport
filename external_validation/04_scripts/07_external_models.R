#!/usr/bin/env Rscript
extra_lib <- Sys.getenv("KIDNEY_R_LIB", unset="")
if(nzchar(extra_lib)) .libPaths(c(extra_lib,.libPaths()))
library(edgeR);library(limma);library(jsonlite)
args<-commandArgs(trailingOnly=TRUE);P<-args[1];assay<-if(length(args)>1)args[2] else "sn"
stopifnot(file.exists(file.path(P,"03_protocol/freeze_receipt.json")))
agg<-readRDS(file.path(P,paste0("02_data/",assay,"_donor_pseudobulk.rds")))
mods<-fromJSON(file.path(P,"03_protocol/modules_v1.json"),simplifyVector=FALSE)
mods<-Filter(function(m)m$family=="primary_transport",mods)
mods<-lapply(mods,function(m){m$genes<-unlist(m$genes);m})
rows<-list();scores<-list();members<-list();genes<-list();logs<-list();norms<-list()
run<-function(ct,disease,scenario="primary",variant="primary",nmin=20,omit=NULL,exclude_covid=FALSE,adjust_sex=FALSE,min_group=4){
 a<-agg[[paste(variant,ct,sep="__")]];d<-a$metadata;x<-a$counts
 k<-d$group %in% c("Ref",disease)&d$n_nuclei>=nmin&!d$patient%in%omit
 if(exclude_covid)k<-k&d$subgroup!="COV-AKI"
 d<-d[k,,drop=FALSE];x<-x[,k,drop=FALSE]
 d$group<-factor(d$group,levels=c("Ref",disease));d$sex<-factor(d$sex)
 n0<-sum(d$group=="Ref");n1<-sum(d$group==disease);can<-min(n0,n1)>=min_group
 relevant<-Filter(function(m)m$celltype==ct,mods)
 if(min(n0,n1)>=2){
  design<-model.matrix(~group,d);ya<-DGEList(x);keep<-filterByExpr(ya,design=design)
  yy<-calcNormFactors(ya[keep,,keep.lib.sizes=TRUE],method="TMM")
  ya$samples$norm.factors<-yy$samples$norm.factors;lcpm<-cpm(ya,log=TRUE,prior.count=2)
 }else{lcpm<-NULL;can<-FALSE}
 sexmodel<-if(adjust_sex&&nlevels(d$sex)>1)model.matrix(~group+sex,d) else NULL
 if(adjust_sex&&(is.null(sexmodel)||qr(sexmodel)$rank<ncol(sexmodel)||nrow(sexmodel)-ncol(sexmodel)<3))can<-FALSE
 for(mm in relevant){
  found<-intersect(mm$genes,rownames(x));covered<-length(found)>=max(2,ceiling(length(mm$genes)*.8))
  status<-if(!covered)"INSUFFICIENT_GENE_COVERAGE" else if(!can)"INSUFFICIENT_DONOR_OR_DESIGN_COVERAGE" else "TESTED"
  effect<-low<-high<-pv<-NA_real_
  if(covered&&!is.null(lcpm)){
   s<-colMeans(lcpm[found,,drop=FALSE]);effect<-mean(s[d$group==disease])-mean(s[d$group=="Ref"])
   scores[[length(scores)+1]]<<-data.frame(assay=assay,scenario=scenario,lineage=ct,contrast=disease,module=mm$module,patient=d$patient,group=as.character(d$group),sex=d$sex,n_nuclei=d$n_nuclei,score=s,altered_fraction=d$altered_fraction)
   if(can){
    if(adjust_sex){
     fit<-lm(s~group+sex,data=d);co<-summary(fit)$coefficients[2,];ci<-confint(fit)[2,];effect<-co[1];low<-ci[1];high<-ci[2];pv<-co[4]
    }else{
     tt<-t.test(s[d$group==disease],s[d$group=="Ref"]);low<-tt$conf.int[1];high<-tt$conf.int[2];pv<-tt$p.value
    }
   }
  }
  rows[[length(rows)+1]]<<-data.frame(assay=assay,scenario=scenario,lineage=ct,contrast=disease,module=mm$module,label=mm$label,n_reference=n0,n_disease=n1,n_genes=length(found),n_fixed_genes=length(mm$genes),effect=effect,ci_low=low,ci_high=high,p_value=pv,status=status)
  if(scenario=="primary")for(g in mm$genes){
   represented<-g%in%rownames(x)
   members[[length(members)+1]]<<-data.frame(assay=assay,contrast=disease,lineage=ct,module=mm$module,gene=g,represented=represented,n_donors=ncol(x),n_detected=if(represented)sum(x[g,]>0) else NA,n_gene_filter=if(represented&&!is.null(lcpm))as.integer(keep[match(g,rownames(x))]) else NA)
  }
 }
 if(scenario=="primary"&&!is.null(lcpm)){
  norms[[paste(ct,disease,sep="__")]]<<-list(logcpm=lcpm,metadata=d)
  logs[[length(logs)+1]]<<-data.frame(assay=assay,contrast=disease,lineage=ct,patient=d$patient,n_nuclei=d$n_nuclei,total_counts=ya$samples$lib.size,TMM_factor=ya$samples$norm.factors)
  if(can){
   v<-voom(yy,design,plot=FALSE);fit<-eBayes(lmFit(v,design),robust=TRUE)
   tt<-topTable(fit,coef=2,number=Inf,sort.by="none",confint=.95)
   genes[[length(genes)+1]]<<-data.frame(gene=rownames(tt),assay=assay,lineage=ct,contrast=disease,n_reference=n0,n_disease=n1,tt,row.names=NULL)
  }
 }
}
for(disease in c("AKI","CKD"))for(ct in c("TAL","DCT")){
 run(ct,disease)
 for(th in c(10,50))run(ct,disease,scenario=paste0("min",th),nmin=th)
 run(ct,disease,scenario="reference_state",variant="reference_state")
 run(ct,disease,scenario="sex_adjusted",adjust_sex=TRUE)
 if(assay=="sn"){
  run(ct,disease,scenario="cortex_any",variant="cortex_any")
  run(ct,disease,scenario="exclude_COVID",exclude_covid=TRUE)
 }
 m<-agg[[paste("primary",ct,sep="__")]]$metadata
 ds<-m$patient[m$group%in%c("Ref",disease)&m$n_nuclei>=20]
 for(d in ds)run(ct,disease,scenario=paste0("LODO_",d),omit=d,min_group=3)
}
res<-do.call(rbind,rows);res$fdr_10<-NA_real_
for(s in unique(res$scenario)){
 k<-res$scenario==s
 if(!startsWith(s,"LODO_"))res$fdr_10[k]<-p.adjust(res$p_value[k],method="BH",n=10)
}
write.table(res,file.path(P,paste0("05_results/",assay,"_program_comparisons.tsv")),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
write.table(do.call(rbind,scores),file.path(P,paste0("05_results/",assay,"_donor_scores.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,members),file.path(P,paste0("05_results/",assay,"_member_coverage.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
write.table(do.call(rbind,logs),file.path(P,paste0("05_results/",assay,"_normalization.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
g<-do.call(rbind,genes);g$FDR_joint<-p.adjust(g$P.Value,method="BH")
write.table(g,gzfile(file.path(P,paste0("05_results/",assay,"_gene_results.tsv.gz"))),sep="\t",quote=FALSE,row.names=FALSE)
fixed<-unique(unlist(lapply(mods,function(m)m$genes)))
write.table(g[g$gene%in%fixed,],file.path(P,paste0("05_results/",assay,"_fixed_gene_results.tsv")),sep="\t",quote=FALSE,row.names=FALSE)
saveRDS(norms,file.path(P,paste0("02_data/",assay,"_normalized_program_inputs.rds")))
writeLines(capture.output(sessionInfo()),file.path(P,paste0("08_qa/",assay,"_R_sessionInfo.txt")))
receipt<-list(status="COMPLETE",assay=assay,n_gene_tests=nrow(g),primary_comparisons=res[res$scenario=="primary",])
write_json(receipt,file.path(P,paste0("08_qa/",assay,"_model_receipt.json")),pretty=TRUE,auto_unbox=TRUE,na="null")
print(res[res$scenario=="primary",c("contrast","module","n_reference","n_disease","effect","ci_low","ci_high","fdr_10","status")],row.names=FALSE)

