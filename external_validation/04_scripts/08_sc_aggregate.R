#!/usr/bin/env Rscript
extra_lib <- Sys.getenv("KIDNEY_R_LIB", unset="")
if(nzchar(extra_lib)) .libPaths(c(extra_lib,.libPaths()))
library(Matrix);library(jsonlite)
P<-commandArgs(trailingOnly=TRUE)[1]
f<-file.path(P,"02_data/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Counts_03282022.RDS.gz")
message("Reading official scCv3 count matrix")
x<-readRDS(gzcon(gzfile(f,"rb")));message("Class ",paste(class(x),collapse=",")," dimensions ",paste(dim(x),collapse=" x "))
stopifnot(inherits(x,"sparseMatrix"),length(dim(x))==2,!is.null(colnames(x)),!is.null(rownames(x)),!anyDuplicated(colnames(x)),!anyDuplicated(rownames(x)))
m<-read.delim(gzfile(file.path(P,"01_sources/GSE183276_Kidney_Healthy-Injury_Cell_Atlas_scCv3_Metadata_03282022.txt.gz")),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
stopifnot(!anyDuplicated(rownames(m)),setequal(rownames(m),colnames(x)))
m<-m[colnames(x),];m$patient<-trimws(as.character(m$patient))
stopifnot(all(vapply(split(m$condition.l1,m$patient),function(z)length(unique(z))==1,logical(1))))
sn<-read.delim(gzfile(file.path(P,"01_sources/GSE183277_Kidney_Healthy-Injury_Cell_Atlas_snCv3_Metadata_03282022.txt.gz")),row.names=1,check.names=FALSE,stringsAsFactors=FALSE)
shared<-intersect(unique(m$patient),unique(trimws(as.character(sn$patient))))
nz<-length(x@x);bad<-0;fractional<-0;total<-0;mn<-Inf;mx<-0
for(i in seq(1,nz,by=2000000)){
 v<-x@x[i:min(nz,i+1999999)];bad<-bad+sum(!is.finite(v)|v<0);fractional<-fractional+sum(abs(v-round(v))>1e-8);total<-total+sum(v);mn<-min(mn,min(v));mx<-max(mx,max(v))
}
stopifnot(bad==0)
metadata_umi_equal<-isTRUE(all.equal(as.numeric(Matrix::colSums(x)),as.numeric(m$nCount_RNA),tolerance=1e-8))
write_json(list(n_cells=ncol(x),nonzero_entries=nz,fractional_count_values=fractional,invalid_count_values=bad,nonzero_min=mn,nonzero_max=mx,library_sum_quantiles=quantile(Matrix::colSums(x)),metadata_umi_equal=metadata_umi_equal,umi_correlation=cor(Matrix::colSums(x),m$nCount_RNA)),file.path(P,"08_qa/sc_representation_audit.json"),pretty=TRUE,auto_unbox=TRUE)
keep<-!m$patient%in%shared
masks<-list(primary=keep,reference_state=keep&m$state.l2=="reference")
aggregates<-list();coverage<-list();audit<-list()
for(variant in names(masks))for(ct in c("TAL","DCT")){
 idx<-which(masks[[variant]]&m$subclass.l1==ct);donors<-sort(unique(m$patient[idx]));j<-match(m$patient[idx],donors)
 design<-sparseMatrix(i=seq_along(idx),j=j,x=1,dims=c(length(idx),length(donors)))
 counts<-as.matrix(x[,idx,drop=FALSE]%*%design);colnames(counts)<-donors
 dat<-do.call(rbind,lapply(donors,function(d){
  ids<-idx[m$patient[idx]==d]
  stopifnot(length(unique(m$condition.l1[ids]))==1,length(unique(m$sex[ids]))==1)
  data.frame(patient=d,group=unique(m$condition.l1[ids]),subgroup=paste(unique(m$condition.l2[ids]),collapse=";"),
   sex=unique(m$sex[ids]),tissue_type=paste(sort(unique(m$tissue_type[ids])),collapse=";"),
   n_nuclei=length(ids),n_libraries=length(unique(m$library[ids])),mean_cortex=NA_real_,
   altered_fraction=mean(m$state.l2[ids]!="reference"),stringsAsFactors=FALSE)
 }))
 stopifnot(identical(dat$patient,colnames(counts)))
 original_sum<-sum(x[,idx,drop=FALSE]);stopifnot(isTRUE(all.equal(sum(counts),original_sum,tolerance=1e-12)))
 key<-paste(variant,ct,sep="__");aggregates[[key]]<-list(counts=counts,metadata=dat);coverage[[key]]<-cbind(variant=variant,lineage=ct,dat)
 audit[[key]]<-list(cells=length(idx),donors=length(donors),aggregated_counts=sum(counts),source_selected_counts=original_sum,absolute_sum_difference=abs(sum(counts)-original_sum))
 message(key,": ",length(idx)," cells, ",length(donors)," reported patient IDs")
}
saveRDS(aggregates,file.path(P,"02_data/sc_donor_pseudobulk.rds"))
write.table(do.call(rbind,coverage),file.path(P,"05_results/sc_donor_coverage.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
write.table(data.frame(gene=rownames(x)),file.path(P,"05_results/sc_measured_genes.tsv"),sep="\t",row.names=FALSE,quote=FALSE)
qa<-list(status=if(fractional>0)"REQUIRES_COUNT_REPRESENTATION_REVIEW" else "PASS",n_genes=nrow(x),n_cells=ncol(x),n_donors=length(unique(m$patient)),nonzero_entries=nz,total_counts=total,
 invalid_count_values=bad,fractional_count_values=fractional,nonzero_min=mn,nonzero_max=mx,library_sum_quantiles=quantile(Matrix::colSums(x)),metadata_umi_equal=metadata_umi_equal,umi_correlation=cor(Matrix::colSums(x),m$nCount_RNA),unique_barcode_linkage=TRUE,
 excluded_shared_patient_ids=shared,n_excluded_shared_patient_ids=length(shared),source_file=basename(f),aggregate_audits=audit,
 donor_interpretation="Source metadata patient IDs; PRE019 technical libraries merged. Shared reported patient IDs excluded; same atlas, not an independent study.")
write_json(qa,file.path(P,"08_qa/sc_aggregation_qa.json"),pretty=TRUE,auto_unbox=TRUE)
message("Completed donor aggregation")
