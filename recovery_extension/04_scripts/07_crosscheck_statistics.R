#!/usr/bin/env Rscript
phase <- normalizePath('.')
meta <- read.delim('02_audit/sample_batches.tsv',check.names=FALSE)
meta <- meta[meta$condition!='Normal',];rownames(meta)<-meta$sample
meta$scan_day <- sub(' .*','',trimws(meta$scan_date))
errors <- list(); tests<-0L
record <- function(label, actual, expected) {
  if(any(!is.finite(actual)) || max(abs(actual-expected))>1e-8) stop(label,': numerical disagreement ',max(abs(actual-expected)))
}
for(kind in c('program','gene')) {
 values<-read.delim(paste0('05_results/',if(kind=='program') 'program_scores.tsv' else 'gene_values.tsv'))
 key<-if(kind=='program') 'module' else 'human_symbol'
 for(analysis in c('primary','detectable')) {
  ref<-read.delim(paste0('05_results/',kind,'_',analysis,'_contrasts.tsv'))
  for(i in which(ref$status=='tested')) {
   r<-ref[i,];v<-values[values$analysis==analysis & values[[key]]==r$entity,]; y<-setNames(v$value,v$sample)
   getv<-function(cond,day) na.omit(y[meta$sample[meta$condition==cond & meta$duration_days==r$duration_days & meta$recovery_days==day]])
   if(r$kind=='same_time') {
    a<-getv('RUUO',r$recovery_days);b<-getv('Sham',r$recovery_days)
    t<-t.test(a,b,var.equal=FALSE);est<-unname(t$estimate[1]-t$estimate[2]);ci<-t$conf.int;p<-t$p.value
   } else if(r$kind=='late_vs_acute') {
    a<-getv('RUUO',28);b<-getv('RUUO',0);t<-t.test(a,b,var.equal=FALSE)
    est<-unname(t$estimate[1]-t$estimate[2]);ci<-t$conf.int;p<-t$p.value
   } else {
    groups<-list(getv('RUUO',28),getv('Sham',28),getv('RUUO',0),getv('Sham',0))
    signs<-c(1,-1,-1,1);est<-sum(signs*vapply(groups,mean,0.))
    vv<-vapply(groups,function(x)var(x)/length(x),0.);nn<-vapply(groups,length,1L)
    se<-sqrt(sum(vv));df<-sum(vv)^2/sum(vv^2/(nn-1));ci<-est+c(-1,1)*qt(.975,df)*se;p<-2*pt(-abs(est/se),df)
   }
   record(paste(kind,analysis,r$contrast,r$entity),c(est,ci,p),c(r$effect,r$ci_low,r$ci_high,r$p));tests<-tests+1L
  }
  pfull<-ifelse(is.na(ref$p),1,ref$p);adjusted<-p.adjust(pfull,'BH')
  record(paste(kind,analysis,'BH'),adjusted[ref$status=='tested'],ref$fdr[ref$status=='tested'])
 }
 ref<-read.delim(paste0('05_results/',kind,'_date_adjusted.tsv'))
 for(i in which(ref$status=='tested')) {
  r<-ref[i,];v<-values[values$analysis=='primary' & values[[key]]==r$entity,];y<-setNames(v$value,v$sample)
  dat<-meta[meta$duration_days==r$duration_days & meta$recovery_days==r$recovery_days,];dat$value<-y[dat$sample];dat<-dat[is.finite(dat$value),]
  shared<-names(which(tapply(dat$condition,dat$scan_day,function(x)length(unique(x)))==2));dat<-dat[dat$scan_day %in% shared,]
  dat$obstruction<-as.integer(dat$condition=='RUUO')
  fit<-lm(value ~ obstruction + factor(scan_day),dat)
  X<-model.matrix(fit);inv<-solve(crossprod(X));res<-residuals(fit);h<-hatvalues(fit)
  meat<-crossprod(X,X*as.vector((res/(1-h))^2));cov<-inv%*%meat%*%inv
  est<-coef(fit)['obstruction'];se<-sqrt(cov['obstruction','obstruction']);df<-df.residual(fit)
  ci<-est+c(-1,1)*qt(.975,df)*se;p<-2*pt(-abs(est/se),df)
  record(paste(kind,'date',r$entity,r$contrast),c(est,se,ci,p),c(r$effect,r$se,r$ci_low,r$ci_high,r$p));tests<-tests+1L
 }
}
writeLines(sprintf('{"status":"PASS_R_STATISTICAL_CROSSCHECK","tested_contrasts":%d,"absolute_tolerance":1e-8,"checks":"Welch t.test, four-group contrasts, BH families, lm HC3 covariance"}',tests),'09_qa/statistical_crosscheck.json')
cat('PASS independent R numerical crosscheck:',tests,'tested comparisons\n')
