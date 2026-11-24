if (.Platform$OS.type == "windows") try(Sys.setlocale("LC_CTYPE","English_United States.utf8"), silent=TRUE)
suppressPackageStartupMessages({library(boot);library(jsonlite)})
script <- normalizePath(sub("^--file=","",commandArgs(FALSE)[grepl("^--file=",commandArgs(FALSE))]),winslash="/")
root <- dirname(script)
args <- commandArgs(trailingOnly=TRUE)
outdir <- if (length(args)) args[1] else file.path(dirname(dirname(dirname(root))),"review-output","china")
dir.create(outdir,recursive=TRUE,showWarnings=FALSE)
d <- read.csv(file.path(root,"species.csv"),check.names=FALSE)
d <- d[order(d$id_no),]
stopifnot(nrow(d)==1820L,!anyDuplicated(d$id_no))
ref <- read.csv(file.path(root,"KruskalWallis_Tests.csv"),check.names=FALSE)
metrics <- c("Q_mean_period_mean","Q_5pct_period_mean","R_low_period_mean","ANCI_2000","ANCI_2020","ANCI_change")
labels <- c("Mean habitat quality, 2000-2020","5th-percentile habitat quality, 2000-2020","Proportion of low-quality habitat, 2000-2020","ANCI in 2000","ANCI in 2020","ANCI change, 2020 minus 2000")
fast_H <- function(x,g,k) {
 n <- length(x)
 r <- rank(x,ties.method="average")
 ties <- rle(sort(x))$lengths
 correction <- 1-sum(ties^3-ties)/(n^3-n)
 stopifnot(correction>0)
 counts <- tabulate(g,nbins=k)
 sums <- rowsum(r,g,reorder=TRUE)[,1]
 (12/(n*(n+1))*sum(sums^2/counts)-3*(n+1))/correction
}
bootstrap_eta <- function(dat,indices,columns,k) {
 z <- dat[indices,]
 g <- as.integer(z$group)
 vapply(columns,function(m) max(0,(fast_H(z[[m]],g,k)-k+1)/(nrow(z)-k)),numeric(1))
}
B <- 10000L
out <- list();diagnostics<-list()
for(j in seq_along(c("category","class"))) {
 grouping<-c("category","class")[j]
 levels_g<-if(grouping=="category")c("CR","EN","VU","NT","LC") else c("Mammals","Amphibians","Reptiles","Plants")
 dat<-d
 dat$group<-factor(dat[[grouping]],levels=levels_g)
 k<-nlevels(dat$group);n<-nrow(dat)
 stopifnot(!anyNA(dat$group))
 estimates<-numeric(length(metrics))
 for(m in seq_along(metrics)){
  kw<-kruskal.test(dat[[metrics[m]]],dat$group)
  h<-unname(kw$statistic)
  old<-ref[ref$metric==labels[m]&ref$grouping==grouping,]
  stopifnot(nrow(old)==1L,n==old$n,k==old$df+1,abs(h-old$H)<1e-8,abs(kw$p.value/old$p_raw-1)<1e-8)
  stopifnot(abs(fast_H(dat[[metrics[m]]],as.integer(dat$group),k)-h)<1e-9)
  estimates[m]<-(h-k+1)/(n-k)
  stopifnot(abs(estimates[m]-old$rank_eta_squared)<1e-10)
 }
 seed<-20261108L+j
 set.seed(seed)
 # Ordinary stratified bootstrap; each row is a complete species record.
 obj<-boot(dat,bootstrap_eta,R=B,strata=dat$group,columns=metrics,k=k)
 stopifnot(all(is.finite(obj$t)),all(obj$t>=0),all(obj$t<=1),max(abs(obj$t0-estimates))<1e-10)
 saveRDS(obj$t,file.path(outdir,paste0("Table_S9_bootstrap_",grouping,".rds")))
 for(m in seq_along(metrics)){
  ci<-boot.ci(obj,conf=.95,type="perc",index=m)$percent[1,4:5]
  old<-ref[ref$metric==labels[m]&ref$grouping==grouping,]
  out[[length(out)+1L]]<-data.frame(metric=labels[m],grouping=grouping,n=n,groups=k,H=old$H,df=old$df,p_raw=old$p_raw,rank_eta_squared=estimates[m],ci_low=ci[1],ci_high=ci[2],bootstrap_replicates=B,seed=seed)
 }
 diagnostics[[grouping]]<-list(seed=seed,counts=as.list(table(dat$group)),bootstrap_replicates=B,nonfinite=0,estimate_match=TRUE)
 cat("Completed ",grouping,": ",B," resamples for 6 metrics.\n",sep="");flush.console()
}
result<-do.call(rbind,out)
key<-paste(ref$metric,ref$grouping)
result<-result[match(key,paste(result$metric,result$grouping)),]
stopifnot(nrow(result)==12L,all(result$ci_low<=result$ci_high))
write.csv(result,file.path(outdir,"Table_S9_effect_sizes_CI.csv"),row.names=FALSE)
write_json(list(R=R.version.string,boot=as.character(packageVersion("boot")),method="Two-sided stratified percentile bootstrap",diagnostics=diagnostics,original_statistics_reproduced=TRUE,formula="max(0, (H-k+1)/(n-k))"),file.path(outdir,"bootstrap_integrity.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
writeLines(capture.output(sessionInfo()),file.path(outdir,"R-sessionInfo.txt"))
print(result[,c("metric","grouping","rank_eta_squared","ci_low","ci_high")],row.names=FALSE)

