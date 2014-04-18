#install.packages("vegan")
#install.packages("HMP")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("gtools")
#install.packages("BiodiversityR")
#install.packages("LiblineaR")
#install.packages("BatchJobs")
#install.packages("Rcmdr")
##base as.Date() method is brittle and does very little error checking
#install.packages("date")
##handle timestamps
#install.packages("timeDate")
#install.packages("pheatmap")

##two alternative implementations of the same stability
##selection paper, and different classification methods
#GLMs
#install.packages("c060")
## used by c060, install here just in case
#install.packages("glmnet")
#combined L1 and L2 penalties
#this will need options(warn=0) on Windows
#install.packages("quadrupen")
##for boxcoxfit
#install.packages("geoR")
## parallel backend that works on Windows with foreach()
#install.packages("doSNOW")
#?gridExtra
#install.packages("elasticnet")
#install.packages("BioMark")
#install.packages("fdrtool")
#install.packages("tikzDevice") #for knitr
#install.packages("markdown") #for knit2html
#source("http://bioconductor.org/biocLite.R")
#biocLite("multtest")
#biocLite("GeneSelector")
#biocLite("RColorBrewer")
## To upgrade Bioconductor packages, use: source("http://bioconductor.org/biocLite.R"); biocLite("BiocUpgrade")
## Normal R-Studio upgrade apparently does not upgrade Bioconductor.

set_trace_options<-function() {
  #setup more verbose error reporting
  #set warn to 2 to convert warinings to errors
  options(warn = 0, keep.source = TRUE, error = 
            quote({ 
              cat("Environment:\n", file=stderr()); 
              
              # TODO: setup option for dumping to a file (?)
              # Set `to.file` argument to write this to a file for post-mortem debugging    
              dump.frames();  # writes to last.dump
              
              #
              # Debugging in R
              #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/index.shtml
              #
              # Post-mortem debugging
              #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/pmd.shtml
              #
              # Relation functions:
              #   dump.frames
              #   recover
              # >>limitedLabels  (formatting of the dump with source/line numbers)
              #   sys.frame (and associated)
              #   traceback
              #   geterrmessage
              #
              # Output based on the debugger function definition.
              
              n <- length(last.dump)
              calls <- names(last.dump)
              cat(paste("  ", 1L:n, ": ", calls, sep = ""), sep = "\n", file=stderr())
              cat("\n", file=stderr())
              
              if (!interactive()) {
                q()
              }
            }))
}



## bind variable to the global environment
## can be used in debugging
make.global <- function(var) {
  assign(deparse(substitute(var)),var,envir=globalenv()) 
}

packages = c(        
  #for quantcut; mask as little as possible (masks permute, 
  "gtools", #pos = "package:base", 
  "reshape2",
  "plyr",
  # for str_split with max splits
  "stringr", 
  "ggplot2",
  "vegan", 
  "BiodiversityR", 
  "LiblineaR", 
  ## glmnet is exported to snow cluster by c060, but it
  ## forgets to load it first, so we do it here
  "glmnet", 
  "c060", 
  "geoR", 
  "foreach", 
  "iterators", 
  "doSNOW", 
  "fdrtool", 
  #"elasticnet", 
  #"BioMark", 
  "HMP", 
  "BatchJobs", 
  "boot",
  "vegan",
  #Bioconductor
  #"multtest",
  "GeneSelector",
  #"knitr",
  "lattice",
  "date",
  "timeDate",
  #for llist
  "Hmisc"
)

for (package in packages) {
  suppressMessages(library(package,character.only=T))
}

## replace oldnames with newnames in vector allnames
## to be used when allnames = names(data.frame)
replace.col.names<-function(allnames,oldnames,newnames) {
  allnames[match(oldnames,allnames)] = newnames
  return(allnames)
}

## Return a copy of data frames with selected columns removed
drop_columns <- function(df,drop_names=c()) {
  mask_col_ignore = names(df) %in% drop_names
  return ( df[,!mask_col_ignore] )
}

count_matr_from_df<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  as.matrix(dat[!mask_col_ignore])
}

split_count_df<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  if(!all(mask_col_ignore)) {
    m = as.matrix(dat[!mask_col_ignore])
  }
  else {
    m = NULL
  }
  attr = dat[mask_col_ignore]
  list(count=m,attr=attr)
}

row_normalize_matr<-function(x) {
  x/rowSums(x)
}


row_normalize<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  matnorm<-row_normalize_matr(x$count)
  datnorm<-as.data.frame(matnorm) 
  res = cbind(x$attr,datnorm)
  #print(apply(res[!(names(res) %in% col_ignore)],1,sum) == 1)
  stopifnot(all(abs(apply(res[!(names(res) %in% col_ignore)],1,sum) - 1)<1e-6))
  res
}

all_normalize<-function(dat,col_ignore=c(),norm.func=NULL) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  if( is.null(norm.func) ) {
    norm.func = ihs
  }
  matnorm<-norm.func(x$count)
  datnorm<-as.data.frame(matnorm) 
  cbind(x$attr,datnorm)
}

## If the other_cnt column is already present, it will be incremented with counts of clades
## relegated to the "other" in this call; otherwise, the new column with this name will be
## created.
## Count columns will be sorted in decreasing order of the column mean frequencies, so that
## you can easily subset the count matrix later to only keep N most abundant columns.
count_filter<-function(dat,col_ignore=c(),min_max_frac=0.0,min_max=10,min_median=0,min_row_sum=100,other_cnt="other") {
  x<-split_count_df(dat,col_ignore)
  row_cnt = rowSums(x$count)
  row_sel = row_cnt >= min_row_sum
  cnt = x$count[row_sel,]
  row_cnt = row_cnt[row_sel]
  attr = x$attr[row_sel,]
  cnt_norm = row_normalize_matr(cnt)
  col_sum = colSums(cnt_norm)
  cnt = cnt[,order(col_sum,decreasing=T)]
  cnt_col_sel = cnt[,apply(cnt_norm,2,max) >= min_max_frac]
  cnt_col_sel = cnt_col_sel[,apply(cnt_col_sel,2,max) >= min_max]
  cnt_col_sel = cnt_col_sel[,apply(cnt_col_sel,2,median) >= min_median]  
  cnt_col_other = as.data.frame(row_cnt - rowSums(cnt_col_sel))
  if (!all(cnt_col_other[,1]==0) && !is.null(other_cnt)) {
    names(cnt_col_other) = c(other_cnt)
    if (other_cnt %in% colnames(cnt_col_sel)) {
      cnt_col_sel[,other_cnt] = cnt_col_sel[,other_cnt] + cnt_col_other[[other_cnt]]
    }
    else {
      return (cbind(attr,as.data.frame(cnt_col_sel),cnt_col_other))
    }
  }
  return (cbind(attr,as.data.frame(cnt_col_sel)))
}

sample.rows<-function(x,size,replace=F,prob=NULL) {
  x[sample(dim(x)[1],size,replace,prob),]
}

list_to_df<-function(x,col_names=NULL,row_names=NULL) {
  y = as.data.frame(do.call(rbind, x),row.names=row_names)
  if (!is.null(col_names)) {
    names(y) = col_names
  }
  y
}

tryCatchAndWarn<-function(expr,catch.warnings=F,catch.errors=T) {
  ok=T
  #somehow if both warning and error arguments are used
  #in a single call to tryCatch, then error handler is not
  #called - probably only the handler for whatever happens first
  #is called
  
  if(catch.warnings) {
    warnHand = function(w) {
      warning(paste("Warning caught in tryCatch; returning NULL: ",w,"\n"))
      ok<<-F
    }
  }
  
  if(catch.errors) {
    errHand = function(e) {
      warning(paste("Error caught in tryCatch; returning NULL: ",e,"\n"))
      ok<<-F
    }
  }
  
  if(catch.warnings && catch.errors) {
    res = tryCatch(tryCatch(expr,warning=warnHand),error=errHand)
  }
  else if(catch.warnings && !catch.errors) {
    res = tryCatch(expr,warning=warnHand)
  }
  else if(!catch.warnings && catch.errors) {
    res = tryCatch(expr,error=errHand)
  }
  else {
    res = tryCatch(expr)
  }
  
  #options(warn=warn_saved)
  if (!ok) {
    #print("Error occured in test, returning NA")
    res = NULL
  }
  return(res)
}

sorted.freq<-function(count.df,col_ignore=c(),col_group="group") {
  freq = row_normalize(count.df,col_ignore=col_ignore)
  freq_m = count_matr_from_df(freq,col_ignore=col_ignore)
  group = count.df[,col_group]
  freq.mean = aggregate(freq_m,list(group),mean,simplify=TRUE)
  row.names(freq.mean) = freq.mean[,"Group.1"]
  freq.mean$Group.1 = NULL
  freq.mean = t(as.matrix(freq.mean))
  return(freq.mean)
}


kelvin.rnames.to.ids<-function(row.names) {
  col_names = c("id_repl","group")
  ids = list_to_df(str_split(row.names,"_",2),row_names=row.names(data),col_names=col_names)
}

read.kelvin.summary.matr<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="Sample")
  stopifnot(all(data$Total==rowSums(subset(data,select=c(-Total)))))
  data$Total = NULL
  ids=kelvin.rnames.to.ids(row.names(data))
  data$id_repl = ids$id_repl
  data$group = ids$group
  data
}

power.kelvin<-function() {
  group_sel = c("Stool","Saliva")
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls")
  
  data = data_all[data_all$group %in% group_sel,]
  
  ##x = Data.filter(data$x, "sample", 1000, 10)
  data = count_filter(data,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)
  cnt_m = count_matr_from_df(data,col_ignore=attr_names)
  
  
  
  x1 = cnt_m[data$group == group_sel[1],]
  x2 = cnt_m[data$group == group_sel[2],]
  
  n_samp = 15
  
  n_iter = 10
  
  alpha = 0.05
  
  p_vals = c()
  for (i_iter in seq(n_iter)) {
    mygroup = list(sample.rows(x1,n_samp),sample.rows(x2,n_samp))
    p_vals = c(p_vals,Xmcupo.sevsample(mygroup,dim(x1)[2])[[2]])
  }
  power = mean(p_vals<=alpha,na.rm = T)
  
  print(power)
  list(x1,x2)
}


read.koren<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="Taxon")
  data = as.data.frame(t(data))
  rnames = row.names(data)
  data$time = as.ordered(as.integer(substr(rnames,3,3)))
  data$id_repl = as.factor(as.integer(substr(rnames,4,300)))
  data
}

wilcox.test.multi <- function(data,resp.vars,group.var,subset) {
  pvals = foreach(resp.var=resp.vars,.combine=c) %do% {
    form = as.formula(paste(resp.var,group.var,sep="~"))
    print(form)
    pval = wilcox.test(form,
                       data=data,
                       subset=subset)$p.value
    print(pval)
    pval
  }
  #p.adjust(pvals,method="BH")
}

wilcox.test.multi.fast <- function(tr.mat,group) {
  RankingWilcoxon(tr.mat,group,type="unpaired",pvalues=T)@pval
}

#test.method argument is for passing function definition to the new
#parallel process on WIndows under Snow. Set it to wilcox.test.multi.fast
#when calling 'boot'
booted.wilcox.test.multi.fast <- function(n,mat,group,indices,test.method) {
  library(GeneSelector)
  ind.n = sample(indices, n, replace = FALSE, prob = NULL)
  d <- t(mat)[,ind.n] # allows boot to select n samples
  g <- group[ind.n]
  #print(summary(g))
  test.method(d,g)
}


# run one test from boot indices
# pass exact=F because boot convertes warnings "cannot compute exact p-value with ties"
#into errors
booted.wilcox.test <- function(n,is.paired,exact,data,indices) {
  d <- data[indices[1:n],] # allows boot to select n sample
  with(d,wilcox.test(freq.x,freq.y,paired=is.paired,exact=exact))$p.value
}


booted.adonis.test.koren <- function(n,data,indices) {
  #nlev = length(levels(data$time))
  d <- data[indices,]
  if(!is.null(n)) {
    d <- sample.rows(d,n*length(levels(d$time)),replace=F) #allows boot to select n sample
  }
  l = count_matr_from_df(d,col_ignore=data_factors)
  ad.res = adonis(l~time,data=d,permutations=1000,method="bray")
  #print (ad.res)
  ad.res$aov.tab$"Pr(>F)"[1]
}

power.koren<-function() {
  
  data = read.koren("GG_OTU_table_L6.txt")
  data_factors = c("time","id_repl")
  #data = count_filter(data,col_ignore=data_factors,min_max_frac=0.1,min_row_sum=0)
  data_col = melt(data,id.vars=data_factors,variable.name="clade",value.name="freq")
  data_time = merge(data_col[data_col$time==1,],data_col[data_col$time==2,],by=c("id_repl","clade"))
  data_time_summary = ddply(data_time,c("clade"),summarize,
                            mean.x=mean(freq.x), sd.x=sd(freq.x), n.x = length(freq.x),
                            mean.y=mean(freq.y), sd.y=sd(freq.y), n.y = length(freq.y),
                            mean_paired=mean(freq.x-freq.y),
                            sd_paired=sd(freq.x-freq.y),n_paired=length(freq.x))
  data_time_summary$effect_t_paired = data_time_summary$mean_paired/data_time_summary$sd_paired
  #get pooled sd (original group ratios)
  data_time_summary$sd = with(data_time_summary,sqrt(((n.x-1)*sd.x**2+(n.y-1)*sd.y**2)/(n.x+n.y)))
  data_time_summary$effect_t = with(data_time_summary,(mean.x-mean.y)/sd)
  
  
  
  #clade_sel = "Root;k__Bacteria;p__Actinobacteria"
  alpha = 0.05
  power.cut = 0.8
  n = 50 # equal sample sizes per group
  eff.div = 1
  is.paired = T
  
  mv.test.res = booted.adonis.test.koren(NULL,data,seq(nrow(data)))
  
  print (paste("adonis p-value: ",mv.test.res))
  
  power.mv.res = mean(boot(data=data,statistic=booted.adonis.test.koren,R=100,n=n,strata=data$time)$t<=alpha)
  
  print(paste("adonis power: ",power.mv.res))
  
  n.tests = length(levels(data_time$clade))
  alpha.bn = alpha/n.tests
  clade = c()
  p.value = c()
  power = c()
  for (clade_sel in levels(data_time$clade)) {
    
    data_sel = data_time[data_time$clade==clade_sel,]
    
    #test.res = with(data_sel,t.test(freq.x,freq.y,paired=is.paired))
    test.res = with(data_sel,wilcox.test(freq.x,freq.y,paired=is.paired,exact=F))
    
    clade = c(clade,clade_sel)
    p.value = c(p.value,test.res$p.value)
    
    data_summ_sel = data_time_summary[data_time_summary$clade==clade_sel,]
    if(is.paired) {
      test.type = "paired"
    }
    else {
      test.type = "two.sample"
    }
    
    sd_power = with(data_summ_sel,sqrt((n-1)/n*(sd.x**2+sd.y**2)/2))
    #power.res = with(data_summ_sel,power.t.test(n=n,delta=(mean.x-mean.y)/eff.div,sd=sd_power,sig.level=alpha.bn,type=test.type))$power
    #power.res = with(data_summ_sel,mean(replicate(1000, wilcox.test(rnorm(n,mean.x,sd_power), rnorm(n,mean.y,sd_power), paired=is.paired,exact=F)$p.value)<=alpha.bn))
    power.res = mean(boot(data=data_sel,statistic=booted.wilcox.test,R=100,n=n,is.paired=is.paired,exact=F)$t<=alpha.bn)
    
    power = c(power,power.res)
  }
  
  p.value.bh = p.adjust(p.value,method="BH")
  p.value.bn = p.value*n.tests
  test.res = data.frame(clade,p.value,p.value.bh,p.value.bn,power)
  
  #we use BN correction to pick those clades which we consider to be significant and
  #for which we measure the power
  test.res.sig = test.res[test.res$p.value.bn<=alpha,]
  test.res.sig$power.ok = test.res.sig$power >= power.cut
  
  test.res.sig = merge(test.res.sig,data_time_summary,by=c("clade"))
  
  #print (test.res.sig)
  
  print(paste("Mean power of the test: ",mean(test.res.sig$power)))
  
}

Xmcupo.effectsize.par <-
  function(par.groups,reads.groups){
    par.groups <- par.groups
    assert.non.zero.dm.par(par.groups)    
    for(i.group in seq(length(par.groups))) {
      par.groups[[i.group]]$reads = reads.groups[[i.group]]
    }
    
    N.group <- length(par.groups)
    Kc <- length(par.groups[[1]]$pi)
    N.total <- do.call(sum, lapply(par.groups, function(x){sum(x$reads)}))
    
    group.parameter.estimated <- lapply(par.groups, function(x){c(x$pi, x$theta, x$reads)})
    Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K=Kc)
    
    if(Kc >= N.group){
      pi.groups <- diag(rep(1, N.group))
    }else{
      stop("The number of taxa must be greater than the number of groups.")
    }
    
    gp.max <- lapply(as.list(1:N.group), function(x,gp.1, pi.groups){
      gpp <- gp.1[[x]]
      ret <- c(pi.groups[x,], gpp$theta, gpp$reads)
      return(ret)
    }, gp.1=par.groups, pi.groups=pi.groups)
    
    Xmcupo.gp <- Xmcupo.statistics(gp.max, K=N.group)
    CramerV<- sqrt(Xmcupo[1]/(N.total*min(Kc-1, N.group-1)))
    Mod.CramerV <- sqrt(Xmcupo[1]/(Xmcupo.gp[1]*min(Kc-1, N.group-1)))
    
    result <- c(Xmcupo[1],CramerV, Mod.CramerV)
    names(result) <- c("Chi-Squared", "Cramer Phi", "Modified-Cramer Phi")
    return(result)
  }

Xmcupo.sevsample.par <-
  function(par.groups, reads.groups, K){
    if(missing(par.groups) || missing(reads.groups) || missing(K))
      stop("par.groups, reads.groups or K is missing.")
    assert.non.zero.dm.par(par.groups)
    n.groups <- length(par.groups)
    index <- as.matrix(seq(1:n.groups))
    group.parameter.estimated <- list()
    
    for(x in index){
      par <- par.groups[[x]]
      nreads.data <- as.matrix(reads.groups[[x]])
      
      group.parameter.estimated[[x]] <- c(par$pi, par$theta, t(nreads.data))
    }
    
    Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K)
    p.value <- 1-pchisq(q=Xmcupo, df=(n.groups-1)*(K-1), ncp=0, lower.tail=TRUE)
    
    sevRAD.mean.test.upo <- list(Xmcupo, p.value)
    names(sevRAD.mean.test.upo) <- c("Xmcupo statistics", "p value")
    
    return(sevRAD.mean.test.upo)
  }


read.nistal<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="id")
  data$group = as.factor(data$group)  
  data
}

power.nistal<-function() {
  
  data_all = read.nistal("children.txt")
  groups = as.factor(c("Celiac","Control"))
  attr_names = c("group")
  all_clades = colnames(data_all)[!(colnames(data_all) %in% attr_names)]
  test_clades = as.factor(c("Prevotella.sp.","Streptococcus.sp."))
  #test_clades = as.factor(c("Comamonas.sp."))
  
  n.samp.grp = c(180,60)
  n.seq.range = seq(1500,2000,by=1)
  #n.seq.range <- rep(50,100)
  #n.samp.grp = c(8,5)
  
  dm.par.orig = list()
  
  #set number of sequences here so that we can estimate
  #effect size and DM test statistics outside of sample generation loop
  n.seq = list()
  
  for (i.group in seq(length(groups))) {
    group = groups[i.group]
    data = data_all[data_all$group == group,]
    cnt_m = count_matr_from_df(data,col_ignore=attr_names)
    dm.par.orig[[i.group]] <- DM.MoM(cnt_m)
    n.seq[[i.group]] <- sample(n.seq.range, size=n.samp.grp[i.group]) 
  }
  
  res.power = power.dirmult.range(
    dm.par.orig=dm.par.orig,
    n.seq=n.seq,
    groups=groups,
    all_clades=all_clades,
    test_clades=test_clades,
    n.samp.grp=n.samp.grp
  )    
  write.csv(res.power,"res.power.csv")
}

dirmult.comp.gamma <- function(pi,theta) {
  # standard formula for gamma
  pi * (1-theta)/theta
}

assert.non.zero.dm.par <- function(dm.par) {
  all.pi = laply(dm.par,function(l){l$pi})
  stopifnot(all(colSums(all.pi)>0))
}

batch.jobs.file.dir <- function(id) {
  tempfile(pattern = paste(id,".",sep=""), tmpdir = getwd(), fileext = ".batch_reg")
}

batch.jobs.chunked.ids <- function(reg,n.chunks) {
  chunk(getJobIds(reg), n.chunks=n.chunks, shuffle=TRUE)
}

power.dirmult.range<-function(
  dm.par.orig,
  n.seq,
  groups,
  all_clades,
  test_clades,
  n.samp.grp,
  alpha=0.05,
  effect.size.range = c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0),
  n.rep = 100,
  n.adonis.perm = 400
  
) {
  
  assert.non.zero.dm.par(dm.par.orig)
  
  n_clades = length(all_clades)
  
  test.ad.res.power.eta = c()
  r2.ad.eta = c()
  test.wil.res.power.sel.eta = c()
  mod.cramer.phi = c()
  test.Xmcupo.res = c()
  
  reg = makeRegistry("onesamp",
                     file.dir=batch.jobs.file.dir("onesamp"),
                     packages=packages,
                     skip=F)
  
  par.list = foreach (eta = effect.size.range) %do% {
    # copying the DM parameters of the alternative hypothesis from one of the parent datasets.
    dm.par <- dm.par.orig
    # interpolate the proportion vectors of the two parent datasets
    dm.par[[2]]$pi <- dm.par[[1]]$pi * (1-eta) + dm.par.orig[[2]]$pi  *(eta)
    dm.par[[2]]$theta <- dm.par[[1]]$theta * (1-eta) + dm.par.orig[[2]]$theta  *(eta)  
    # standard formula for gamma
    dm.par[[2]]$gamma <- dirmult.comp.gamma(dm.par[[2]]$pi,dm.par[[2]]$theta)
    
    mod.cramer.phi = Xmcupo.effectsize.par(dm.par,n.seq)["Modified-Cramer Phi"]
    test.Xmcupo.res = Xmcupo.sevsample.par(dm.par,n.seq,n_clades)$"p value"
    
    list(dm.par=dm.par,eta=eta,mod.cramer.phi=mod.cramer.phi,test.Xmcupo.res=test.Xmcupo.res)
  }    
  one.sample.action = function(par,i.rep,all_clades,groups,n.seq,n.samp.grp,n.adonis.perm) {
    dm.par = par$dm.par
    ### Generate a random vector of number of reads per sample
    dm.counts = matrix(nrow=0,ncol=length(all_clades),dimnames=list(NULL,all_clades))
    dm.group = c()
    for (i.group in seq(length(groups))) {
      dm.counts <- rbind(dm.counts,Dirichlet.multinomial(n.seq[[i.group]], dm.par[[i.group]]$gamma))
      dm.group <- c(dm.group,rep(groups[i.group],n.samp.grp[i.group]))
    }
    dm.freq = row_normalize_matr(dm.counts)
    dm.group = as.factor(dm.group)
    ad.res = adonis(dm.freq~dm.group,permutations=n.adonis.perm,method="bray")
    #print (ad.res)
    
    test.wil.res = foreach (clade = all_clades,.combine=c) %do% {
      #should use groups[1] and groups[2], but rep(factor,...) loses labels
      dm.freq.1 = dm.freq[dm.group==1,clade]
      dm.freq.2 = dm.freq[dm.group==2,clade]
      wilcox.test(dm.freq.1,dm.freq.2,paired=FALSE,exact=F)$p.value
    }
    
    return(list(ad.res=ad.res,test.wil.res=test.wil.res,
                eta=par$eta,
                mod.cramer.phi=par$mod.cramer.phi,
                test.Xmcupo.res=par$test.Xmcupo.res
    ))
  }
  
  batchExpandGrid(reg,one.sample.action,
                  par=par.list,i.rep=seq(n.rep),
                  more.args=list(all_clades,groups,n.seq,n.samp.grp,n.adonis.perm))
  
  submitJobs(reg, ids=batch.jobs.chunked.ids(reg,n.chunks=10))
  
  #reduceResultsDataFrame loses column names, so we use ldply on a list:
  res.all = ldply(reduceResultsList(reg,fun=function(job,res,all_clades,test_clades,alpha) {
    test.ad.res = res$ad.res$aov.tab$"Pr(>F)"[1]
    r2.ad = res$ad.res$aov.tab$R2[1]
    test.wil.res = res$test.wil.res
    names(test.wil.res) = all_clades
    test.wil.res = t(test.wil.res)
    test.ad.res.power = mean(test.ad.res<=alpha)
    for (i.col in seq(ncol(test.wil.res))) {
      test.wil.res[,i.col] = p.adjust(test.wil.res[,i.col],method="BH")
    }
    test.wil.res.power = rowMeans(test.wil.res<=alpha)
    test.wil.res.power.sel = test.wil.res.power[all_clades %in% test_clades]
    cbind(eta=res$eta,mod.cramer.phi=res$mod.cramer.phi,r2.ad,test.Xmcupo.res=res$test.Xmcupo.res,test.ad.res,test.wil.res)
  },
  test_clades=test_clades,
  all_clades=all_clades,
  alpha=alpha)
  )
  make.global(res.all)
  res = ddply(res.all,c("eta"),
              function(x,alpha) { 
                first = 
                  summarize(x[,-which(names(x) %in% all_clades)],
                            mod.cramer.phi = mean(mod.cramer.phi),
                            test.Xmcupo.res = mean(test.Xmcupo.res),
                            r2.ad = mean(r2.ad),
                            test.ad.res.power = mean(test.ad.res<=alpha))
                #use dfrm[col_ind] notation because dfrm[,col_ind] will
                #return a vector if length(col_ind) == 1
                second = 
                  as.data.frame(colMeans(x[which(names(x) %in% test_clades)]<=alpha))
                cbind(first,t(second))
              },
              alpha = alpha
  )
  res = res[order(res$eta),]
  col.names.clades = laply(names(res)[names(res) %in% test_clades],function(x){paste(x,"Mann-Whitney U power")})
  colnames(res) = c(c("Parameter scaling coeff.", "Modified Cramer Phi", 
                      "Gen. Wald-type p-value", "Adonis R2", "Adonis power"),
                    col.names.clades)
  
  return(t(res))
}


sorted.freq.kelvin<-function() {
  group_sel = c("Stool","Saliva")
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls")
  
  data = data_all[data_all$group %in% group_sel,]
  data = count_filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
  
  freq.mean = sorted.freq(data,col_ignore=attr_names)
  
  write.csv(as.matrix(sort(freq.mean[,"Saliva"],decreasing=TRUE)),"freq.top.saliva.csv")
  
  return(freq.mean)
}

sorted.freq.nistal<-function() {
  
  attr_names = c("group")
  data_all = read.nistal("children.txt")
  data = data_all
  #data = count_filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
  #return(data_all)
  
  freq.mean = sorted.freq(data,col_ignore=attr_names)
  
  write.csv(as.matrix(sort(freq.mean[,"Control"],decreasing=TRUE)),"freq.top.control.csv")
  
  return(freq.mean)
}

#freq.mean = sorted.freq.kelvin()
#freq.mean = sorted.freq.nistal()
#res = power.nistal()
#x = res[[1]]
#y = res[[2]]
#z = res[[3]]

read.humann.summary<-function(file_name) {
  #skip first 4 rows with diversity indices
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  data = data[5:nrow(data),]
  data$ID = NULL
  row.names(data) = data$NAME
  data$NAME = NULL
  data = t(data)
  row.names(data) = unlist(lapply(strsplit(row.names(data), "\\."),function(x) x[1]))
  return (as.data.frame(data))
}

annot.levels <- function(annot.type) {
  if( annot.type == "humann" ) {
    return (c("level.3"))
  }
  levs = switch(annot.type,
                cog=c("level.2"),
                c("level.2","level.3")
  )
  levs = c(levs,c("function.","gene"))
  return (levs)
}

read.mgrast.summary<-function(file_name,file_name.id.map=NULL) {
  #if quote!=NULL, resulting dataset is cropped in mid-file with a warning about
  #EOF inside a quoted string. Quotes probably do not match in pairs.
  data = read.delim(file_name, header=T,stringsAsFactors=T,quote=NULL)
  #there are more unique id values than function values, which means
  #a mapping id -> function is 1->many. We create field "gene" that gives
  #a descriptive name to the id
  data$gene = as.factor(paste(data$function.,data$id,sep="."))
  if (!is.null(file_name.id.map)) {
    id.map = read.delim(file_name.id.map, header=T,stringsAsFactors=T,row.names="mgrast_metagenome_id")
  }
  return (merge(data,id.map,by.x="metagenome",by.y="row.names"))
}


read.mothur.taxa.summary <- function(file_name) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  #where this X == NA column comes from, I have no idea. I do not see
  #it in the text file. Maybe R bug?
  data$X = NULL
  data$clade = paste(data$taxon,data$rankID,sep="_")
  return (data)
}


multi.mothur.to.abund.df <- function(data,level) {
  data.level = data[data$taxlevel==level,]
  attr = c("taxlevel","rankID","taxon","daughterlevels","total","clade")
  x = split_count_df(data.level,col_ignore=attr)
  row.names(x$count) = x$attr$clade
  return (as.data.frame(t(x$count)))
}

mgrast.to.abund.df <- function(data,level) {
  x = data[,c("project_sample_id",level,"abundance")]
  x = ddply(x,c("project_sample_id",level),summarize,abundance=sum(abundance))
  x = dcast(x,list("project_sample_id",level),value.var="abundance",fill=0.)
  row.names(x) = x$project_sample_id
  x$project_sample_id = NULL
  return (x)
}


#Merge counts data frame with meta data frame.
#Merge in this method is always done on row.names. If you need to
#merge on a component of row.names(counts), first
#get row.names(counts), split them into components, merge with 
#meta data frame, and then use the resulting frame as a new meta data
merge.counts.with.meta <- function(x,y,suffixes=c("","meta")) {
  mrg = merge(x,y,by="row.names",suffixes=suffixes)
  #Row.names column is generated by merge() when by="row.names"
  #the assignment below also serves as assertion that count records were
  #not duplicated as a result of merge
  row.names(mrg) = mrg$Row.names
  mrg$Row.names = NULL
  all.names = names(mrg)
  attr.names = all.names[!(all.names %in% names(x))]
  return (list(data=mrg,attr.names=attr.names))
}

aggregate_by_meta_data <- function(meta_data,
                                   group_col,
                                   col_ignore=c(),
                                   count_aggr=sum,
                                   attr_aggr=NULL,
                                   group_col_result_name="SampleID") {
  x = split_count_df(meta_data,col_ignore=col_ignore)
  #This is my way of assigning list attribute name from a variable
  groups = list(x$attr[,group_col])
  names(groups) = group_col_result_name
  #We need to drop resulting grouping name from the input,
  #otherwise we will get duplicated names (e.g. two SampleID columns. 
  #The assumption here is that if such name is supplied for the grouping
  #output column, then the existing column with the same name is not needed.
  if(group_col_result_name %in% names(x$attr)) {
    x$attr = drop_columns(x$attr,c(group_col_result_name))
  }
  ##TODO: make the default attr_aggr to drop all columns
  ## with more that one value in any group. 
  if(is.null(attr_aggr)) {
    attr_aggr = function(y) {y[1]}
  }
  x$attr = aggregate(x$attr,groups,attr_aggr)
  row.names(x$attr) = x$attr[,group_col_result_name]
  if(!is.null(x$count)) {
    count_names = names(x$count)
    x$count = aggregate(x$count,groups,count_aggr)  
    row.names(x$count) = x$count[,group_col_result_name]
    x$count = drop_columns(x$count,c(group_col_result_name))
    return (merge.counts.with.meta(x$count,x$attr))
  }
  else {
    return (x$attr)
  }
  
}

load.meta.choc <- function(file_name,counts.row.names) {
  meta = read.delim(file_name,header=T,stringsAsFactors=T)
  
  allnames = replace.col.names(names(meta),
                               c("Subject.ID..blinded.","SampleID","Subject.s.Gender"),
                               c("SubjectID","SampleID","Gender"))
  
  names(meta) = allnames
  
  ## Format of SubjectID from Raja:
  ## The third letter in the letter set describes whether the sample is 
  ## from the patient (P) or sibling (S). The first set of numbers is 
  ## the date of consent without year. The second set of two numbers 
  ## links patient and sibling together. For example, ARP-0328-01 and 
  ## CRS-0329-01 are family set 01
  
  Subject.ID.split = aaply(as.character(meta$SubjectID),c(1),function(x){strsplit(x,"-")[[1]]})
  
  meta$FamilyID = as.factor(Subject.ID.split[,3])
  #meta$SubjectIndFamily = as.factor(Subject.ID.split[,3])
  
  meta$Sample.type = gsub(" ",".",meta$Sample.type)
  
  #Therapy.Status
  meta$Sample.type.1 = meta$Sample.type
  
  meta$Sample.type = as.factor(unlist(apply(meta,
                                            1,
                                            function(row) {switch(row["Sample.type.1"],
                                                                  sibling ="sibling",
                                                                  patient.after.chemo="patient",
                                                                  patient.before.chemo="patient")})))
  meta$TherapyStatus = as.factor(unlist(apply(meta,
                                              1,
                                              function(row) {switch(row["Sample.type.1"],
                                                                    sibling ="before.chemo",
                                                                    patient.after.chemo="after.chemo",
                                                                    patient.before.chemo="before.chemo")})))
  
  meta$Antibiotic = as.factor(meta$Antibiotic.treatment.within.the.last.1.months. != "No")
  
  meta$SampleID.1 = as.factor(paste(meta$SubjectID,meta$Sample.type.1,sep="."))
  meta$Sample.type = as.factor(meta$Sample.type)
  meta$Sample.type.1 = as.factor(meta$Sample.type.1)
  meta$TherapyStatus = as.factor(meta$TherapyStatus)
  meta$SampleID = as.factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  make.global(meta)  
  return (meta)
}


melt.abund.meta <- function(data,id.vars,attr.names,value.name="abundance") {
  clade.names = get.clade.names(data,attr.names)
  return (melt(data,id.vars=id.vars,measure.vars=clade.names,variable.name="clade",value.name=value.name))
}

order.factor.by.total <- function(factor.val,sort.val) {
  
  o = aggregate(sort.val, list(factor.val), sum)
  o = o[order(o$x, decreasing=TRUE),1]  
  
  return (factor(factor.val, levels = o, ordered = T))
}

order.levels <- function(lev,keys) {
  lev[order(keys,decreasing=T)]
}

plot.abund.meta <- function(data,id.vars,attr.names,value.name="abundance",
                            file_name=NULL,ggp.comp=NULL,facet_grid.margins=FALSE,
                            clades.order=NULL) {
  dat = melt.abund.meta(data,id.vars=id.vars,attr.names=attr.names,value.name=value.name)
  
  if (is.null(clades.order)) {
    dat$clade = order.factor.by.total(dat$clade,dat[[value.name]])
  }
  else {
    dat$clade = factor(dat$clade,levels=clades.order,ordered=T)
  }
  
  #show only 20 top
  dat = dat[dat$clade %in% levels(dat$clade)[1:20],]
  
  if (length(id.vars) == 1) {
    wr = facet_wrap(as.formula(paste("~",id.vars[1],sep="")))
  }
  else {
    wr = facet_grid(as.formula(paste(id.vars[2],id.vars[1],sep="~")),drop=T,margins=facet_grid.margins)
  }
  
  
  clade.names = get.clade.names(data,attr.names)
  #this will be used to label each facet with number of cases in it
  facet.cnt <- ddply(.data=data, id.vars, function(x,clade.names) { c(.n=nrow(x),
                                                                      .y=mean(colMeans(as.matrix(x[,clade.names])),
                                                                              names=F)) },clade.names)
  facet.cnt$.n = paste("n =", facet.cnt$.n)
  #facet.cnt$y = facet.cnt$V2
  
  gp = ggplot(dat, aes(x=clade,y=ihs(abundance),
                       fill = clade,color = clade))+
    #geom_violin(scale= "count", trim=TRUE, adjust=1)+
    #stat_summary(fun.y="mean", geom="point",size=8,color="black",shape=4)+
    stat_summary(fun.y="mean", geom="bar")+
    #stat_summary(fun.data = mean_cl_boot, geom = "pointrange",color="black")+
    coord_flip()+
    #geom_boxplot(color="black")+
    #geom_boxplot(fill=NA,na.value=NA,outlier.size = 0)+
    #geom_point(position = "jitter")+
    #stat_identity(geom="bar")+
    #geom_bar(stat="identity")
    #scale_fill_brewer(type = "seq", palette = 1)
    wr+
    #labels facet with number of cases
    geom_text(data=facet.cnt, aes(x=1.5, y=.y, label=.n, size=1), colour="black", inherit.aes=FALSE, parse=FALSE)+
    theme(legend.position = "none",
          axis.title=element_blank(),
          axis.text.y=element_text(color=c("black","black")))
  
  if (!is.null(ggp.comp)) {
    for (g.c in ggp.comp) {
      gp = gp + g.c
    }
  }
  if (!is.null(file_name)) {
    ggsave(file_name)
  }
  return (gp)
}


std.plots <- function(abund.meta.data,
                      abund.meta.attr.names,
                      id.vars.list,
                      extra.attr=NULL,
                      label="",
                      do.diversity=F,
                      file.name.root="std.plot",
                      file.name.sfx="png",
                      res.tests=NULL) {
  
  abund.meta.data.norm = row_normalize(abund.meta.data,col_ignore=abund.meta.attr.names)
  file.name.root = paste(file.name.root,label,sep=".")
  
  #write files w/o row names because Excel will incorrectly shift column names due to first empty name
  write.table(abund.meta.data,paste(file.name.root,".count.tab"),sep="\t",row.names = FALSE)
  write.table(abund.meta.data.norm,paste(file.name.root,".freq.tab"),sep="\t",row.names = FALSE)
  
  make.graph.file.name <- function(kind) {
    paste(file.name.root.id.vars,kind,file.name.sfx,sep=".")
  }
  
  m_a = split_count_df(abund.meta.data,col_ignore=abund.meta.attr.names)
  #make.global(m_a)
  
  is.count.data = any(m_a$count > 1)
  
  #to avoid dependency of diversity indices on sequencing depth, we rarefy
  #all samples to the total count taken from the smallest sample,
  #and average over the number of replications of this
  
  if (do.diversity) {
    
    #DEBUG:
    n.rar.rep = 400
    n.rar = min(rowSums(m_a$count))
    
    #S.ACE & all se.* give NaN often
    if (is.count.data) {
      est.r = replicate(n.rar.rep,estimateR(rrarefy(m_a$count,n.rar))[c("S.obs","S.chao1"),])
      #est.r is a 3-d array with (num_indices,num_samples,num_replications); find mean over last dim:
      rich.mat = t(aaply(est.r,1:2,mean))
    }
    else {
      rich.mat = estimateR(m_a$count)[c("S.obs","S.chao1"),]
    }
    rich.meta = cbind(rich.mat,m_a$attr)
    
    #make.global(rich.meta)
    
    f.div = function(m) { cbind(div.shannon=diversity(m,index="shan"),div.simpson = diversity(m,index="simp")) }
    if(is.count.data) {
      div = aaply(replicate(n.rar.rep,f.div(rrarefy(m_a$count,n.rar))),1:2,mean)
    }
    else {
      div = f.div(m_a$count)
    }
    
    if(is.count.data) {
      #this is already unbiased:
      div.unb.simpson = as.matrix(rarefy(m_a$count,2)-1)
      #div.fisher.alpha = fisher.alpha(m_a$count)
    }
    else {
      div.unb.simpson = div$div.simpson
    }
    div.meta = cbind(div,div.unb.simpson,m_a$attr)
    #make.global(div.meta)
  }
  
  for (id.vars in id.vars.list) {
    file.name.root.id.vars = paste0(c(file.name.root,id.vars),collapse="_")
    #coord_trans(y="sqrt"),
    #scale_y_sqrt(),    
    clades.order = NULL
    if (!is.null(names(res.tests)) && !is.null(res.tests$stab.feat)) {
      clades.order = res.tests$stab.feat$labels[res.tests$stab.feat$mean.order]
    }
    pl.hist = plot.abund.meta(abund.meta.data.norm,
                              id.vars=id.vars,
                              attr.names=abund.meta.attr.names,
                              file_name=make.graph.file.name("abund"),
                              ggp.comp=list(
                                scale_y_sqrt(),
                                theme(
                                  plot.title = element_text(size = rel(2)),
                                  axis.text.x = element_text(size = rel(0.85),angle=90),
                                  axis.text.y = element_text(size = rel(0.85))
                                )
                              ),
                              clades.order=clades.order
    )
    #env=as.environment(as.list(environment(), all.names=TRUE))
    #print(names(as.list(env)))
    #print(evals("pl.hist",env=env))
    report$add(pl.hist,
               caption=paste("Abundance profile grouped by ",
                             paste(id.vars,collapse=","),
                             if(!is.null(clades.order)) 
                             {"sorted according to GLMnet feature stability order"} 
                             else {""})
    )
    
    if ( do.diversity ) {
      pl.hist = plot.abund.meta(rich.meta,id.vars=id.vars,attr.names=abund.meta.attr.names,
                                file_name=make.graph.file.name("richness"))
      report$add(pl.hist)
      pl.hist = plot.abund.meta(div.meta,id.vars=id.vars,attr.names=abund.meta.attr.names,
                                file_name=make.graph.file.name("diversity"))
      report$add(pl.hist)
    }
  }
}


dirmult.from.df <- function(data,col_ignore=c()) {
  cnt_m = count_matr_from_df(data,col_ignore=col_ignore)
  #We do not use dirmult() from dirmult package directly
  #because it drops pi elements that are zero, which makes
  #further work with the result very clumsy.
  #However, keep in mind that various methods from HMP package
  #will return NaNs if you pass dirmult parameters that contain
  #zero elements. You will have to filter the parameter object
  #before passing it to such methods.
  dm.par = DM.MoM(cnt_m)
  return (dm.par)
}


dirmult.kelvin <- function(file_name,group_sel) {
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr(file_name)
  data_all = count_filter(data_all,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)  
  dm.par = list()
  for (gr in group_sel) {
    data = data_all[data_all$group == gr,]
    dm.par[[gr]] = dirmult.from.df(data,col_ignore=attr_names)
  }
  return (dm.par)
}

get.clade.names <- function(data,attr.names) {
  colnames(data)[!(colnames(data) %in% attr.names)]
}

power.choc<-function(taxa.meta.data,taxa.meta.attr.names) {
  
  
  dm.par.kelv = dirmult.kelvin("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls",c("Stool"))[[1]]
  
  
  abund_matr_norm = row_normalize_matr(count_matr_from_df(taxa.meta.data,col_ignore=taxa.meta.attr.names))
  
  #assign("abund_matr_norm",abund_matr_norm,envir=globalenv())
  
  #We have only one paired sample. 
  #set pi as proportions from a patient before therapy. Keep theta form HMP dataset.
  dm.par.1 <- dm.par.kelv
  dm.par.1$pi = abund_matr_norm["LK3_91",]
  dm.par.1$gamma = dirmult.comp.gamma(dm.par.1$pi,dm.par.1$theta)
  
  #set pi from the same patient after therapy
  dm.par.2 <- dm.par.kelv
  dm.par.2$pi = abund_matr_norm["LK4_92",]
  dm.par.2$gamma = dirmult.comp.gamma(dm.par.2$pi,dm.par.2$theta)
  
  dm.par.orig = list(dm.par.1,dm.par.2)
  
  groups = as.factor(c("Before","After"))
  
  all_clades = get.clade.names(taxa.meta.data,taxa.meta.attr.names)
  
  #test_clades = as.factor(c("Actinobacteria_0.1.1.1","Clostridia_0.1.6.2"))
  test_clades = as.factor(c("Lachnospiracea_incertae_sedis_0.1.6.2.1.3.7","Faecalibacterium_0.1.6.2.1.5.4","unclassified_0.1.6.2.1.5.9"))
  
  
  n.samp.grp = c(50,60)
  n.seq.range = seq(1500,2000,by=1)
  #n.seq.range <- rep(50,100)
  #n.samp.grp = c(8,5)
  
  #set number of sequences here so that we can estimate
  #effect size and DM test statistics outside of sample generation loop
  n.seq = list()
  
  for (i.group in seq(length(groups))) {
    group = groups[i.group]
    n.seq[[i.group]] <- sample(n.seq.range, size=n.samp.grp[i.group]) 
  }
  
  res.power = power.dirmult.range(
    dm.par.orig=dm.par.orig,
    n.seq=n.seq,
    groups=groups,
    all_clades=all_clades,
    test_clades=test_clades,
    n.samp.grp=n.samp.grp,
    effect.size.range = c(0, 0.1, 0.2, 0.3, 0.4, 1.0),
    n.rep = 100
  )    
  write.csv(res.power,"res.power.csv")
}

read.choc <- function(taxa.level=3) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  moth.taxa <- read.mothur.taxa.summary("d3226.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count_filter(taxa.lev.all,col_ignore=c(),min_max_frac=0.001,min_row_sum=50,other_cnt="other")
  meta = load.meta.choc("JCVI_ALL_sample_information_Feb_2013.AT.txt")
  return (merge.counts.with.meta(taxa.lev,meta))
}

proc.choc.nov_2013_grant_proposal <- function() {
  #DEBUG:
  do.power = F
  do.plots = T
  do.tests = T
  taxa.levels = c(2,3,4,5,6)
  for (taxa.level in taxa.levels) {
    label = paste("16s",taxa.level,sep=".")
    
    taxa.meta = read.choc(taxa.level)
    taxa.meta.data = taxa.meta$data
    taxa.meta.attr.names = taxa.meta$attr.names
    
    #(optionally) subsample just the rows we need and then filter out clades that are all zero
    #taxa.meta.data[taxa.meta.data$Sample.type=="patient",]
    taxa.meta.data = count_filter(taxa.meta.data,
                                  col_ignore=taxa.meta.attr.names,
                                  min_max_frac=0.001,min_row_sum=50,other_cnt="other")
    
    if(do.power) {
      power.choc(taxa.meta.data,taxa.meta.attr.names)
    }
    if (do.tests) {
      res.tests = try(
        test.counts.choc(taxa.meta.data,taxa.meta.attr.names,
                         label=label,
                         stability.transform.counts="ihs",
                         do.stability=T,
                         do.tests=F)
      )
    }
    if(do.plots) {
      try(
        std.plots(taxa.meta.data,taxa.meta.attr.names,id.vars.list=
                    list(
                      c("Sample.type","Therapy.Status"),
                      c("Subject.ID","Therapy.Status"),
                      c("Sample.type.1"),
                      c("Sample.ID.1")
                    ),
                  label=label,
                  res.tests=res.tests
        )
      )
    }
  }
}

proc.choc <- function() {
  taxa.levels = c(2,3,4,5,6)
  #taxa.levels = c(3)
  do.std.plots = T
  do.tests = T
  
  report$add.descr("Largely identical set of analysis routines is applied
                   in loops over different 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out the taxonomic level.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture.")
  
  for (taxa.level in taxa.levels) {
    taxa.meta = read.choc(taxa.level)
    
    label = paste("16s","l",taxa.level,sep=".",collapse=".")
    report$add.p(paste("Tag:",label))
    report$set.tag(label)
    
    report$add.header(paste("Taxonomic level:",taxa.level),2)
    
    make.global(taxa.meta)
    #taxa.meta$data = taxa.meta$data[taxa.meta$data$Sample.type.1 != "sibling",]
    
    taxa.meta.aggr = taxa.meta
    
    aggr_var = "SubjectID"
    
    taxa.meta.aggr = aggregate_by_meta_data(taxa.meta$data,
                                            aggr_var,
                                            taxa.meta$attr.names)
    report$add.p(paste("After aggregating samples by ", aggr_var, ":",nrow(taxa.meta.aggr$data)))
    
    report$add.p(paste("Number of samples:",nrow(taxa.meta.aggr$data)))      
    
    taxa.meta.aggr$data = count_filter(taxa.meta.aggr$data,
                                       col_ignore=taxa.meta.aggr$attr.names,
                                       min_max=10)
    #taxa.meta.data = count_filter(taxa.meta.data,col_ignore=taxa.meta.attr.names,
    #                              min_median=0.002,
    #                              min_max_frac=0.1,min_max=10,
    #                              min_row_sum=500,other_cnt="other")
    
    report$add.p(paste("After count filtering,",
                       (ncol(taxa.meta.aggr$data)-length(taxa.meta.aggr$attr.names)),"clades left."))
    
    xtabs.formulas = list(~Sample.type+TherapyStatus,~FamilyID,~Sample.type.1,~SubjectID)
    for(xtabs.formula in xtabs.formulas) {
      fact.xtabs = xtabs(xtabs.formula,data=taxa.meta.aggr$data,drop.unused.levels=T)
      report$add.table(fact.xtabs,show.row.names=T,caption=paste("Sample cross tabulation",xtabs.formula))
      report$add.printed(summary(fact.xtabs))
    }
    
    res.tests = NULL
    if (do.tests) {
      res.tests = tryCatchAndWarn(
        test.counts.choc(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,
                         label=label,
                         stability.transform.counts="ihs",
                         do.stability=T,
                         do.tests=T,
                         do.genesel=T,do.glmnet=T,
                         do.glmer=T,do.adonis=T,
        )
      )
    }
    if (do.std.plots) {
      tryCatchAndWarn({
        plot.group = list(
          c("Sample.type","TherapyStatus"),
          #c("SubjectID","TherapyStatus"),
          c("Sample.type.1")
          #c("SampleID.1")
        )            
        std.plots(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,id.vars.list=
                    plot.group,
                  label=label,
                  res.tests=res.tests
        )
      })
      tryCatchAndWarn({
        heatmap.choc(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,label=label)
      })
    }
  }
}


#Function for symbolic derivative of arbitrary order
#Copied from help page for 'deriv()'
#Use as:
#dd.expr = DD(expression(log(x+(x**2+1)**0.5)),"x",2)
#To create a function from the 'expression' output:
#f = function(x) eval(dd.expr)
DD <- function(expr, name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}

## Return a function that computes the derivative of expr
make.DD <- function(expr,name,order = 1) {
  dn.expr = DD(expr,name,order=order)
  return (function(x) {eval(dn.expr)})
}

##Delta method to approximate mean and variance of
##a function of random variable
##(http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables)
##See ref for the caviats
mom.f <- function(var.x,mean.x,f,f.d1,f.d2) {
  mean.f = f(mean.x) + f.d2(mean.x)/2*var.x
  var.f = f.d1(mean.x)**2*var.x
  c(mean.f=mean.f,var.f=var.f,sd.f=var.f**0.5)
}

cohens.d.from.mom <- function(mean.gr,var.gr,n.gr) {
  lx <- n.gr[1] - 1
  ly <- n.gr[2] - 1
  md  <- abs(mean.gr[1] - mean.gr[2])        ## mean difference (numerator)
  csd <- lx * var.gr[1] + ly * var.gr[2]
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  md/csd                        ## cohen's d
}

## IHS (inverse hyperbolic sign) transform
## This is the same as log(x+(x**2+1)**0.5)
ihs <- function(x,theta=1) {
  asinh(theta*x)/theta
}

ihs.d1 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",1)
ihs.d2 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",2)

## Boxcox transform
boxcox <- function(x,lambda1,lambda2=0) {
  #print(paste("l1=",lambda1,"l2=",lambda2))
  if (lambda1 != 0) {
    ((x+lambda2)**lambda1-1)/lambda1
  }
  else {
    log(x+lambda2)
  }
}

## Fit boxcox transform to the data and transform the data
boxcox.transform.vec <- function(x) {
  l.2 = T
  if(!all(x>0)) {
    b = boxcoxfit(x,lambda2=T)
  }
  else {
    b = boxcoxfit(x)
    l.2 = F
  }
  #b = try(boxcoxfit(x,lambda2=T),TRUE)
  #for some reason the above fails if all(b>0)
  #we just repeat it with lambda2 NULL
  #if (inherits(b, "try-error")) {
  #  b = boxcoxfit(x)
  #  l.2 = T
  #  print("err")
  #  print(b)
  #}
  lambda1 = NA
  lambda2 = 0
  if (l.2) { 
    lambda1 = b$lambda["lambda"] 
    lambda2 = b$lambda["lambda2"] 
  }
  else {
    lambda1 = b$lambda
  }
  return (list(x=boxcox(x,lambda1=lambda1,lambda2=lambda2),
               boxcox=b))
}

boxcox.transform.mat <- function(m) {
  foreach(x=iter(m,by="col"),.combine=cbind) %do% (boxcox.transform.vec(x)$x)
}

show.distr <- function(x) {
  #because aes leaves its expressions unevaluated, we need
  #to bind the value of x as data frame parameter of ggplot
  ggplot(data.frame(x=x),aes(x))+
    geom_histogram()
  #geom_density(aes(y = ..density..))
  #stat_density()
  #hist(x,
  #     freq = F,
  #     breaks = "FD",      # For more breaks than the default
  #     col = "darkslategray4", border = "seashell3")
  #lines(density(x),
  #      col = "firebrick2", lwd = 3)
}

show.distr.group <- function(x,group) {
  #because aes leaves its expressions unevaluated, we need
  #to bind the value of x as data frame parameter of ggplot
  ggplot(data.frame(x=x,group=group),aes(x=x,fill=group))+
    #geom_histogram(alpha=0.2, position="identity")+
    geom_density(alpha=0.2,aes(y = ..density..))
}

show.trend <- function(meta.data,x,y,group,title="Group trends") {
  df = meta.data[,c(x,y,group)]
  names(df) = c("x","y","group")
  #df$x = as.Date(timeDate(df$x))  
  df$y = boxcox.transform.vec(df$y)$x
  print(summary(lm(y~x,df)))
  smooth_method = "loess" #"lm"
  print(
    ggplot(df, aes(x=x, y=y,color=group)) +
      geom_point() +
      #geom_line(alpha=0.3, linetype=3) + 
      #geom_smooth(aes(group=group,color=group), method='lm', formula=y~x+I(x^2)+I(x^3)) + 
      stat_smooth(method=smooth_method, se = T,degree=1) + 
      #scale_x_date() +
      labs(title=title)
  )
}

##Code adapted from Jyoti Shankar. Find best alpha through cross-validation as
##per the recipe from cv.glmnet help page
cv.glmnet.alpha <- function(y, x, family, seed=NULL, standardize=T) {
  # This seed is taken from LASSO's example. 
  if (!is.null(seed) ) { set.seed(seed) }
  # Setting the number of folds to the number of samples (leave one out)
  #is not recommended by cv.glmnet help page
  numfolds <- min(15,dim(x)[1])
  # Grid for alpha crossvalidation
  alphas <- c(0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999)
  # In a k fold crossvalidation, fold ID gives the iteration in which the sample would be in the test set.
  
  foldid <- sample(rep(seq(numfolds),length=dim(x)[1]))
  # Go through alpha grid
  # Run crossvalidation for lambda.
  # Each model for each alpha is run by a parallel core and put into a list called lassomodels
  lassomodels <- foreach(i = c(1:length(alphas))) %dopar%
{
  #re-import required for windows parallelism with doSNOW
  library(glmnet)
  # set.seed(seed)
  # the function finds the best lambda for a given alpha
  # within each model there is cross-validation happening for lambda for each alpha.
  # lambda1 = lambda*alpha 
  # lambda2 = lambda*(1-alpha)
  model <- try(cv.glmnet(x=x, y=y, family=family,
                         nfolds=numfolds, 
                         type.measure="deviance", 
                         foldid=foldid,
                         standardize=standardize, 
                         alpha=alphas[i]))
}
# there are two lambdas per model
# minimum lamda
# lambda within 1 standard deviation of the error
# find the best alpha
best_alpha_index <- 0
lowest_error <- 0
for (i in c(1:length(alphas))) {
  # if a model fails, "try-error" will return true. "try-error" is an object that traps errors
  # inherits is a function that will be true if try-error has collected an error from the model
  # we want to avoid any errors recorded in "try-error" in the list of models we just generated
  # example too small a dataset
  if (!inherits(lassomodels[[i]], "try-error")) {
    # First we will find the index of the lambda corresponding to the lambda.min
    index <- which(lassomodels[[i]]$lambda.min == lassomodels[[i]]$lambda)
    # high lambda means more penalty.
    # lambdas are arranged from highest to lowest
    # alpha = 1 ==> lasso
    # alph = 0 ==> ridge
    
    #cvm is the cross-validated error, in this case, deviance.
    error <- lassomodels[[i]]$cvm[index]
    #print(error)
    if (best_alpha_index == 0 || error < lowest_error) {
      best_alpha_index <- i # picks an alpha from the grid of alphas
      lowest_error <- error # picks the lowest deviance from the grid
    }
  }
}
#print(best_alpha_index)
out <- list(c(), c(), c())
if (best_alpha_index != 0) {
  # print the lassomodel at the best_alpha_index
  lasso_model <- lassomodels[[best_alpha_index]]
  alpha <- alphas[best_alpha_index]
  # Use lambda which gives the lowest cross validated error
  lambda <- lasso_model$lambda.min
  out <- list(lasso_model, alpha, lambda)
}
names(out) <- c("glmnet.model", "alpha", "lambda")
out
}

stability.selection.c060.at <- function (x, fwer, pi_thr = 0.6) 
{
  stopifnot(pi_thr > 0.5, pi_thr < 1)
  if (class(x$fit)[1] == "multnet") {
    p <- dim(x$fit$beta[[1]])[1]
  }
  else {
    p <- dim(x$fit$beta)[1]
  }
  qv <- ceiling(sqrt(fwer * (2 * pi_thr - 1) * p))
  lpos <- which(x$qs > qv)[1]
  if (!is.na(lpos)) {
    stable <- which(x$x[, lpos] >= pi_thr)
    stable.p <- x$x[stable, lpos]
    lpos.order <- order(x$x[, lpos],decreasing = TRUE)
    lpos.order.p <- x$x[lpos.order,lpos]
  }
  else {
    stable <- NULL
    lpos.order <- NULL
  }
  mean.p = rowMeans(x$x)
  mean.order = order(mean.p,decreasing = TRUE)
  mean.order.p = mean.p[mean.order]
  #Can use factor(z,levels=labels[mean.order],ordered=T)
  #to order some z by the ranking of variables
  #returned by this method
  out <- list(stable = stable, 
              stable.p = stable.p,
              lpos.order = lpos.order,
              lpos.order.p = lpos.order.p,
              mean.order = mean.order,
              mean.order.p = mean.order.p,
              labels = rownames(x$x),
              lambda = x$fit$lambda[lpos], 
              lpos = lpos, fwer = fwer,
              stab.path=x)
  return(out)
}

plot.stability.selection.c060.at <- function (y, 
                                              annot = TRUE, 
                                              main = "Stability path", 
                                              log.scale = TRUE, 
                                              nvar = 8,
                                              xvar="lambda",
                                              rank="lpos") 
{
  x = y$stab.path
  prob = t(x$x)
  nzeros <- which(colSums(prob) != 0)  
  p = ncol(prob)
  
  rank.order = switch(rank,lpos=y$lpos.order,mean=y$mean.order)
  selection <- rep("unselected", ncol(prob))
  names(selection) = colnames(prob)
  labels = names(selection)
  rank.sel = paste("Top",nvar,"by rank")
  fwer.sel = paste("FWER <",y$fwer)
  selection[rank.order[1:nvar]] = rank.sel
  if(!is.null(y$stable)) {
    selection[y$stable] = fwer.sel
  }
  
  xv <- switch(xvar, fraction = x$fit$lambda/max(x$fit$lambda), 
               x$fit$lambda)
  if (xvar == "lambda") {
    if (log.scale == T) {
      xv <- log10(xv)
    }
  }
  data.coef <- melt(data.frame(xvar = xv, prob = t(x$x)), 
                    id="xvar")
  data.coef$selection <- factor(rep(selection, each = length(xv)))
  data.coef$labels <- factor(rep(names(selection), each = length(xv)),
                             levels=labels[rank.order],
                             ordered=T)
  colnames(data.coef) <- c("xvar", "var", "prob", "Selection", 
                           "Variables")
  linetype.values = c(unselected = "dotted")
  linetype.values[rank.sel] = "dashed"
  linetype.values[fwer.sel] = "solid"
  #data.coef$prob = log10(data.coef$prob)
  #make.global(data.coef)
  d <- ggplot() + 
    #geom_line(aes(x = xvar, 
    #              y = prob)) + 
    geom_line(data=data.coef[data.coef$Selection == "unselected",],
              aes(x = xvar, y = prob, group=var),
              color="black")+    
    geom_line(data=data.coef[data.coef$Selection != "unselected",],
              aes(x = xvar, y = prob, group=var,
                  linetype = Selection, 
                  colour = Variables),
              size=1) +
    scale_linetype_manual(values = linetype.values)+
    labs(x = switch(xvar, fraction = expression(lambda[1]/max[lambda]), 
                    ifelse(log.scale, expression(log[10](lambda)), 
                           expression(lambda[1]))), y = "Selection probabilities") + 
    ggtitle(main) +
    coord_fixed()
  #d <- d + scale_x_reverse()
  if (is.null(labels)) {
    d <- d + theme(legend.position = "none")
  }
  else {
    if (length(labels[nzeros]) != length(nzeros)) {
      d <- d + theme(legend.position = "none")
    }
  }
  return (d)
}


test.counts.t1d <- function(data,attr.names,label,alpha=0.05,
                            do.tests=T,do.stability=T,
                            do.genesel=T,do.glmnet=T,
                            do.glmer=T,do.adonis=T,
                            stability.transform.counts="ihs") {
  
  n.batch.levels = sum(table(data$Batch) > 0)
  stability.resp.attr = "T1D"
  stability.model.family = "binomial"
  n.adonis.perm = 4000
  res = list()
  #data = data[data$Batch %in% c(1,2,3),]
  
  ##only families with 1 sibling and 1 patient
  #tbl = with(data,table(Family,T1D))
  #fams = unique(row.names(tbl[(tbl[,"T1D"]>0 & tbl[,"Control"]>0),]))
  #data = data[data$Family %in% fams,]
  
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  
  data.norm = row_normalize(data,col_ignore=attr.names)
  
  #make.global(tbl)
  #make.global(fams)
  #make.global(data.norm)
  
  all.clades = get.clade.names(data,attr.names)
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  
  m = m_a$count
  
  #m = (m_a.abs$count > 0)
  #storage.mode(m) = "integer"
  
  count = switch(stability.transform.counts,
                 boxcox=boxcox.transform.mat(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m_a.abs$count > 0))
  
  if (do.stability) {
    
    if(do.genesel) {
      ## We need to use row-normalized counts here, because we perform Wilcox test
      ## inside, and could false significance due to different depth of sequencing
      res.stab_sel_genesel = stab_sel_genesel(m_a$count,m_a$attr[,stability.resp.attr])
      report$add.header("GeneSelector stability ranking",3)
      report$add.package.citation("GeneSelector")
      report$add.descr("Univariate test (Wilcox) is applied to each clade on random
                     subsamples of the data. Consensus ranking is found with a
                     Monte Carlo procedure. The top clades of the consensus ranking
                     are returned, along with the p-values computed on the full
                     original dataset (with multiple testing correction).")
      report$add.table(res.stab_sel_genesel$stab_feat,
                       caption=
                         paste("GeneSelector stability ranking for response ",stability.resp.attr)
      )
      report$add.package.citation("vegan")
      report$add(
        plot.features.mds(m_a,species.sel=(colnames(m_a$count) %in% 
                                             res.stab_sel_genesel$stab_feat$name),
                          sample.group=m_a$attr[,stability.resp.attr]),
        caption="metaMDS plot. 'x' marks top selected clades, 'o' marks samples"
      )
    }
    if(do.glmnet) {
      standardize.glm = T
      if(standardize.glm) {
        count = decostand(count,method="standardize",MARGIN=2)
      }
      
      cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
      registerDoSNOW(cl)  
      cv.res = cv.glmnet.alpha(m_a$attr[,stability.resp.attr],count,family=stability.model.family,standardize=F)
      stopCluster(cl)
      
      penalty.alpha = cv.res$alpha
      #alpha = 0.8
      stab.res.c060 = stabpath(m_a$attr[,stability.resp.attr],count,weakness=0.8,
                               family=stability.model.family,steps=600,
                               alpha=penalty.alpha,standardize=F)
      #stab.res.c060 = stabpath(m_a$attr$A1C[m_a$attr$T1D=="T1D"],
      #                               count[m_a$attr$T1D=="T1D",],
      #                               weakness=0.9,
      #                               family="gaussian",steps=600,
      #                               alpha=penalty.alpha,standardize=F)
      
      fwer = alpha
      pi_thr = 0.6
      stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
      #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
      
      stab.path.file = paste("stability.path",label,"png",sep=".")
      
      pl.stab = tryCatchAndWarn({
        #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
        p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
        ggsave(stab.path.file)
        p
      })
      
      #report$add.printed(stab.res.c060$fit,
      #                   caption=paste(
      #                     "Glmnet stability path analysis for response (",
      #                     stability.resp.attr,
      #                     ")"
      #                     )
      #)
      report$add.header(paste(
        "Glmnet stability path analysis for response (",
        stability.resp.attr,
        ")"
      ),3
      )
      report$add.package.citation("c060")
      report$add.descr("This multivariate feature selection method builds glmnet models 
                     for a given response variable at 
                     varying strengths of L1 norm
                     regularization parameter and on multiple random subsamples at each
                     level of the parameter. The features (e.g. taxonomic clades)
                     are ranked according to their probability to be included in the
                     model")
      report$add.printed(stab.feat.c060$stable,
                         caption="Features that passed FWER in stability analysis:")
      
      if(!is.null(pl.stab)) {
        report$add(pl.stab,caption="Stability path. Probability of each variable 
                 to be included in the model as a function of L1 regularization
                 strength. Paths for top ranked variables are colored. The variables
                 (if any) that passed family wide error rate (FWER) cutoff are
                 plotted as solid lines.")
      }
      
      res$stab.feat=stab.feat.c060
    }
  }
  
  if (do.tests) {
    
    if(do.glmer) {
      ##intercept varying among families and among individuals within families (nested random effects)
      formula_rhs = "T1D + (1|FamilyID/SubjectID)"
      ##if there are repeated samples per individual, add sample random effect for overdispersion
      ##otherwise, overdispersion is taken care of by SubjectID random effect(?)
      if(anyDuplicated(m_a.abs$attr$SubjectID)) {
        formula_rhs = paste(formula_rhs,"+(1|SampleID)",sep="")
      }
      
      if(n.batch.levels>1) {
        ## different slope and intercept for T1D across batches as random effect
        formula_rhs = paste(formula_rhs,"+(1+T1D|Batch)",sep="")
        linfct=c("T1DT1D = 0")
        #linfct=c("T1D1 = 0")
      }
      else {
        formula_rhs = paste(formula_rhs,"+TimestampMonth",sep="")
        linfct=c("T1DT1D = 0")
      }
      res$glmer = test_taxa_count_glmer(m_a.abs,alpha=alpha,
                                        formula_rhs=formula_rhs,
                                        linfct=linfct)
      report$add.header("Binomial mixed model analysis",3)
      report$add.package.citation("lme4")
      report$add.descr("The random effect terms in the model formula describe that: 
                         - intercept varies among families 
                           and among individuals within families (nested random effects);
                         - when there are multiple observations (samples) per
                           individual, we add a random effect for each observation to account
                           for the overdispersion;
                       The fixed effect is T1D status. The binomial family
                       is used to build a set of univariate models, with each
                       model describing the observed counts of one clade.
                       P-values are estimated from the model under a null hypothesis
                       of zero coefficients and a two-sided alternative. 
                       Benjamini & Hochberg (1995) method is used 
                       for multiple testing correction, and the significant clades
                       are reported.")
      report.taxa_count_glmer(report,res$glmer)
    }
    if(do.adonis) {
      ##Negative values break bray-curtis and jaccard distances; we standardize to "range" to reduce
      ##the influence of very abundant species:
      count = decostand(m_a$count,method="range",MARGIN=2)
      #or, no standartization:
      #count = m_a$count
      adonis.dist = "bray" #"jaccard" #"bray"
      #ad.res = adonis(count~T1D + Batch,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
      #print(ad.res)
      report$add.header("PermANOVA (adonis) analysis of taxonomic profiles",3)
      report$add.package.citation("vegan")
      report$add.descr("Non-parametric multivariate test for association between
                     taxonomic profiles and meta-data variables. Profile is normalized
                     to proportions across clades and then to range across samples.")
      
      if(n.batch.levels > 1) {
        ad.res.batch.t1d = adonis(count~Batch + T1D,data=m_a$attr,
                                  permutations=n.adonis.perm,method=adonis.dist) 
        report$add.printed(ad.res.batch.t1d,caption="Association with Batch")
      }
      ad.res.t1d = adonis(count~T1D,data=m_a$attr,
                          permutations=n.adonis.perm,method=adonis.dist)
      
      report$add.printed(ad.res.t1d,"Association with T1D unpaired")
      
      ad.res.t1d.paired = adonis(count~T1D,data=m_a$attr,strata=m_a$attr$FamilyID,
                                 permutations=n.adonis.perm,method=adonis.dist)
      
      report$add.printed(ad.res.t1d.paired,"Association with T1D matched within families")
      
      report$add.printed(adonis(count~age,data=m_a$attr,
                                permutations=n.adonis.perm,method=adonis.dist),
                         "Association with age unpaired")
      
      report$add.printed(adonis(count~Timestamp,data=m_a$attr,
                                permutations=n.adonis.perm,method=adonis.dist),
                         "Association with collection date unpaired")      
      
      report$add.printed(adonis(count~age+T1D,data=m_a$attr,
                                permutations=n.adonis.perm,method=adonis.dist),
                         "Association with age and T1D unpaired")            
      
      report$add.printed(adonis(count~T1D+age,data=m_a$attr,
                                permutations=n.adonis.perm,method=adonis.dist),
                         "Association with T1D and age unpaired")            
      
      
      if(F) {
        ## We filter rows and need to do range transformation on the subset again
        ##TODO: we also need to filter first clades with zero count
        count_a1c = m_a$count[!is.na(m_a$attr$A1C),]
        count_a1c = decostand(count_t1d,method="range",MARGIN=2)
        meta_a1c = m_a$attr[!is.na(m_a$attr$A1C),]
        ad.res.a1c = adonis(count_a1c~A1C,
                            data=meta_a1c,
                            permutations=n.adonis.perm,method=adonis.dist)
        print ("Results of adonis:")
        print(ad.res.a1c)
      }
    }
    #print (ad.res)
    #test.ad.res = ad.res$aov.tab$"Pr(>F)"[1]
    #r2.ad = ad.res$aov.tab$R2[1]
    if(F) {
      gr.abs = m_a.abs$attr$T1D
      gr.abs.cnt = list(x=m_a.abs$count[gr.abs=="T1D",],y=m_a.abs$count[gr.abs!="T1D",])
      print(Xdc.sevsample(gr.abs.cnt))
      K = dim(gr.abs.cnt[[1]])[2]-1
      print (K)
      #make.global(gr.abs.cnt)
      print(Xmcupo.sevsample(gr.abs.cnt,K))
      
      ##this is a rank based pairwise test; column-wise standartization does not
      ##affect it, but row-wise - does. We do need row-normalized data
      gr = m_a$attr$T1D
      gr.cnt = list(x=m_a$count[gr=="T1D",],y=m_a$count[gr!="T1D",])
      
      test.wil.res = list()
      
      for (clade in all.clades) {
        test.wil.res[[clade]] = wilcox.test(gr.cnt[[1]][,clade],gr.cnt[[2]][,clade],paired=FALSE,exact=F)$p.value
      }
      test.wil.res = unlist(test.wil.res)
      make.global(test.wil.res)    
      names(test.wil.res) = all.clades
      test.wil.res.adj = p.adjust(test.wil.res,method="BH")
      fdr.res = fdrtool(test.wil.res,"pvalue",plot=F,verbose=F)
      names(test.wil.res.adj) = names(test.wil.res)
      print("Unpaired Wilcox test for clade-T1D")
      print("Before corection:")
      print(test.wil.res[test.wil.res<=alpha])
      print("After BH corection:")
      print(test.wil.res.adj[test.wil.res.adj<=alpha])
      make.global(test.wil.res.adj)
      test.wil.res.adj.q = fdr.res$qval
      names(test.wil.res.adj.q) = names(test.wil.res)
      print("Q-values after FDR corection:")
      print(test.wil.res.adj.q[test.wil.res.adj.q<=alpha])
      make.global(test.wil.res.adj.q)
      
      if(F) {
        for (clade in all.clades) {
          test.wil.res[[clade]] = cor.test(count_a1c[,clade],meta_a1c$A1C,
                                           method="spearman")$p.value
        }
        test.wil.res = unlist(test.wil.res)
        make.global(test.wil.res)    
        names(test.wil.res) = all.clades
        test.wil.res.adj = p.adjust(test.wil.res,method="BH")
        fdr.res = fdrtool(test.wil.res,"pvalue",plot=F,verbose=F)
        names(test.wil.res.adj) = names(test.wil.res)
        print("Spearman RHO test for clade-A1C")
        print("Before corection:")
        print(test.wil.res[test.wil.res<=alpha])
        print("After BH corection:")
        print(test.wil.res.adj[test.wil.res.adj<=alpha])
      }
    }
  }
  return (res)
}


test.counts.choc.Nov_2013_gran_proposal <- function(data,attr.names,label,alpha=0.05,
                                                    do.tests=T,do.stability=T,
                                                    stability.transform.counts="ihs") {
  
  n.adonis.perm = 400
  res = list()
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  
  data.norm = row_normalize(data,col_ignore=attr.names)
  
  #make.global(tbl)
  #make.global(fams)
  #make.global(data.norm)
  
  all.clades = get.clade.names(data,attr.names)
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  make.global(m_a)
  
  m = m_a$count
  
  #m = (m_a.abs$count > 0)
  #storage.mode(m) = "integer"
  
  count = switch(stability.transform.counts,
                 boxcox=boxcox.transform.mat(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m_a.abs$count > 0))
  make.global(count)
  
  if (do.stability) {
    
    standardize.glm = T
    
    cl<-makeCluster(getOption("mc.cores", 4L)) #number of CPU cores
    registerDoSNOW(cl)
    stab.resp.var = m_a$attr$Therapy.Status
    #stab.resp.var = m_a$attr$Sample.type
    cv.res = cv.glmnet.alpha(stab.resp.var,count,family="binomial",standardize=standardize.glm)
    stopCluster(cl)
    
    make.global(cv.res)
    
    penalty.alpha = cv.res$alpha
    #alpha = 0.8
    stab.res.c060 = stabpath(stab.resp.var,count,weakness=0.6,
                             family="binomial",steps=600,
                             alpha=penalty.alpha,standardize=standardize.glm)
    make.global(stab.res.c060)
    fwer = alpha
    pi_thr = 0.6
    stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    make.global(stab.feat.c060)
    report$add.printed(names(stab.feat.c060$stable))
    
    stab.path.file = paste("stability.path",label,"png",sep=".")
    
    #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
    p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
    ggsave(stab.path.file)
    #png(stab.path.file)
    #plot(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    #dev.off()
    res$stab.feat=stab.feat.c060
  }
  
  if (do.tests) {
    ##Negative values break bray-curtis and jaccarda= distances; we standardize to "range" to reduce
    ##the influence of very abundant species
    count = decostand(m_a$count,method="range",MARGIN=2)
    adonis.dist = "jaccard" #"bray"
    #ad.res = adonis(count~T1D + Batch,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    #print(ad.res)
    ad.res.unpaired = adonis(count~Therapy.Status,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    print(ad.res.unpaired)
    ad.res.paired = adonis(count~Therapy.Status,data=m_a$attr,strata=m_a$attr$Family,permutations=n.adonis.perm,method=adonis.dist)  
    print (ad.res.paired)
    #make.global(ad.res)
    
    #print (ad.res)
    #test.ad.res = ad.res$aov.tab$"Pr(>F)"[1]
    #r2.ad = ad.res$aov.tab$R2[1]
    
    gr.abs = m_a.abs$attr$Therapy.Status
    gr.abs.cnt = list(x=m_a.abs$count[gr.abs=="before.chemo",],y=m_a.abs$count[gr.abs!="before.chemo",])
    print(Xdc.sevsample(gr.abs.cnt))
    K = dim(gr.abs.cnt[[1]])[2]-1
    print (K)
    #make.global(gr.abs.cnt)
    print(Xmcupo.sevsample(gr.abs.cnt,K))
    
    ##this is a rank based pairwise test; column-wise standartization does not
    ##affect it, but rwo-wise - does
    gr = m_a$attr$Therapy.Status
    gr.cnt = list(x=m_a$count[gr=="before.chemo",],y=m_a$count[gr!="before.chemo",])
    
    test.wil.res = list()
    
    for (clade in all.clades) {
      test.wil.res[[clade]] = wilcox.test(gr.cnt[[1]][,clade],gr.cnt[[2]][,clade],paired=FALSE,exact=F)$p.value
    }
    test.wil.res = unlist(test.wil.res)
    make.global(test.wil.res)    
    names(test.wil.res) = all.clades
    test.wil.res.adj = p.adjust(test.wil.res,method="BH")
    fdr.res = fdrtool(test.wil.res,"pvalue")
    
    names(test.wil.res.adj) = names(test.wil.res)
    print("Significant unadjusted p-values:")
    print(test.wil.res[test.wil.res<=alpha])
    print("Significant BH adjusted p-values:")    
    print(test.wil.res.adj[test.wil.res.adj<=alpha])
    make.global(test.wil.res.adj)
    test.wil.res.adj.q = fdr.res$qval
    names(test.wil.res.adj.q) = names(test.wil.res)
    print("Q-values below cutoff:")
    print(test.wil.res.adj.q[test.wil.res.adj.q<=alpha])
    make.global(test.wil.res.adj.q)
    
    
  }
  return (res)
}

test.counts.choc <- function(data,attr.names,label,alpha=0.05,
                             do.tests=T,do.stability=T,
                             do.genesel=T,do.glmnet=T,
                             do.glmer=T,do.adonis=T,
                             stability.transform.counts="ihs") {
  
  if(F)
    stability.resp.attr = "TherapyStatus"
  else
    stability.resp.attr = "Sample.type"
  stability.model.family = "binomial"
  n.adonis.perm = 4000
  if(do.adonis) {
    if(F) {
      adonis.tasks = list(
        list(formula_rhs="TherapyStatus",
             strata=NULL,
             descr="Association with the therapy status unpaired"),
        list(formula_rhs="TherapyStatus",
             strata="FamilyID",
             descr="Association with the therapy status paired")
      )
    }
    else {
      adonis.tasks = list(
        list(formula_rhs="Sample.type",
             strata=NULL,
             descr="Association with the patient/sibling status unpaired"),
        list(formula_rhs="Sample.type",
             strata="FamilyID",
             descr="Association with the patient/sibling status paired")
      )
    }
  }
  if(do.glmer) {
    ##intercept varying among families and among individuals within families (nested random effects)
    if(F) {
      formula_rhs = "TherapyStatus + (1|FamilyID/SubjectID)"
      linfct=c("TherapyStatusbefore.chemo = 0")
    }
    else {
      formula_rhs = "Sample.type + (1|FamilyID/SubjectID)"
      linfct=c("Sample.typesibling = 0")
    }
    ##if there are repeated samples per individual, add sample random effect for overdispersion
    ##otherwise, overdispersion is taken care of by SubjectID random effect(?)
    if(anyDuplicated(data$SubjectID)) {
      formula_rhs = paste(formula_rhs,"+(1|SampleID)",sep="")
    }
    
    glmer_descr = "The random effect terms in the model formula describe that: 
                       - intercept varies among families 
                       and among individuals within families (nested random effects);
                       - when there are multiple observations (samples) per
                       individual, we add a random effect for each observation to account
                       for the overdispersion;
                       The fixed effect is T1D status. The binomial family
                       is used to build a set of univariate models, with each
                       model describing the observed counts of one clade.
                       P-values are estimated from the model under a null hypothesis
                       of zero coefficients and a two-sided alternative. 
                       Benjamini & Hochberg (1995) method is used 
                       for multiple testing correction, and the significant clades
                       are reported."
  }  
  res = list()
  
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  make.global(m_a.abs)
  
  row_summ = summary(rowSums(m_a.abs$count))
  report$add.printed(row_summ,caption="Summary of absolute counts")
  
  data.norm = row_normalize(data,col_ignore=attr.names)
  
  #make.global(data.norm)
  
  all.clades = get.clade.names(data,attr.names)
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  
  make.global(m_a)
  
  m = m_a$count
  
  #m = (m_a.abs$count > 0)
  #storage.mode(m) = "integer"
  
  count = switch(stability.transform.counts,
                 boxcox=boxcox.transform.mat(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m_a.abs$count > 0))
  
  if (do.stability) {
    
    if(do.genesel) {
      ## We need to use row-normalized counts here, because we perform Wilcox test
      ## inside, and could false significance due to different depth of sequencing
      res.stab_sel_genesel = stab_sel_genesel(m_a$count,m_a$attr[,stability.resp.attr])
      report$add.header("GeneSelector stability ranking",3)
      report$add.package.citation("GeneSelector")
      report$add.descr("Univariate test (Wilcox) is applied to each clade on random
                       subsamples of the data. Consensus ranking is found with a
                       Monte Carlo procedure. The top clades of the consensus ranking
                       are returned, along with the p-values computed on the full
                       original dataset (with multiple testing correction).")
      report$add.table(res.stab_sel_genesel$stab_feat,
                       caption=
                         paste("GeneSelector stability ranking for response ",stability.resp.attr)
      )
      report$add.package.citation("vegan")
      report$add(
        plot.features.mds(m_a,species.sel=(colnames(m_a$count) %in% 
                                             res.stab_sel_genesel$stab_feat$name),
                          sample.group=m_a$attr[,stability.resp.attr]),
        caption="metaMDS plot. 'x' marks top selected clades, 'o' marks samples"
      )
    }
    if(do.glmnet) {
      standardize.glm = T
      if(standardize.glm) {
        count = decostand(count,method="standardize",MARGIN=2)
      }
      
      cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
      registerDoSNOW(cl)  
      cv.res = cv.glmnet.alpha(m_a$attr[,stability.resp.attr],count,family=stability.model.family,standardize=F)
      stopCluster(cl)
      
      penalty.alpha = cv.res$alpha
      #alpha = 0.8
      stab.res.c060 = stabpath(m_a$attr[,stability.resp.attr],count,weakness=0.8,
                               family=stability.model.family,steps=600,
                               alpha=penalty.alpha,standardize=F)
      #stab.res.c060 = stabpath(m_a$attr$A1C[m_a$attr$T1D=="T1D"],
      #                               count[m_a$attr$T1D=="T1D",],
      #                               weakness=0.9,
      #                               family="gaussian",steps=600,
      #                               alpha=penalty.alpha,standardize=F)
      
      fwer = alpha
      pi_thr = 0.6
      stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
      #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
      
      stab.path.file = paste("stability.path",label,"png",sep=".")
      
      pl.stab = tryCatchAndWarn({
        #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
        p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
        ggsave(stab.path.file)
        p
      })
      
      #report$add.printed(stab.res.c060$fit,
      #                   caption=paste(
      #                     "Glmnet stability path analysis for response (",
      #                     stability.resp.attr,
      #                     ")"
      #                     )
      #)
      report$add.header(paste(
        "Glmnet stability path analysis for response (",
        stability.resp.attr,
        ")"
      ),3
      )
      report$add.package.citation("c060")
      report$add.descr("This multivariate feature selection method builds glmnet models 
                       for a given response variable at 
                       varying strengths of L1 norm
                       regularization parameter and on multiple random subsamples at each
                       level of the parameter. The features (e.g. taxonomic clades)
                       are ranked according to their probability to be included in the
                       model")
      report$add.printed(stab.feat.c060$stable,
                         caption="Features that passed FWER in stability analysis:")
      
      if(!is.null(pl.stab)) {
        report$add(pl.stab,caption="Stability path. Probability of each variable 
                   to be included in the model as a function of L1 regularization
                   strength. Paths for top ranked variables are colored. The variables
                   (if any) that passed family wide error rate (FWER) cutoff are
                   plotted as solid lines.")
      }
      
      res$stab.feat=stab.feat.c060
    }
  }
  
  if (do.tests) {
    if(do.glmer) {  
      res$glmer = test_taxa_count_glmer(m_a.abs,alpha=alpha,
                                        formula_rhs=formula_rhs,
                                        linfct=linfct)
      report$add.header("Binomial mixed model analysis",3)
      report$add.package.citation("lme4")
      report$add.descr(glmer_descr)
      report.taxa_count_glmer(report,res$glmer)
    }
    if(do.adonis) {
      ##Negative values break bray-curtis and jaccard distances; we standardize to "range" to reduce
      ##the influence of very abundant species:
      count = decostand(m_a$count,method="range",MARGIN=2)
      #or, no standartization:
      #count = m_a$count
      adonis.dist = "bray" #"jaccard" #"bray"
      #ad.res = adonis(count~T1D + Batch,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
      #print(ad.res)
      report$add.header("PermANOVA (adonis) analysis of taxonomic profiles",3)
      report$add.package.citation("vegan")
      report$add.descr("Non-parametric multivariate test for association between
                       taxonomic profiles and meta-data variables. Profile is normalized
                       to proportions across clades and then to range across samples.")
      
      for(adonis.task in adonis.tasks) {
        ad.res = adonis(
          as.formula(paste("count",adonis.task$formula_rhs,sep="~")),
          data=m_a$attr,
          strata=if(!is.null(adonis.task$strata)) m_a$attr[,adonis.task$strata] else NULL,
          permutations=n.adonis.perm,
          method=adonis.dist)
        
        report$add.printed(ad.res,adonis.task$descr)
      }
      
    }
  }
  return (res)
}

## annHeatmap2 from Heatplus, with fixed reordering of labels
annHeatmap2AT <-
  function (x, dendrogram, annotation, cluster, labels, scale = c("row", 
                                                                  "col", "none"), breaks = 256, col = g2r.colors, legend = FALSE) 
  {
    if (!is.matrix(x) | !is.numeric(x)) 
      stop("x must be a numeric matrix")
    nc = ncol(x)
    nr = nrow(x)
    if (nc < 2 | nr < 2) 
      stop("x must have at least two rows/columns")
    def = list(clustfun = hclust, distfun = dist, status = "yes", 
               lwd = 3, dendro = NULL)
    dendrogram = extractArg(dendrogram, def)
    def = list(data = NULL, control = list(), asIs = FALSE, inclRef = TRUE)
    annotation = extractArg(annotation, def)
    def = list(cuth = NULL, grp = NULL, label = NULL, col = RainbowPastel)
    cluster = extractArg(cluster, def)
    def = list(cex = NULL, nrow = 3, side = NULL, labels = NULL, status = "yes")
    labels = extractArg(labels, def)
    if (is.logical(legend)) {
      if (legend) 
        leg = NULL
      else leg = 0
    }
    else {
      if (!(legend %in% 1:4)) 
        stop("invalid value for legend: ", legend)
      else leg = legend
    }
    layout = heatmapLayout(dendrogram, annotation, leg.side = leg)
    x2 = x
    scale = match.arg(scale)
    if (scale == "row") {
      x2 = sweep(x2, 1, rowMeans(x, na.rm = TRUE))
      sd = apply(x2, 1, sd, na.rm = TRUE)
      x2 = sweep(x2, 1, sd, "/")
    }
    else if (scale == "column") {
      x2 = sweep(x2, 2, colMeans(x, na.rm = TRUE))
      sd = apply(x2, 2, sd, na.rm = TRUE)
      x2 = sweep(x2, 2, sd, "/")
    }
    breaks = niceBreaks(range(x2, na.rm = TRUE), breaks)
    col = breakColors(breaks, col)
    dendrogram$Row = within(dendrogram$Row, if (!inherits(dendro, 
                                                          "dendrogram")) {
      dendro = clustfun(distfun(x))
      dendro = reorder(as.dendrogram(dendro), rowMeans(x, na.rm = TRUE))
    })
    dendrogram$Col = within(dendrogram$Col, if (!inherits(dendro, 
                                                          "dendrogram")) {
      dendro = clustfun(distfun(t(x)))
      dendro = reorder(as.dendrogram(dendro), colMeans(x, na.rm = TRUE))
    })
    rowInd = with(dendrogram$Row, if (status != "no") 
      order.dendrogram(dendro)
      else 1:nr)
    colInd = with(dendrogram$Col, if (status != "no") 
      order.dendrogram(dendro)
      else 1:nc)
    x2 = x2[rowInd, colInd]
    labels$Row = within(labels$Row, {
      if (is.null(cex)) 
        cex = 0.2 + 1/log10(nr)
      if (is.null(side)) 
        side = if (is.null(annotation$Row$data)) 
          4
      else 2
      if(status!="no") {
        if (is.null(labels)) 
          labels = rownames(x2)
        else
          ## need to reorder - plot.annHeatmap does not do it
          labels = labels[rowInd]
      }
      else {
        labels = NULL
      }
    })
    labels$Col = within(labels$Col, {
      if (is.null(cex)) 
        cex = 0.2 + 1/log10(nc)
      if (is.null(side)) 
        side = if (is.null(annotation$Col$data)) 
          1
      else 3
      if(status!="no") {
        if (is.null(labels)) 
          labels = colnames(x2)
        else
          ## need to reorder - plot.annHeatmap does not do it
          labels = labels[colInd]
      }
      else {
        labels = NULL
      }
    })
    cluster$Row = within(cluster$Row, if (!is.null(cuth) && (cuth > 
                                                               0)) {
      grp = cutree.dendrogram(dendrogram$Row$dendro, cuth)[rowInd]
    })
    cluster$Col = within(cluster$Col, if (!is.null(cuth) && (cuth > 
                                                               0)) {
      grp = cutree.dendrogram(dendrogram$Col$dendro, cuth)[colInd]
    })
    annotation$Row = within(annotation$Row, {
      data = convAnnData(data, asIs = asIs, inclRef = inclRef)
    })
    annotation$Col = within(annotation$Col, {
      data = convAnnData(data, asIs = asIs, inclRef = inclRef)
    })
    ret = list(data = list(x = x, x2 = x2, rowInd = rowInd, colInd = colInd, 
                           breaks = breaks, col = col), dendrogram = dendrogram, 
               cluster = cluster, annotation = annotation, labels = labels, 
               layout = layout, legend = legend)
    class(ret) = "annHeatmap2AT"
    ret
  }

environment(annHeatmap2AT) <- asNamespace('Heatplus')

plot.annHeatmap2AT <-
  function (x, widths, heights, ...) 
  {
    if (!missing(widths)) 
      x$layout$width = widths
    if (!missing(heights)) 
      x$layout$height = heights
    with(x$layout, layout(plot, width, height, respect = TRUE))
    nc = ncol(x$data$x2)
    nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels)
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) 
      mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) 
      mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
      par(mar = mmar)
      image(1:nc, 1:nr, t(x2), axes = FALSE, xlim = c(0.5, 
                                                      nc + 0.5), ylim = c(0.5, nr + 0.5), xlab = "", ylab = "", 
            col = col, breaks = breaks, ...)
    })
    with(x$labels, {
      if (doRlab) 
        axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, 
             labels = Row$labels, cex.axis = Row$cex)
      if (doClab) 
        axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, 
             labels = Col$labels, cex.axis = Col$cex)
    })
    with(x$dendrogram$Col, if (status == "yes") {
      par(mar = c(0, mmar[2], 3, mmar[4]))
      cutplot.dendrogram(dendro, h = x$cluster$Col$cuth, cluscol = x$cluster$Col$col, 
                         horiz = FALSE, axes = FALSE, xaxs = "i", leaflab = "none", 
                         lwd = x$dendrogram$Col$lwd)
    })
    with(x$dendrogram$Row, if (status == "yes") {
      par(mar = c(mmar[1], 3, mmar[3], 0))
      cutplot.dendrogram(dendro, h = x$cluster$Row$cuth, cluscol = x$cluster$Row$col, 
                         horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", 
                         lwd = x$dendrogram$Row$lwd)
    })
    if (!is.null(x$annotation$Col$data)) {
      par(mar = c(1, mmar[2], 0, mmar[4]), xaxs = "i", yaxs = "i")
      picketPlot(x$annotation$Col$data[x$data$colInd, , drop = FALSE], 
                 grp = x$cluster$Col$grp, grpcol = x$cluster$Col$col, 
                 control = x$annotation$Col$control, asIs = TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
      par(mar = c(mmar[1], 0, mmar[3], 1), xaxs = "i", yaxs = "i")
      picketPlot(x$annotation$Row$data[x$data$rowInd, , drop = FALSE], 
                 grp = x$cluster$Row$grp, grpcol = x$cluster$Row$col, 
                 control = x$annotation$Row$control, asIs = TRUE, 
                 horizontal = FALSE)
    }
    if (x$legend) {
      if (x$layout$legend.side %in% c(1, 3)) {
        par(mar = c(2, mmar[2] + 2, 2, mmar[4] + 2))
      }
      else {
        par(mar = c(mmar[1] + 2, 2, mmar[3] + 2, 2))
      }
      doLegend(x$data$breaks, col = x$data$col, x$layout$legend.side)
    }
    invisible(x)
  }

environment(plot.annHeatmap2AT) <- asNamespace('Heatplus')

## taxa count columns in meta.data must be already sorted by some abundance metrics
heatmap.counts <- function(meta.data,attr.names,attr.annot.names,attr.row.labels,label,
                           max.species.show=30, stand.clust="range", stand.show="range",
                           trans.show=sqrt,
                           attr.order=NULL,
                           agglo.fun.order=sum,
                           cluster.row.cuth=2) {
  
  library(RColorBrewer)
  library(Heatplus)
  library(vegan)
  
  data.norm = row_normalize(meta.data,col_ignore=attr.names)
  
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  
  ##permute samples to make sure that our dendrogram
  ##clustering is not influenced by the original order
  perm.ind = sample(nrow(m_a$count))
  
  count.src = m_a$count[perm.ind,]
  
  if(!is.null(stand.clust)) {
    count = decostand(count.src,method=stand.clust,MARGIN=2)
  }
  else {
    count = count.src
  }
  
  attr = m_a$attr[perm.ind,]
  data.dist.samp <- vegdist(count, method = "bray")
  row.clus <- hclust(data.dist.samp, "ward")
  row.dendro = as.dendrogram(row.clus)
  
  if(is.null(attr.order)) {
    wgts = rowMeans(count, na.rm = TRUE)
  }
  else {
    wgts = attr[,c(attr.order)]
  }
  row.dendro = reorder(row.dendro, wgts, agglo.FUN=agglo.fun.order)
  #row.ind = order.dendrogram(row.dendro)
  
  attr.annot = attr[,attr.annot.names,drop=F]
  
  ## cluster on normalized columns
  count.sub = count[,1:max.species.show]
  # you have to transpose the dataset to get the taxa as rows
  data.dist.taxa <- vegdist(t(count.sub), method = "bray")
  col.clus <- hclust(data.dist.taxa, "ward")
  
  ## go back to un-normalized
  count.sub = count.src[,1:max.species.show]
  if(!is.null(stand.show)) {
    count.sub = decostand(count.sub,method=stand.show,MARGIN=2)
  }
  if(!is.null(trans.show)) {
    count.sub = trans.show(count.sub)
  }
  
  pl.heat = annHeatmap2AT(count.sub,
                          col = colorRampPalette(c("lightyellow", "red"), space = "Lab")(50),
                          #col = heat.colors(50),
                          breaks = niceBreaks(c(min(count.sub),max(count.sub)),50),
                          scale = "none", # we get false bands with default standartization
                          legend = F,
                          dendrogram = list(status="yes",
                                            Row = list(dendro = row.dendro), 
                                            Col = list(dendro = as.dendrogram(col.clus))),
                          labels = list(Col = list(nrow = 20),
                                        Row = if(is.null(attr.row.labels)) list(status="no") 
                                        else list(nrow=6,labels=attr[,attr.row.labels])), #cex=0.8
                          ann = list(Row = list(data = attr.annot)),
                          cluster = list(Row = list(cuth = cluster.row.cuth, 
                                                    col = function(n) {brewer.pal(n, "Set2")})) 
                          # cuth gives the height at which the dedrogram should be cut to form 
                          # clusters, and col specifies the colours for the clusters
  )
  report$add(plot(pl.heat))
}

heatmap.choc <- function(meta.data,attr.names,label) {
  meta.data = cbind(meta.data,Group=meta.data$Sample.type.1)
  heatmap.counts(meta.data,
                 attr.names,
                 attr.annot.names=c("Group","Antibiotic"),
                 attr.row.labels="SubjectID",
                 label=label,
                 stand.clust=NULL)
}

heatmap.t1d <- function(meta.data,attr.names,label) {
  label.ini = label
  meta.data.ini = meta.data
  for(i in 1:6) {
    meta.data = meta.data.ini
    label = label.ini
    if(i>1) {
      meta.data$TimestampMonth = sample(meta.data$TimestampMonth)
      label = paste(label.ini,"randomization ",i,"of TimestampMonth")
    }
  heatmap.counts(meta.data,
                 attr.names,
                 attr.annot.names=c("T1D","TimestampMonth"),
                 attr.row.labels=NULL,
                 label=label,
                 stand.show="max",
                 trans.show=sqrt,
                 attr.order="TimestampMonth",
                 agglo.fun.order=mean)
  }
}


load.meta.t1d.sebastian <- function(file_name,cleared.samp.file=NULL,as.merged=F) {
  
  #Sebastian's code from AnUnivariate.r
  #meta =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  meta =read.csv(file_name, header=TRUE,stringsAsFactors=T)
  meta = meta[!duplicated(meta$SampleID),]
  meta$gender = unlist( lapply( meta$gender, substr, 1,1  ))
  names(meta)[12:14] = c("Family","Person","Reservoir")
  
  meta$BatchRepeats=meta$Batch
  meta$BatchRepeats = ifelse (regexpr("_1$", meta$SampleID)>-1, paste(meta$BatchRepeats, "_1", sep=""), meta$BatchRepeats)
  meta$BatchRepeats = ifelse (regexpr("_2$", meta$SampleID)>-1, paste(meta$BatchRepeats, "_2", sep=""), meta$BatchRepeats)
  meta$BatchRepeats = ifelse (regexpr("_3$", meta$SampleID)>-1, paste(meta$BatchRepeats, "_3", sep=""), meta$BatchRepeats)
  
  meta$T1D[meta$T1D=="Unknown"] = "Control"
  
  x = ddply(meta, .(Family),  function(x) 
  { 
    C= unique(x$gender[x$T1D=="Control"]); 
    T= unique(x$gender[x$T1D=="T1D"]); 
    if (length(C) == 1 & length(T) == 1 ) 
    {
      return (paste(C,T,sep=""))
    } 
    else
    {
      return ("x")
    } 
  } )
  names(x)[2]= "genders"
  meta = merge(meta, x, all.x = TRUE, all.y = TRUE)
  meta$age.quant = quantcut(meta$age)
  row.names(meta) = meta$SampleID
  
  if (as.merged) {
    #Extract unique records for SampleID_merged
    #meta = meta[,.(SampleID_merged,T1D,gender,age,Family,Reservoir,genders,
    #               age.quant,Person,yob,A1C,A1C_bin,DIAGmonthssince,DIAGyearsince)]
    meta = meta[!duplicated(meta$SampleID_merged),]
    meta$SampleID = NULL
    row.names(meta) = meta$SampleID_merged
  }
  
  if(!is.null(cleared.samp.file)) {
    cleaned_yap_input <- read.csv(cleared.samp.file)
    meta$is.cleared = meta$SampleID %in% cleaned_yap_input$SampleID
    #DEBUG:
    if(FALSE) {
      meta = cleaned_yap_input
      meta = meta[!duplicated(meta$SampleID),]
      meta$T1D = meta$status
      meta$Batch = as.integer(substring(meta$Batch,6))
      meta$is.cleared = TRUE
      row.names(meta) = meta$SampleID
      print("Used cleared samples file")
    }
  }
  else {
    meta$is.cleared = TRUE
  }
  return (meta)
}

load.meta.t1d <- function(file_name,batch=NULL,aggr_var=NULL) {
  
  meta =read.csv(file_name, header=TRUE,stringsAsFactors=T)
  
  names(meta) =
    replace.col.names(
      names(meta),
      c("YAP_Aliquot_ID","Family.ID..blinded.", "Subject.ID..blinded.", 
        "Subject.s.Gender", "Autoantibody.Status", "T1D.status",
        "Aliquot_ID", "Subject.s.YEAR.of.birth"),
      c("SampleID",      "FamilyID",            "SubjectID",            
        "Gender",           "AA",                  "T1D",
        "AliquotID",  "YearOfBirth")
    )
  
  meta = meta[!duplicated(meta$SampleID),]
  
  meta$gender = unlist( lapply( meta$gender, substr, 1,1  ))
  
  SampleID_Splits = str_split_fixed(meta$SampleID,"_",3)
  meta$Batch = as.factor(as.numeric(substring(SampleID_Splits[,2],2)))
  ## We set Batch=0 downstream as a special value, it should not be
  ## present already
  stopifnot(sum(meta$Batch==0) == 0)
  meta$BatchRepeat=as.factor(SampleID_Splits[,3])
  
  meta$Specimen.Collection.Date = 
    as.date(as.character(meta$Specimen.Collection.Date),order="mdy")
  meta$age = as.numeric((
    meta$Specimen.Collection.Date
    - 
      as.date(paste(meta$YearOfBirth,"-06-01",sep=""),order="ymd")
  )/365)
  
  meta$T1D[meta$T1D=="Unknown"] = "Control"
  
  meta$Timestamp = 
    strptime(as.character(meta$Timestamp),format = "%m/%d/%Y %H:%M")
  
  meta$BioSampleID = paste(meta$SubjectID,meta$Specimen.Type,as.numeric(meta$Timestamp),sep="_")
  
  make.global(meta)
  
  meta$age.quant = quantcut(meta$age)
  meta$A1C = as.double(as.character(meta$A1C))
  #meta$A1C[is.na(meta$A1C)] = 6.2
  meta$A1C.quant = quantcut(as.double(as.character(meta$A1C)))
  
  meta = arrange(meta,SubjectID,Timestamp)
  
  ## Filtering of rows has to be done before we start marking repeats for throwing
  ## out repeats
  if(!is.null(batch)) {
    report$add.p(paste(c("Filtering metadata by batch:",batch),collapse=" "))
    batch.mask = (meta$Batch %in% batch)
    meta = meta[batch.mask,]
  }
  
  ## That would leave only AA defined samples or T1D patients (under
  ## an assumption that T1D patients are always AA positive even if AA is not measured)
  #meta = meta[meta$AA!="Unknown" | meta$T1D == "T1D",]
  
  ## Subject is at least first repeat, but bio sample is first occurence, in time order 
  ## (e.g. second visit)
  isBioRepeatFirst = duplicated(meta$SubjectID) & ! duplicated(meta$BioSampleID)
  
  BioSampleIDRep = meta$BioSampleID[isBioRepeatFirst]
  
  meta$isBioRepeat = meta$BioSampleID %in% BioSampleIDRep
  
  ## Now do the same as above, but ignore failed samples when looking
  ## for bio repeats
  meta_qa = meta[meta$Sample.QA=="PASS" & meta$Subject.QA == "PASS",]
  
  isBioRepeatFirst = duplicated(meta_qa$SubjectID) & ! duplicated(meta_qa$BioSampleID)
  
  BioSampleIDRep = meta_qa$BioSampleID[isBioRepeatFirst]
  
  ## mark all failed samples or bio repeats following at least one good sample
  meta$isBioRepeatOrFailed = (meta$BioSampleID %in% BioSampleIDRep) | 
    ! (meta$Sample.QA=="PASS" & meta$Subject.QA == "PASS")
  
  ##
  ## This section builds a difference in months between next samples and first sample for every subject
  ##
  meta_keys = meta[,c("SubjectID","SampleID","BioSampleID","Timestamp","isBioRepeat","isBioRepeatOrFailed")]  
  meta_keys_non_rep = meta_keys[!meta$isBioRepeatOrFail,]
  #join() will keep duplicate names, so we rename the names, except SubjectID merge key
  names(meta_keys_non_rep) = paste(names(meta_keys_non_rep),".first",sep="")
  names(meta_keys_non_rep)[1] = "SubjectID"
  meta_keys_1 = join(meta_keys,meta_keys_non_rep,by=c("SubjectID"),match="first",type="left")
  stopifnot(meta$SampleID==meta_keys_1$SampleID)
  stopifnot(is.na(meta_keys_1$Timestamp.first) |
              (meta_keys_1$Timestamp >= meta_keys_1$Timestamp.first) |
              ! (meta$Sample.QA=="PASS" & meta$Subject.QA == "PASS")) #should not see otherwise
  meta_keys_1$MonthsAfterFirstBioSample = as.numeric(meta_keys_1$Timestamp - meta_keys_1$Timestamp.first,units="days")/30
  ## NA will be at records that did not pass QA. Set to 0.
  meta_keys_1$MonthsAfterFirstBioSample[is.na(meta_keys_1$MonthsAfterFirstBioSample)] = 0
  ##
  ## uncheck those bio repeats that passed QA and were made less than 9 mon after the first sampling - they
  ## are testing intentional repeats and should be aggregated just like batch repeats when analyzing 
  ## unpaired group differences
  ##
  meta_keys_1$isAnnualRepeatOrFailed = meta_keys_1$isBioRepeatOrFailed & 
    ! (
      (meta$Sample.QA=="PASS" & meta$Subject.QA == "PASS") & 
        meta_keys_1$MonthsAfterFirstBioSample < 9
    )
  make.global(meta_keys_1)
  meta$isAnnualRepeatOrFailed = meta_keys_1$isAnnualRepeatOrFailed
  meta$MonthsAfterFirstBioSample = meta_keys_1$MonthsAfterFirstBioSample
  meta$TimestampMonth = as.numeric(meta$Timestamp - min(meta$Timestamp),units="days")/30  
  meta$Timestamp = as.numeric(meta$Timestamp)
  meta$TimestampDate = as.Date(timeDate(meta$Timestamp))
  
  row.names(meta) = meta$SampleID
  
  if(!is.null(aggr_var)) {
    
    meta = aggregate_by_meta_data(meta,
                                  aggr_var,
                                  names(meta))  
  }
  
  return (meta)
}


read.t1d <- function(taxa.level=3,batch=NULL,aggr_var_meta = NULL) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("43aa6.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary.txt")
  ##Sebastian's run:
  #moth.taxa <- read.mothur.taxa.summary("ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("acc5a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  ##Andrey's run
  #taxa.file = "dec1a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  ##MiSeq run
  taxa.file = "147936b4aacf1a450a0dfe2d97c3dd5a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  moth.taxa <- read.mothur.taxa.summary(taxa.file)
  ##Vishal's run
  #moth.taxa <- read.mothur.taxa.summary("vishal_all_e61bb.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  #moth.taxa <- read.mothur.taxa.summary("vishal_aa_isdef_3b4ea.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  make.global(taxa.lev.all)
  #taxa.lev = count_filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
  taxa.lev = taxa.lev.all
  make.global(taxa.lev)
  meta = load.meta.t1d("seq_id_to_metadata.csv",batch=batch,aggr_var = aggr_var_meta)
  make.global(meta)
  #meta = load.meta.t1d("annotation20130819.csv")
  #meta = load.meta.t1d("annotation20130819.csv","all_batches_cleaned_yap_input.csv")
  report$add.p(sprintf("Loaded %i records from taxonomic assignment file %s for taxonomic level %i",
                       nrow(taxa.lev),taxa.file,taxa.level))
  taxa.meta = merge.counts.with.meta(taxa.lev,meta)
  report$add.p(sprintf("After merging with metadata, %i records left",
                       nrow(taxa.meta)))
  return(taxa.meta)
}

plot.features.mds <- function(m_a,species.sel,sample.group) {
  ##https://stat.ethz.ch/pipermail/r-help/2008-April/159351.html
  ##http://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
  m = m_a$count
  ##default distance is bray
  m = decostand(m,method="range",MARGIN=2)
  sol = metaMDS(m,autotransform = F,trymax=40)
  site.sc <- scores(sol, display = "sites")
  species.sc <- scores(sol, display = "species")
  plot(site.sc)
  points(sol,display="sites",col=sample.group)
  points(sol,display="species",pch="x",select=species.sel)
  text(sol,display="species", cex=0.7, col="blue",select=species.sel)
  #points(site.sc,col=m_a$attr$T1D)
  #points(species.sc)
  #points(species.sc["Streptococcus_0.1.11.1.2.6.2",1:2,drop=F],pch="+")
  legend(-0.5,1,unique(sample.group),col=1:length(sample.group),pch=1)
}

proc.t1d.som <- function() {
  library(kononen)
  m.som = som(m, grid = somgrid(5, 5, "hexagonal"))
  m.pred = predict(m.som,newdata=m,trainX=m,trainY=m_a$attr$T1D)
  table(m_a$attr$T1D,m.pred$prediction)
  #plot assignment of original data
  plot(m.som,type="mapping",col=as.numeric(m_a$attr$T1D=="T1D")+1,pch=1)
  plot(m.som,type="mapping",col=as.numeric(m_a$attr$Batch),pch=2)
}


## Build mixed effects model for count data of a single clade
## as a function of patient-control classification
## (Binomial regression model with a subject specific random effect)
## Value: p value of testing against a null hypothesis that 
## model coefficients are zero (ignoring intercept; two-sided)
test_taxa_count_glmer_col <- function(taxa.count,
                                      attr,
                                      taxa.count_rowsum,
                                      formula_rhs,
                                      linfct) {
  
  
  ## T1D groups samples into patient and control groups (fixed effect).
  ## FamilyID groups patients and their control siblings into 
  ## matched groups (random effect).
  ## We add the unique ID of each observation (SampleID) as a random effect, in order 
  ## to account for overdispersion.
  ## This approach has been used in a T1D study [10.1371/journal.pone.0025792]:
  ## (http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0025792)
  ## Specific implementations are described in
  ## (http://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html)
  ## and
  ## (http://depts.washington.edu/cshrb/wordpress/wp-content/uploads/2013/04/Tutorials-Tutorial-on-Count-Regression-R-code.txt)
  ## 
  ##ANOVA comparison for models with and without the overdispersion term in the
  ## case of Streptococcus is below:
  ##> anova(dat.glm.0,dat.glm)
  ##Data: dat
  ##Models:
  ##  dat.glm.0: y ~ T1D + (1 | FamilyID)
  ##dat.glm: y ~ T1D + (1 | FamilyID) + (1 | SampleID)
  ##Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
  ##dat.glm.0  3 5837.3 5846.1 -2915.66   5831.3                             
  ##dat.glm    4 1287.0 1298.7  -639.49   1279.0 4552.4      1  < 2.2e-16 ***
  ##  ---
  ##  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1  
  response = cbind(taxa.count,taxa.count_rowsum-taxa.count)
  test_glmer_inner <- function() {
    
    ##default optimizer was not converging often, 
    ##the developer recommended using bobyga:
    ##http://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
    
    dat.glmer = glmer(as.formula(paste("response",formula_rhs,sep="~")),
                      data=attr,family="binomial",
                      control=glmerControl(optimizer="bobyqa") 
                      ## need to load package "optimx" for the optimizer below
                      #control=glmerControl(optimizer="optimx",
                      #                    optCtrl=list(method="L-BFGS-B") #"nlminb"
                      #)
    )
    ## For some reason, T1D name soemtimes becomes T1DT1D in the output - it looks like
    ## there is a bug in glmer, and glmer corrupts something in its internal state,
    ## and generates dummy indicator variables for the factors with names build from
    ## contrasts values rather than the original factor levels, when the whole
    ## process is executed again in the same R session. Run the whole analysis each time
    ## in a fresh R session.
    
    #p_val = summary(dat.glmer)$coefficients[,"Pr(>|z|)"]["T1DT1D"]
    ## If we have multiple fixed effects factors, then the right way to obtain
    ##p-vals is below, using multiple testing correction
    #library(multcomp)
    dat.glht = glht(dat.glmer,linfct=linfct)
    dat.glht.summ = summary(dat.glht)
    p_val = as.numeric(dat.glht.summ$test$pvalues)
    p_val
  }
  
  #warn_saved = options()$warn
  #options(warn=2)
  ok=T
  #somehow if both warning and error arguments are used
  #in a single call to tryCatch, then error handler is not
  #called
  p_val = tryCatch(tryCatch(test_glmer_inner(),
                            
                            warning=function(w) {
                              warning(paste("Warning caught in glmer; results converted to NA: ",w,"\n"))
                              ok<<-F
                            }
  ),                   error=function(e) {
    warning(paste("Error caught in gmler; results converted to NA: ",e,"\n"))
    ok<<-F
  }
  )
  
  #options(warn=warn_saved)
  
  if (!ok) {
    
    print("Warnings or errors in glmer, returning NA, making data available for debugging")
    
    #cat(xtabs(~(response[,"taxa.count"]>0)+attr$Batch+attr$T1D))
    
    #make.global(dat.glmer)
    make.global(attr)
    make.global(response)
    make.global(formula_rhs)
    
    p_val = NA
    
  }
  #print(paste("p_val=",p_val))
  p_val
}

## Test all clades pairwise using Poisson regression model with 
## a subject specific random effect
test_taxa_count_glmer <- function(m_a,
                                  formula_rhs,
                                  linfct,
                                  alpha=0.05,
                                  p_adjust_method="BH") {
  ##http://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html
  library(lme4)
  library(multcomp)
  cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
  registerDoSNOW(cl)
  count = m_a$count
  if(F) {
    count = m_a$count[,c(1:7)]
  }
  if(F) {
    count = m_a$count[,c("Bifidobacterium_0.1.1.1.2.1.2",
                         "Asaccharobacter_0.1.1.1.3.1.2",
                         "Escherichia_Shigella_0.1.5.5.1.1.1",
                         "Eggerthella_0.1.1.1.3.1.6",
                         "Streptococcus_0.1.11.1.2.6.2")]
  }
  count_rowsum = rowSums(count)
  #attr = m_a$attr[,c("T1D","FamilyID","SampleID","Batch")]
  attr = m_a$attr
  # Not working properly with .parallel=T in the presence
  # of warnings or errors - getting all NAs no matter what
  # I try for error handling in tryCatch
  p_vals = aaply(count,
                 2,
                 test_taxa_count_glmer_col,
                 attr,
                 count_rowsum,
                 formula_rhs=formula_rhs,
                 linfct=linfct,
                 .inform=F,
                 .parallel=F,
                 .paropts=list(.packages=c("lme4","multcomp")))
  stopCluster(cl)
  names(p_vals) = colnames(count)
  p_vals = p_vals[order(p_vals,na.last=T)]
  p_vals_signif = p_vals[!is.na(p_vals) & p_vals<=alpha]
  p_vals_adj = p.adjust(p_vals,method=p_adjust_method)
  p_vals_adj_signif = p_vals_adj[!is.na(p_vals_adj) & p_vals_adj<=alpha]
  ret = llist(p_vals,
              p_vals_adj,
              p_vals_signif,
              p_vals_adj_signif,
              formula_rhs,
              linfct,
              alpha,
              p_adjust_method)
  class(ret) = "taxa_count_glmer"
  
  return(ret)
}

print.taxa_count_glmer <- function(x,...) {
  cat ("Test results from fitting mixed effects binomial model for each taxa\n")
  cat ("Right-hand side of formula:\n")
  cat (x$formula_rhs,"\n")
  cat ("Specification of the linear hypotheses (see glht):\n")
  cat (x$linfct,"\n")
  cat ("p-values that pass significance cutoff before multiple testing correction:\n")
  print(x$p_vals_signif)
  cat ("p-values that pass significance cutoff after BH multiple testing correction:\n")
  print(x$p_vals_adj_signif)
}

report.taxa_count_glmer <- function(report,x,...) {
  report$add.p("Test results from fitting mixed effects binomial model for each taxa")
  report$add.p("Right-hand side of formula:")
  report$add.p(x$formula_rhs)
  report$add.p("Specification of the linear hypotheses (see glht):")
  report$add.p(x$linfct)
  
  report$add.vector(x$p_vals_signif,"p_value",
                    caption="p-values that pass significance cutoff before multiple testing correction.")
  report$add.vector(x$p_vals_adj_signif,"p_value",
                    caption="p-values that pass significance cutoff after BH multiple testing correction.")
  report$add.vector(x$p_vals_adj[1:min(10,length(x$p_vals_adj))],"p_value",
                    caption="Top 10 p-values after BH multiple testing correction.")
  failed.names = names(x$p_vals[is.na(x$p_vals)])
  if(length(failed.names)>0) {
    report$add(list(failed.names),
               caption="Clades for which the model could not be built (e.g. too many all-zeros samples)")
  }
}

## It is assumed that the count matrix x is transformed/normalized already
## e.g. with decostand(ihs(count),method="standardize",MARGIN=2),
## otherwise set tran_norm to TRUE and the command above will be
## applied
stab_sel_genesel <- function(x,
                             y,
                             samp_in_col=F,
                             tran_norm=F,
                             type="unpaired",
                             replicates=400,
                             fold_ratio=0.5,
                             maxrank=10,
                             samp_filter=NULL) {
  library(GeneSelector)
  if(!samp_in_col) {
    x = t(x)
  }
  if(!is.null(samp_filter)) {
    x = x[,samp_filter]
    y = y[samp_filter]
  }
  if(tran_norm) {
    x = decostand(ihs(x),method="standardize",MARGIN=1)
  }
  
  ranking_methods = c(RankingWilcoxon,RankingLimma,RankingFC)
  ranking_method = RankingWilcoxon
  n_feat = nrow(x)
  n_samp = ncol(x)
  #contrary to the help page, does not work when y is a factor - need to convert
  fold_matr = GenerateFoldMatrix(y = as.numeric(y), k=trunc(n_samp*(1-fold_ratio)),replicates=replicates)
  rnk = ranking_method(x,y,type=type,pvalues=T)
  pvals.adjusted = p.adjust(rnk@pval,method="BH")
  #make.global(rnk)
  #make.global(pvals.adjusted)
  rep_rnk = RepeatRanking(rnk,fold_matr,scheme="subsampling",pvalues=T)
  #make.global(rep_rnk)
  #toplist(rep_rnk,show=F)
  stab_ovr = GetStabilityOverlap(rep_rnk, scheme = "original", decay = "linear")
  ### for a short summary
  #summary(stab_ovr, measure = "intersection", display = "all", position = 10)
  #summary(stab_ovr, measure = "overlapscore", display = "all", position = 10)
  ### for a graphical display
  #plot(stab_ovr)
  #aggr_rnk = AggregateSimple(rep_rnk, measure="mode")
  aggr_rnk = AggregateMC(rep_rnk, maxrank=n_feat)
  #make.global(aggr_rnk)
  #toplist(aggr_rnk)
  #gsel = GeneSelector(list(aggr_rnk), threshold = "BH", maxpval=0.05)
  gsel = GeneSelector(list(aggr_rnk), threshold = "user", maxrank=maxrank)
  #show(gsel)
  #str(gsel)
  #toplist(gsel)
  selected = SelectedGenes(gsel)
  #pvals are always NA somehow, remove the field
  selected$pvals=NULL
  selected = cbind(name=rownames(x)[selected$index],
                   selected,
                   pvals.orig=rnk@pval[selected$index],
                   pvals.adjusted=pvals.adjusted[selected$index])
  #print(selected)
  #make.global(gsel)
  return (list(stab_feat=selected,gsel=gsel))
}

feat_sel_samr <- function(m_a.abs) {
  library(samr) #Tibshirani's package for feature selection in microarrays and RNASeq
  ##TODO: study assumptions of this method on sequence counts. The help page mentions only RNASeq.
  ##Something is probably not right because it only reports "genes down" and empty for "genes up" for
  ##T1D genus count data
  samfit = SAMseq(t(m_a.abs$count),m_a.abs$attr$T1D,resp.type="Two class unpaired",geneid=colnames(m_a.abs$count))
  print(samfit)
  plot(samfit)
}

## Taken from (http://depts.washington.edu/cshrb/wordpress/wp-content/uploads/2013/04/Tutorials-Tutorial-on-Count-Regression-R-code.txt)
### Small utility function to get (conditional) rate ratios and
### 95% CI from fitted glmer() object -- only appropriate for 
### a Poisson model (though, we might do similar things with a
### binomial outcome)
lmerCI <- function(obj, rnd = 2) {
  cc <- fixef(obj)
  se <- sqrt(diag(vcov(obj)))
  upper <- cc + 1.96*se
  lower <- cc - 1.96*se
  out <- data.frame(RR = round(exp(cc), rnd), 
                    upper = round(exp(upper), rnd), 
                    lower = round(exp(lower), rnd))
  out
}


proc.t1d <- function() {
  #taxa.levels = c(2,3,4,5,6)
  taxa.levels = c(6)
  #batches = list(c(1,2,3),c(1),c(2),c(3),c(2,3),c(1,3))
  batches = list(c(1,2,3))
  do.std.plots = T
  do.tests = T
  
  report$add.descr("Largely identical set of analysis routines is applied
                   in nested loops over combinations of batches and 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out the batch combination and taxonomic level.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture.")
  
  for (batch in batches) {
    
    report$add.header(paste(c("Batch combination:",batch),collapse=" "),1)
    label_batch = paste("16s","b",paste(batch,collapse="-"),sep=".")
    report$set.tag(label_batch)
    
    for (taxa.level in taxa.levels) {
      taxa.meta = read.t1d(taxa.level,batch=batch,aggr_var_meta = "AliquotID")
      report$add.header(paste("Taxonomic level:",taxa.level),2)
      report$add.p(paste("Before filtering for QAed and redundant samples:",nrow(taxa.meta$data)))
      
      label = paste(label_batch,"l",taxa.level,sep=".",collapse=".")
      report$add.p(paste("Tag:",label))
      report$set.tag(label)
      
      aggrBySubject = F
      if(!aggrBySubject) {
        taxa.meta$data = taxa.meta$data[taxa.meta$data$Sample.QA=="PASS" & 
                                          taxa.meta$data$Subject.QA=="PASS",]
        taxa.meta.aggr = taxa.meta
        report$add.p(paste("After filtering for QAed samples:",nrow(taxa.meta$data)))      
      }
      else {
        taxa.meta$data = taxa.meta$data[!taxa.meta$data$isAnnualRepeatOrFailed,]
        report$add.p(paste("After filtering for QAed and non-annual repeat samples:",nrow(taxa.meta$data)))
        
        #aggr_var = "AliquotID"
        aggr_var = "SubjectID"
        
        taxa.meta.aggr = aggregate_by_meta_data(taxa.meta$data,
                                                aggr_var,
                                                taxa.meta$attr.names)
        ##When we aggregate across SubjectID, Batch is no longer valid
        ##we reset it here to switch off batch effect analysis downstream
        taxa.meta.aggr$data$Batch = as.factor(0)
        report$add.p(paste("After aggregating samples by ", aggr_var, ":",nrow(taxa.meta.aggr$data)))
      }
      ## no batches in MiSeq run
      taxa.meta.aggr$data$Batch = as.factor(0)
      n.batches = sum(table(taxa.meta.aggr$data$Batch)>0)
      
      if(n.batches>1) {
        xtabs.formula = ~T1D+Batch
      }
      else {
        xtabs.formula = ~T1D
      }
      fact.xtabs = xtabs(~T1D+Batch,data=taxa.meta.aggr$data,drop.unused.levels=T)
      report$add.table(fact.xtabs,show.row.names=T,caption="Sample cross tabulation")
      report$add.printed(summary(fact.xtabs))
      
      taxa.meta.aggr$data = count_filter(taxa.meta.aggr$data,
                                         col_ignore=taxa.meta.aggr$attr.names,
                                         min_max=10)
      #taxa.meta.data = count_filter(taxa.meta.data,col_ignore=taxa.meta.attr.names,
      #                              min_median=0.002,
      #                              min_max_frac=0.1,min_max=10,
      #                              min_row_sum=500,other_cnt="other")
      
      report$add.p(paste("After count filtering,",
                         (ncol(taxa.meta.aggr$data)-length(taxa.meta.aggr$attr.names)),"clades and",
                         nrow(taxa.meta.aggr$data),"rows left."))
      
      
      
      with(taxa.meta.aggr$data,{
        print(levels(T1D))
        print(contrasts(T1D))
      }
      )
      
      if (do.tests) {
        res.tests = try(
          test.counts.t1d(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,
                          label=label,
                          stability.transform.counts="ihs",
                          do.stability=F,
                          do.tests=T,
                          do.genesel=F,do.glmnet=F,
                          do.glmer=T,do.adonis=F,
          )
        )
      }
      if (do.std.plots) {
        tryCatchAndWarn({
          heatmap.t1d(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,label=label)
        })
        
        try({
          plot.group = list(c("T1D"))
          if(n.batches>1) {
            plot.group[[2]] = c("T1D","Batch")
          }
          std.plots(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,id.vars.list=
                      plot.group,
                    label=label,
                    res.tests=res.tests
          )
        })
      }
    }
  }
}

read.t1d.mg <-function(annot.type,level) {
  annot.dir = "../BATCH_01_02_META"
  mgrast.dir = paste(annot.dir,"BATCH_01-02_METAGENOMICS_MGRAST",sep="/")
  if (annot.type == "humann") {
    counts = read.humann.summary(paste(annot.dir,"humann/output/04b-keg-mpt-cop-nul-nve-nve.txt",sep="/"))
    counts = count_filter(counts,col_ignore=c(),min_max_frac=0,min_max=0,min_row_sum=0,other_cnt="other")      
  }
  else if (annot.type %in% c("cog","kegg","subsys")) {
    make.global(level)
    counts = read.mgrast.summary(paste(mgrast.dir,paste(annot.type,"tsv",sep="."),sep="/"),
                                 file_name.id.map=paste(mgrast.dir,"mgrast_to_samp_id.tsv",sep="/"))
    make.global(counts)
    print(level)
    counts = mgrast.to.abund.df(counts,level)
    counts = count_filter(counts,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=0,other_cnt="other")
  }
  #make.global(counts)
  
  meta = load.meta.t1d("annotation20130819.csv",as.merged=T)
  return (merge.counts.with.meta(counts,meta))
}

proc.t1d.mg <- function() {
  annot.types = c("humann","cog","subsys","kegg")
  #annot.types = c("subsys")
  do.std.plots = T
  do.tests = T
  for (annot.type in annot.types) {
    
    levels = annot.levels(annot.type)
    if(annot.type=="humann") {
      do.diversity = F
    }
    else {
      do.diversity = T
    }
    do.diversity = F
    
    for (level in levels) {
      
      label = paste(annot.type,level,sep=".")
      print (paste("Working on",label))
      
      count.meta = read.t1d.mg(annot.type,level)
      count.meta.data = count.meta$data
      count.meta.attr.names = count.meta$attr.names
      make.global(count.meta.data)
      
      if (do.tests) {
        
        res.tests = try(test.counts.t1d(count.meta.data,count.meta.attr.names,
                                        label=label,
                                        stability.transform.counts="ihs",
                                        do.stability=T,
                                        do.tests=T))
      }
      if (do.std.plots) {
        try(
          std.plots(count.meta.data,count.meta.attr.names,id.vars.list=
                      list(
                        c("T1D","age.quant"),
                        #c("Family","T1D"),
                        #c("SampleID","Batch"),
                        c("T1D"),
                        c("T1D","Batch")
                      ),
                    label=label,
                    do.diversity=do.diversity,
                    res.tests=res.tests
          )
        )
      }
    }
  }
}


test.counts.mr_oralc <- function(data,attr.names,label,alpha=0.05,
                                 do.tests=T,do.stability=T,
                                 stability.transform.counts="ihs") {
  
  n.adonis.perm = 4000
  res = list()
  #data = data[data$Batch %in% c(1,2,3),]
  
  ##only families with 1 sibling and 1 patient
  #tbl = with(data,table(Family,T1D))
  #fams = unique(row.names(tbl[(tbl[,"T1D"]>0 & tbl[,"Control"]>0),]))
  #data = data[data$Family %in% fams,]
  
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  
  data.norm = row_normalize(data,col_ignore=attr.names)
  
  #make.global(tbl)
  #make.global(fams)
  #make.global(data.norm)
  
  
  all.clades = get.clade.names(data,attr.names)
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  make.global(m_a)
  
  m = m_a$count
  
  #m = (m_a.abs$count > 0)
  #storage.mode(m) = "integer"
  
  count = switch(stability.transform.counts,
                 boxcox=boxcox.transform.mat(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m_a.abs$count > 0))
  make.global(count)
  
  if (do.stability) {
    
    standardize.glm = T
    
    cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
    registerDoSNOW(cl)  
    print(levels(m_a$attr$sample.type))
    cv.res = cv.glmnet.alpha(m_a$attr$sample.type,count,family="binomial",standardize=standardize.glm)
    stopCluster(cl)
    
    make.global(cv.res)
    penalty.alpha = cv.res$alpha
    #alpha = 0.8
    stab.res.c060 = stabpath(m_a$attr$sample.type,count,weakness=0.9,
                             family="binomial",steps=600,
                             alpha=penalty.alpha,standardize=standardize.glm)
    make.global(stab.res.c060)
    fwer = alpha
    pi_thr = 0.6
    stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    make.global(stab.feat.c060)
    print ("Features that meet FWER cutoff in stability analysis:")
    print (names(stab.feat.c060$stable))
    
    stab.path.file = paste("stability.path",label,"png",sep=".")
    
    #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
    p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
    ggsave(stab.path.file)
    #png(stab.path.file)
    #plot(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    #dev.off()
    res$stab.feat=stab.feat.c060
  }
  
  if (do.tests) {
    ##Negative values break bray-curtis and jaccarda= distances; we standardize to "range" to reduce
    ##the influence of very abundant species
    count = decostand(m_a$count,method="range",MARGIN=2)
    adonis.dist = "jaccard" #"bray"
    #ad.res = adonis(count~sample.type + Batch,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    #print(ad.res)
    ad.res.sample.type = adonis(count~sample.type,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    print(ad.res.sample.type)
    ad.res.paired = adonis(count~sample.type,data=m_a$attr,strata=m_a$attr$Subject,permutations=n.adonis.perm,method=adonis.dist)  
    print (ad.res.paired)
    #make.global(ad.res)
    
    #print (ad.res)
    #test.ad.res = ad.res$aov.tab$"Pr(>F)"[1]
    #r2.ad = ad.res$aov.tab$R2[1]
    
    gr.abs = m_a.abs$attr$sample.type
    gr.abs.cnt = list(x=m_a.abs$count[gr.abs=="Tumor Tissue",],y=m_a.abs$count[gr.abs=="Healthy Tissue",])
    print(Xdc.sevsample(gr.abs.cnt))
    K = dim(gr.abs.cnt[[1]])[2]-1
    print (K)
    #make.global(gr.abs.cnt)
    print(Xmcupo.sevsample(gr.abs.cnt,K))
    
    ##this is a rank based pairwise test; column-wise standartization does not
    ##affect it, but row-wise - does
    gr = m_a$attr$sample.type
    gr.cnt = list(x=m_a$count[gr=="Tumor Tissue",],y=m_a$count[gr=="Healthy Tissue",])
    
    test.wil.res = list()
    
    for (clade in all.clades) {
      test.wil.res[[clade]] = wilcox.test(gr.cnt[[1]][,clade],gr.cnt[[2]][,clade],paired=FALSE,exact=F)$p.value
    }
    test.wil.res = unlist(test.wil.res)
    make.global(test.wil.res)    
    names(test.wil.res) = all.clades
    test.wil.res.adj = p.adjust(test.wil.res,method="BH")
    #does not make sense on so few samples, exits with errors:
    if(F) {
      fdr.res = fdrtool(test.wil.res,"pvalue")
      test.wil.res.adj.q = fdr.res$qval
      names(test.wil.res.adj.q) = names(test.wil.res)
      print(test.wil.res.adj.q[test.wil.res.adj.q<=alpha])
      make.global(test.wil.res.adj.q)
    }
    names(test.wil.res.adj) = names(test.wil.res)
    print("Significant unadjusted p-values:")
    print(test.wil.res[test.wil.res<=alpha])
    print("Significant BH adjusted p-values:")    
    print(test.wil.res.adj[test.wil.res.adj<=alpha])
    make.global(test.wil.res.adj)
    
  }
  return (res)
}


load.meta.mr_oralc <- function(file_name) {
  
  #Sebastian's code from AnUnivariate.r
  #meta =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  meta =read.delim(file_name, header=TRUE,stringsAsFactors=T)
  row.names(meta) = meta$ID
  return (meta)
}


read.mr_oralc <- function(taxa.level=3) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("43aa6.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary.txt")
  moth.taxa <- read.mothur.taxa.summary("5e2c2.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count_filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
  #taxa.lev = taxa.lev.all
  meta = load.meta.mr_oralc("UCFTumorTissueSampleInventory.cleaned.txt")
  return (merge.counts.with.meta(taxa.lev,meta))
}


proc.mr_oralc <- function() {
  taxa.levels = c(2,3,4,5,6)
  #taxa.levels = c(6)
  do.std.plots = T
  do.tests = T
  for (taxa.level in taxa.levels) {
    taxa.meta = read.mr_oralc(taxa.level)
    make.global(taxa.meta)
    taxa.meta.data = taxa.meta$data[taxa.meta$data$sample.type %in% c("Tumor Tissue","Healthy Tissue"),]
    taxa.meta.data$sample.type = as.factor(as.character(taxa.meta.data$sample.type))
    print(levels(taxa.meta.data$sample.type))
    print(sort(unique(as.character(taxa.meta.data$sample.type))))
    make.global(taxa.meta.data)    
    taxa.meta.attr.names = taxa.meta$attr.names
    label = paste("16s",taxa.level,sep=".")
    
    print (paste("Working on",label))
    
    if (do.tests) {
      res.tests = try(
        test.counts.mr_oralc(taxa.meta.data,taxa.meta.attr.names,
                             label=label,
                             stability.transform.counts="ihs",
                             do.stability=T,
                             do.tests=T)
      )
    }
    if (do.std.plots) {
      try(
        std.plots(taxa.meta.data,taxa.meta.attr.names,id.vars.list=
                    list(
                      c("Subject","sample.type"),
                      c("sample.type"),
                      c("Subject")
                    ),
                  label=label,
                  res.tests=res.tests
        )
      )
    }
  }
}

read.pieper.t1d <- function() {
  
  file_name = "aim3/September 28 analysis_T1D.txt"
  #data =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  data =read.delim(file_name, header=F,sep="\t",stringsAsFactors=FALSE)
  #row.names(meta) = meta$ID
  group = c(t(data[1,][-(1:4)]))
  #rownames(group) = NULL
  make.global(group)
  data = t(data[-1,-c(1,(3:4))])
  col_names = data[1,]
  #col_names = sub("UniRef.*\\_([0-9A-Z]{4,})","\\1",col_names)
  data = data[,!duplicated(col_names)]
  colnames(data) = data[1,]
  rownames(data) = data[,1]
  colnames(data)[1] = "ID"
  data = as.data.frame(data[-1,])
  data$ID = as.factor(data$ID)
  data[,-1] = as.numeric(as.matrix(data[,-1]))
  data$group = group
  #data = cbind(data,group)
  data = data[data$group != 0,]
  data$group = as.factor(data$group)  
  data$id.person = as.factor(sub(".* ([A-Z]{5,10}).*","\\1",data$ID))
  #data = ddply(data,.(id.person),colwise(function(y) {if(is.numeric(y)) sum(y) else y[1]} ))
  data$ID = NULL
  if(F) {
    id.group = data[,c("id.person","group")]
    #id.group = id.group[!duplicated(id.group$id.person),]
    #rownames(id.group)=id.group$id.person
    #id.group$id.person = NULL
    #id.group = id.group[order(rownames(id.group)),]
    data$group = NULL
    data = ddply(data,.(id.person),function(x) {id.person=x$id.person[1]
                                                x$id.person=NULL
                                                y = colSums(x)
                                                c(id.person,y)})
    rownames(data) = data[,1]
    make.global(data)
    make.global(id.group)
    #data = data[order(rownames(data)),]
    #stopifnot(all(rownames(id.group)==rownames(data)))
    #data$group = id.group$group
    nrow.before = nrow(data)
    data = join(data,id.group,by=c("id.person"),type="left",match="first")
    stopifnot(!any(aaply(c(data$group),1,is.null)))
    stopifnot(nrow(data)==nrow.before)
  }
  data = data[!duplicated(data$id.person),]
  return (list(data=data,attr.names=c("id.person","group")))
}

power.pieper.t1d <- function() {
  alpha = 0.05
  use.fdrtool = T
  taxa.meta = read.pieper.t1d()
  make.global(taxa.meta)
  taxa.meta.attr.names = taxa.meta$attr.names    
  taxa.meta.data.raw = count_filter(taxa.meta$data,col_ignore=taxa.meta.attr.names,
                                    min_max_frac=0,min_max=0,min_median=100,min_row_sum=0,other_cnt=NULL)
  taxa.meta.data = all_normalize(taxa.meta.data.raw,col_ignore=taxa.meta.attr.names,norm.func=ihs)
  make.global(taxa.meta.data)
  clade.names = get.clade.names(taxa.meta.data,taxa.meta.attr.names)
  #print(clade.names)
  #pvals = wilcox.test.multi(data=taxa.meta.data,resp.vars=clade.names,group.var="group",subset=NULL)
  m = count_matr_from_df(taxa.meta.data,taxa.meta.attr.names)
  #mtp.res = MTP(X=t(m),Y=taxa.meta.data$group,get.adjp=F)
  tr.m = t(m)
  group = taxa.meta.data$group
  make.global(m)
  make.global(group)
  m.raw = count_matr_from_df(taxa.meta.data.raw,taxa.meta.attr.names)
  pvals = wilcox.test.multi.fast(tr.m,group)  
  make.global(pvals)
  if(use.fdrtool) {
    pvals.adj = fdrtool(pvals,statistic="pvalue",plot=F,verbose=F)$qval
  }
  else {
    pvals.adj = p.adjust(pvals,method="BH")
  }
  make.global(pvals.adj)
  ind.sig = which(pvals.adj<=alpha)
  group.mean = ddply(cbind(as.data.frame(m),group),.(group),colwise(mean))
  group.sd = ddply(cbind(as.data.frame(m),group),.(group),colwise(sd))
  make.global(ind.sig)
  group.mean.sig = group.mean[,ind.sig]
  #print(group.mean.sig)
  group.sd.sig = group.sd[,ind.sig]
  #print(group.sd.sig)
  group.n = c(table(group))[rownames(group.mean.sig)]
  cohens.d = foreach(i.var=seq(ncol(group.mean.sig)),.combine=c) %do% {
    cohens.d.from.mom(mean.gr=group.mean.sig[,i.var],var.gr=group.sd.sig[,i.var]**2,n.gr=group.n)
  }
  #print(cohens.d)
  pvals.boot = boot(data=m,
                    statistic=booted.wilcox.test.multi.fast,
                    R=4000,
                    n=50,
                    strata=group,
                    group=group,
                    test.method=wilcox.test.multi.fast)$t
  make.global(pvals.boot)
  if(use.fdrtool) {
    pvals.boot.adj = aaply(pvals.boot,1,function(x) fdrtool(x,statistic="pvalue",plot=F,verbose=F)$qval)
  }
  else {
    pvals.boot.adj = aaply(pvals.boot,1,p.adjust,method="BH")
  }
  make.global(pvals.boot.adj)
  power.sig = colMeans(pvals.boot.adj[,ind.sig] <= alpha)
  make.global(power.sig)
  #print(power.sig)
  #print(mean(power.sig))
  return(list(cohens.d=cohens.d,
              mean.cohens.d=mean(cohens.d),
              group.mean.sig=group.mean.sig,
              group.sd.sig=group.sd.sig,group.n=group.n,
              power.sig=power.sig,
              mean.power.sig=mean(power.sig),
              pvals.adj.sig=pvals.adj[ind.sig]))
}

pediatric.cancer.2013.aim2 <- function() {
  doc = "(http://www.ncbi.nlm.nih.gov/pubmed/17414136)
  
  Fecal specimens were collected from 148 consecutive pediatric patients 
  (79 with Crohn disease, 62 with ulcerative colitis, and 7 with irritable 
  bowel syndrome) and 22 healthy control individuals.
  
  Lactoferrin levels were significantly higher in patients with ulcerative 
  colitis (1880 +/- 565 microg/mL) (mean +/- SE) or Crohn disease 
  (1701 +/- 382 microg/mL) than in healthy control individuals under 21 
  years of age (1.17 +/- 0.47 microg/mL, P < 0.001). 
  "
  
  mean.gr = c(1880,1.17)
  se.gr = c(565,0.47)
  n.gr = c(79,22)
  
  var.gr = se.gr**2*n.gr
  sd.gr = var.gr**0.5
  
  mom.f.gr = foreach(i.gr=1:2,.combine=rbind) %do% {
    mom.f(var.x=var.gr[i.gr],mean.x=mean.gr[i.gr],f=ihs,f.d1=ihs.d1,f.d2=ihs.d2)
  }
  cohens.d = cohens.d.from.mom(mean.gr=mom.f.gr[,"mean.f"],var.gr=mom.f.gr[,"var.f"],n.gr=n.gr)
  return(list(mom.f.gr=mom.f.gr,mean.gr=mean.gr,
              se.gr=se.gr,var.gr=var.gr,sd.gr=sd.gr,
              cohens.d=cohens.d))
}

power.pediatric.cancer.2013<-function() {
  #proc.choc()
  power.res = power.pieper.t1d()
  print("Power results Aim 3:")
  print(power.res)
  mom = pediatric.cancer.2013.aim2()
  print("Power results Aim 2:")
  print(paste("Cohen's d Lactoferrin:",mom$cohens.d))
  print(paste("Mean Cohen's d for significant proteins from Aim 3:",power.res$mean.cohens.d))
}

