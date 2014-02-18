#install.packages("vegan")
#install.packages("HMP")
#install.packages("reshape2")
#install.packages("ggplot2")
#install.packages("gtools")
#install.packages("BiodiversityR")
#install.packages("LiblineaR")
#install.packages("BatchJobs")

##two alternative implementations of the same stability
##selection paper, and different classification methods
#GLMs
#install.packages("c060")
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
#source("http://bioconductor.org/biocLite.R")
#biocLite("multtest")
#biocLite("GeneSelector")

if(F) {
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
  "c060", 
  "quadrupen", 
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
  "knitr",
  "lattice"
)

for (package in packages) {
  library(package,character.only=T)
}


count_matr_from_df<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  as.matrix(dat[!mask_col_ignore])
}

split_count_df<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  m = as.matrix(dat[!mask_col_ignore])
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
## created
count_filter<-function(dat,col_ignore=c(),min_max_frac=0.0,min_max=30,min_median=0,min_row_sum=500,other_cnt="other") {
  x<-split_count_df(dat,col_ignore)
  row_cnt = rowSums(x$count)
  row_sel = row_cnt >= min_row_sum
  cnt = x$count[row_sel,]
  row_cnt = row_cnt[row_sel]
  attr = x$attr[row_sel,]
  cnt_col_sel = cnt[,apply(row_normalize_matr(cnt),2,max) >= min_max_frac]
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

load.meta.choc <- function(file_name,counts.row.names) {
  meta = read.delim(file_name,header=T,stringsAsFactors=T)
  meta$Family = as.factor(meta$Family)
  
  meta$Sample.type.1 = as.factor(unlist(apply(meta,
                                              1,
                                              function(row) {switch(paste(row["Sample.type"],row["Therapy.Status"],sep="."),
                                                                    sibling.Before="sibling",
                                                                    patient.After="patient.after",
                                                                    patient.Before="patient.before")})))
  meta$Sample.ID.1 = as.factor(paste(meta$Subject.ID..blinded.,meta$Sample.type.1,sep="."))
  row.names(meta) = meta$sample_id_data
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
                      do.diversity=T,
                      file.name.root="std.plot",
                      file.name.sfx="png",
                      res.tests=NULL) {
  
  abund.meta.data.norm = row_normalize(abund.meta.data,col_ignore=abund.meta.attr.names)
  file.name.root = paste(file.name.root,label,sep=".")
  
  write.table(abund.meta.data,paste(file.name.root,".count.tab"),sep="\t")
  write.table(abund.meta.data.norm,paste(file.name.root,".freq.tab"),sep="\t")
  
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
    if (!is.null(res.tests$stab.feat)) {
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
    if ( do.diversity ) {
      pl.hist = plot.abund.meta(rich.meta,id.vars=id.vars,attr.names=abund.meta.attr.names,
                                file_name=make.graph.file.name("richness"))
      pl.hist = plot.abund.meta(div.meta,id.vars=id.vars,attr.names=abund.meta.attr.names,
                                file_name=make.graph.file.name("diversity"))
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
  moth.taxa <- read.mothur.taxa.summary("70cbd.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count_filter(taxa.lev.all,col_ignore=c(),min_max_frac=0.001,min_row_sum=50,other_cnt="other")
  meta = load.meta.choc("JCVI_ALL_sample_information.AT.txt")
  return (merge.counts.with.meta(taxa.lev,meta))
}

proc.choc <- function() {
  #DEBUG:
  do.power = F
  do.plots = T
  level = 6
  taxa.meta = read.choc(level)
  taxa.meta.data = taxa.meta$data
  taxa.meta.attr.names = taxa.meta$attr.names
  
  #subsample just the rows we need and then filter out clades that are all zero
  taxa.meta.data = count_filter(taxa.meta.data[c("LK3_91","LK4_92"),],
                                col_ignore=taxa.meta.attr.names,
                                min_max_frac=0.001,min_row_sum=50,other_cnt="other")
  
  if(do.power) {
    power.choc(taxa.meta.data,taxa.meta.attr.names)
  }
  
  if(do.plots) {
    
    std.plots(taxa.meta.data,taxa.meta.attr.names,id.vars.list=
                list(
                  c("Sample.type","Therapy.Status"),
                  c("Subject.ID..blinded.","Therapy.Status"),
                  c("Sample.type.1"),
                  c("Sample.ID.1")
                ),
              label=level
    )
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
    stable <- NA
    lpos.order <- NA
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
  selection[y$stable] = fwer.sel
  make.global(selection)
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
    ggtitle(main)
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
                            stability.transform.counts="ident") {
  
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
    cv.res = cv.glmnet.alpha(m_a$attr$T1D,count,family="binomial",standardize=standardize.glm)
    stopCluster(cl)
    
    make.global(cv.res)
    if(F) {
      stab.res.qdp = stability(m_a$count,m_a$attr$T1D=="T1D",weakness=0.7,mc.cores=1)
      make.global(stab.res.qdp)
      #stab.feat.qdp = plot(stab.res.qdp, sel.mode="PFER", cutoff=0.75, PFER=1, nvar=5, plot=FALSE)
      stab.feat.qdp = plot(stab.res.qdp, cutoff=0.75, PFER=1, nvar=5, plot=FALSE)
      make.global(stab.feat.qdp)
      print("")
      print (colnames(m_a$count)[stab.feat.qdp$selected])
    }
    penalty.alpha = cv.res$alpha
    #alpha = 0.8
    stab.res.c060 = stability.path(m_a$attr$T1D,count,weakness=0.9,
                                   family="binomial",steps=600,
                                   alpha=penalty.alpha,standardize=standardize.glm)
    make.global(stab.res.c060)
    fwer = alpha
    pi_thr = 0.6
    stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
    make.global(stab.feat.c060)
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
    #ad.res = adonis(count~T1D + Batch,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    #print(ad.res)
    ad.res.batch.t1d = adonis(count~Batch + T1D,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    print(ad.res.batch.t1d)
    ad.res.t1d = adonis(count~T1D,data=m_a$attr,permutations=n.adonis.perm,method=adonis.dist)
    print(ad.res.t1d)
    ad.res.paired = adonis(count~T1D,data=m_a$attr,strata=m_a$attr$Family,permutations=n.adonis.perm,method=adonis.dist)  
    print (ad.res.paired)
    #make.global(ad.res)
    
    #print (ad.res)
    #test.ad.res = ad.res$aov.tab$"Pr(>F)"[1]
    #r2.ad = ad.res$aov.tab$R2[1]
    
    gr.abs = m_a.abs$attr$T1D
    gr.abs.cnt = list(x=m_a.abs$count[gr.abs=="T1D",],y=m_a.abs$count[gr.abs!="T1D",])
    print(Xdc.sevsample(gr.abs.cnt))
    K = dim(gr.abs.cnt[[1]])[2]-1
    print (K)
    #make.global(gr.abs.cnt)
    print(Xmcupo.sevsample(gr.abs.cnt,K))
    
    ##this is a rank based pairwise test; column-wise standartization does not
    ##affect it, but rwo-wise - does
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
    fdr.res = fdrtool(test.wil.res,"pvalue")
    names(test.wil.res.adj) = names(test.wil.res)
    print(test.wil.res[test.wil.res<=alpha])
    print(test.wil.res.adj[test.wil.res.adj<=alpha])
    make.global(test.wil.res.adj)
    test.wil.res.adj.q = fdr.res$qval
    names(test.wil.res.adj.q) = names(test.wil.res)
    print(test.wil.res.adj.q[test.wil.res.adj.q<=alpha])
    make.global(test.wil.res.adj.q)
    
  }
  return (res)
}

load.meta.t1d <- function(file_name,counts.row.names,as.merged=F) {
  
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
  return (meta)
}


read.t1d <- function(taxa.level=3) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("43aa6.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary.txt")
  moth.taxa <- read.mothur.taxa.summary("ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count_filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
  #taxa.lev = taxa.lev.all
  meta = load.meta.t1d("annotation20130819.csv")
  return (merge.counts.with.meta(taxa.lev,meta))
}


proc.t1d <- function() {
  #taxa.levels = c(2,3,4,5,6)
  taxa.levels = c(6)
  do.std.plots = T
  do.tests = T
  for (taxa.level in taxa.levels) {
    taxa.meta = read.t1d(taxa.level)
    make.global(taxa.meta)
    taxa.meta.data = taxa.meta$data
    taxa.meta.attr.names = taxa.meta$attr.names
    label = paste("16s",taxa.level,sep=".")
    
    for (batch in list(c(1,2,3))) {
      batch.mask = (taxa.meta$data$Batch %in% batch)
      taxa.meta.data = taxa.meta$data[batch.mask,]
      taxa.meta.data = count_filter(taxa.meta.data,col_ignore=taxa.meta.attr.names,
                                    min_max_frac=0,min_max=1,min_row_sum=500,other_cnt="other")    
      make.global(taxa.meta.data)
      label = paste("16s",taxa.level,"b",paste(batch,collapse="-"),sep=".",collapse=".")
      
      print (paste("Working on",label))
      
      if (do.tests) {
        res.tests = try(
          test.counts.t1d(taxa.meta.data,taxa.meta.attr.names,
                          label=label,
                          stability.transform.counts="ident",
                          do.stability=F,
                          do.tests=T)
        )
      }
      if (do.std.plots) {
        try(
          std.plots(taxa.meta.data,taxa.meta.attr.names,id.vars.list=
                      list(
                        c("T1D","age.quant"),
                        #c("Family","T1D"),
                        #c("SampleID","Batch"),
                        c("T1D"),
                        c("T1D","Batch")
                      ),
                    label=label,
                    res.tests=res.tests
          )
        )
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
                                        stability.transform.counts="ident",
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
                                 stability.transform.counts="ident") {
  
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
    stab.res.c060 = stability.path(m_a$attr$sample.type,count,weakness=0.9,
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
                             stability.transform.counts="ident",
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

options(mc.cores=4)
options(boot.ncpus=4)
options(boot.parallel="snow")
#cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
#registerDoSNOW(cl)  

#power.nistal()
proc.choc()
#meta = load.meta.t1d("annotation20130819.csv")
#mgrast.dir = "../BATCH_01_02_META/BATCH_01-02_METAGENOMICS_MGRAST"
#mgr = read.mgrast.summary(paste(mgrast.dir,"cog.tsv",sep="/"),file_name.id.map=paste(mgrast.dir,"mgrast_to_samp_id.tsv",sep="/"))
#mgr.cnt = mgrast.to.abund.df(mgr,"level.2")
#taxa.meta = read.t1d(3)
#taxa.meta.data = taxa.meta$data
#taxa.meta.attr.names = taxa.meta$attr.names

#sink("analysis.log",split=T)

#proc.t1d()
#proc.t1d.mg()

#sink(NULL)

#print(plot.stability.selection.c060.at(s),rank="mean")
#spca.res = spca(count,K=6,type="predictor",sparse="varnum",trace=T,para=c(7,4,4,1,1,1))
##This has selected a single variable (Gordonibacter) from 16S level 3 with normalization "ident".
##Note that by default this function does not scale the predictor variables.
#hc.res <- get.biom(X = m_a$count, Y = m_a$attr$T1D, fmethod = c("studentt", "pls", "vip"), type = "HC")
#taxa.meta$data[taxa.meta$data$Batch==3,c("Gordonibacter_0.1.1.1.3.1.6","T1D")]
#taxa.meta$data[taxa.meta$data$Batch==3,c("Akkermansia_0.1.8.2.1.1.1","T1D")]

#taxa.meta = read.mr_oralc(3)
#proc.mr_oralc()
#power.pediatric.cancer.2013()
#stopCluster(cl)
