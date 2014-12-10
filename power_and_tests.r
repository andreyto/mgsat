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
#install.packages("kernlab")
#install.packages("ROCR")
#install.packages("caret")
#source("http://bioconductor.org/biocLite.R")
#biocLite("multtest")
#biocLite("GeneSelector")
#biocLite("RColorBrewer")
## To upgrade Bioconductor packages, use: source("http://bioconductor.org/biocLite.R"); biocLite("BiocUpgrade")
## Normal R-Studio upgrade apparently does not upgrade Bioconductor.

set_trace_options<-function(try.debug=T) {
  #tell our custom tryCatchAndWarn not to catch anything
  options(try.debug=try.debug)
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
  #options(error=recover)
}



## bind variable to the global environment
## can be used in debugging
make.global <- function(var) {
  assign(deparse(substitute(var)),var,envir=globalenv()) 
}


take_first<-function(x,n) {
  return (x[1:min(length(x),n)])
}

str_blank <- function(x) {
  return (nchar(gsub("\\s","",x)))
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

## Return row names as orderd factor (in the original order of rows)
rownames.ordered <- function(x) {
  rn = rownames(x)
  ordered(rn,levels=rn)
}

## Return column names as orderd factor (in the original order of columns)
colnames.ordered <- function(x) {
  rn = colnames(x)
  ordered(rn,levels=rn)
}

## Call quantcut and return its result as ordered factor
quantcut.ordered <- function(...) {
  ordered(quantcut(...))
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

## operation opposite to split_count_df()
join_count_df<-function(m_a) {
  
  if(!is.null(m_a$count)) {
    cbind(as.data.frame(m_a$count),m_a$attr)
  }
  else {
    m_a$attr
  }
  
}

as.dds.m_a <- function(m_a,formula.rhs,force.lib.size=T) {
  dds <- DESeqDataSetFromMatrix(countData = t(m_a$count),
                                colData = m_a$attr,
                                design = as.formula(paste("~",formula.rhs)))  
  
  if(force.lib.size) {
    ## from phyloseq vignette at 
    ## http://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  }
  
  return(dds)
}

norm.prop <- function(x, ...) {
  UseMethod('norm.prop', x)
}

norm.prop.default <- function(x.f) {
  x.f/sum(x.f)
}

norm.prop.matrix<-function(x.f,mar=1) {
  if(mar==1)
    x.f/rowSums(x.f)
  else if(mar==2)
    x.f/colSums(x.f)
  else
    stop("mar can be only 1 or 2")
}

norm.prop.m_a<-function(m_a,mar=1) {
  m_a$count = norm.prop(m_a$count,mar=mar)
  return (m_a)
}

norm.prop.data.frame<-function(dat,col_ignore=c(),mar=1) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  matnorm<-norm.prop(x$count,mar=mar)
  datnorm<-as.data.frame(matnorm) 
  res = cbind(x$attr,datnorm)
  #print(apply(res[!(names(res) %in% col_ignore)],1,sum) == 1)
  stopifnot(all(abs(apply(res[!(names(res) %in% col_ignore)],1,sum) - 1)<1e-6))
  res
}

norm.meta.data<-function(dat,col_ignore=c(),norm.func=NULL,...) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  if( is.null(norm.func) ) {
    norm.func = ihs
  }
  matnorm<-norm.func(x$count,...)
  datnorm<-as.data.frame(matnorm) 
  cbind(x$attr,datnorm)
}

# CLR transform is copied from SpiecEasi package
# https://github.com/zdk123/SpiecEasi
# If data is non-normalized count OTU/data table with samples on rows
# and features/OTUs in columns, then the transform is applied as
# clr(data)
# By default, it adds an offset of 1 before applying the transform, and acts
# on rows, returning matrix in the same order as input. The defaul base=2 in
# order to produce fold change between columns:
# (m[,k]/m[,l] == 2**(clr.mgsat(m)[,k]-clr.mgsat(m)[,l]))
# Note that in the original SpiecEasi implementation the transform should
# be applied as
# t(clr(data+1, 1, base=2))
# because it uses apply(mar=1) which transposes the result, and default base=exp(1)
# See `compositions` package for other Aitchison transforms

#' Centered log-ratio functions
#' @export
norm.clr <- function(x, ...) {
  UseMethod('norm.clr', x)
}

#' @method clr default
#' @export
norm.clr.default <- function(x.f, offset=1, base=2, tol=.Machine$double.eps) {
  x.f = x.f + offset
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
norm.clr.matrix <- function(x.f, mar=1, ...) {
  y = aaply(x.f, mar, norm.clr, ...)
  if(mar==2) {
    y = t(y)
  }
  return(y)
}

#' @method clr data.frame
#' @export
norm.clr.data.frame <- function(x.f, mar=1, ...) {
  as.data.frame(norm.clr(as.matrix(x.f), mar, ...))
}

#' @method clr on results of DESeq2 Rlog transform. x.f is not used
norm.clr.rlog.dds <- function(x.f,dds,...) {
  aaply(rlog.dds(dds,...),1,function(x) (x - mean(x)))
}

#' @method results of DESeq2 Rlog transform. x.f is not used
norm.rlog.dds <- function(x.f,dds,...) {
  rlog.dds(dds,...)
}

#' @method extract Rlog transformed values from DESeq2 object as sample-row matrix
rlog.dds <- function(dds,blind=T,fast=T,fitType="local") {
  t(assay(rlog(dds,blind=blind,fast=fast,fitType=fitType)))
}


## IHS (inverse hyperbolic sign) transform
## This is the same as log(x+(x**2+1)**0.5)
ihs <- function(x,theta=1) {
  asinh(theta*x)/theta
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

ihs.d1 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",1)
ihs.d2 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",2)

norm.ihs = ihs

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

norm.boxcox <- function(x, ...) {
  UseMethod('norm.boxcox', x)
}

norm.boxcox.default <- function(x.f) {
  boxcox.transform.vec(x.f)$x
}

norm.boxcox.matrix <- function(x.f, mar=2, ...) {
  y = aaply(x.f, mar, norm.boxcox, ...)
  if(mar==2) {
    y = t(y)
  }
  return(y)
}

norm.ihs.prop <- function(x,theta=1,mar=1) {
  norm.ihs(norm.prop(x,mar=mar),theta=theta)
}

## normalize raw count data according to one of the
## methods defined above.

norm.count <- function(x, ...) {
  UseMethod('norm.count', x)
}

norm.count.matrix <- function(count,method,drop.features=c("other"),method.args=list()) {
  count = do.call(method,c(
    list(count),
    method.args
  )
  )
  if(!is.null(drop.features)) {
    count = count[,!colnames(count) %in% drop.features]
  }
  return (count)
}

norm.count.m_a <- function(m_a,...) {
  m_a.norm = m_a
  m_a.norm$count = norm.count(m_a.norm$count,...)
  return(m_a.norm)
}

## subset method that will use the same subset argument on all data objects in m_a
subset.m_a <- function(m_a,subset,select.count=NULL,select.attr=NULL) {
  if(is.null(select.count)) select.count = T
  if(is.null(select.attr)) select.attr = T
  m_a$count = m_a$count[subset,select.count,drop=F]
  m_a$attr = m_a$attr[subset,select.attr,drop=F]
  return(m_a)
}

## If the other_cnt column is already present, it will be incremented with counts of clades
## relegated to the "other" in this call; otherwise, the new column with this name will be
## created.
## Count columns will be sorted in decreasing order of the column mean frequencies, so that
## you can easily subset the count matrix later to only keep N most abundant columns.
count.filter.m_a<-function(m_a,
                           min_max_frac=0.0,
                           min_max=0,
                           min_mean=0,
                           min_mean_frac=0.0,
                           min_row_sum=0,
                           max_row_sum=.Machine$integer.max,
                           other_cnt="other") {
  ##Note that filtering columns by a median value would not be a good idea - if we have a slightly
  ##unbalanced dataset where one group is 60% of samples and has zero presence in some column,
  ##and another group is 40% and has a large presence, then median filter will always throw this
  ##column away.
  #m.call = match.call()
  #make.global(m.call)
  #stop("DEBUG")
  x = m_a
  row_cnt = rowSums(x$count)
  row_sel = row_cnt >= min_row_sum & row_cnt < max_row_sum
  cnt = x$count[row_sel,]
  row_cnt = row_cnt[row_sel]
  attr = x$attr[row_sel,]
  cnt_norm = norm.prop(cnt)
  ind_col_ord = order(colSums(cnt_norm),decreasing=T)
  cnt = cnt[,ind_col_ord]
  cnt_norm = cnt_norm[,ind_col_ord]
  ind_col_sel = apply(cnt_norm,2,max) >= min_max_frac
  cnt = cnt[,ind_col_sel]
  cnt_norm = cnt_norm[,ind_col_sel]
  ind_col_sel = apply(cnt_norm,2,mean) >= min_mean_frac
  cnt = cnt[,ind_col_sel]
  cnt_norm = cnt_norm[,ind_col_sel]  
  cnt = cnt[,apply(cnt,2,max) >= min_max]
  cnt = cnt[,apply(cnt,2,mean) >= min_mean]
  
  cnt_col_other = as.matrix(row_cnt - rowSums(cnt))
  
  if (!all(cnt_col_other[,1]==0) && !is.null(other_cnt)) {
    colnames(cnt_col_other) = c(other_cnt)
    if (other_cnt %in% colnames(cnt)) {
      cnt[,other_cnt] = cnt[,other_cnt] + cnt_col_other[,other_cnt]
    }
    else {
      cnt = cbind(cnt,cnt_col_other)
    }
  }
  x$attr = attr
  x$count = cnt
  return(x)
}

count.filter<-function(dat,
                       col_ignore=c(),
                       ...) {
  ##Note that filtering columns by a median value would not be a good idea - if we have a slightly
  ##unbalanced dataset where one group is 60% of samples and has zero presence in some column,
  ##and another group is 40% and has a large presence, then median filter will always throw this
  ##column away.
  #m.call = match.call()
  #make.global(m.call)
  #stop("DEBUG")
  x<-split_count_df(dat,col_ignore)
  x = count.filter.m_a(x,...)
  return (join_count_df(x))
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
  if(getOption("try.debug",F)) {
    message("In debug mode, not catching anything in tryCatchAndWarn")
    return (eval.parent(expr))
  }
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
  freq = norm.prop(count.df,col_ignore=col_ignore)
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
  data = count.filter(data,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)
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
  ##replace=TRUE to select more samples than in original dataset
  ##TODO: apply strata to keep group count ratio
  ind.n = sample(indices, n, replace = TRUE, prob = NULL)
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
  #data = count.filter(data,col_ignore=data_factors,min_max_frac=0.1,min_row_sum=0)
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
    dm.freq = norm.prop(dm.counts)
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
  data = count.filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
  
  freq.mean = sorted.freq(data,col_ignore=attr_names)
  
  write.csv(as.matrix(sort(freq.mean[,"Saliva"],decreasing=TRUE)),"freq.top.saliva.csv")
  
  return(freq.mean)
}

sorted.freq.nistal<-function() {
  
  attr_names = c("group")
  data_all = read.nistal("children.txt")
  data = data_all
  #data = count.filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
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

read.mothur.otu.shared <- function(file_name) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  #where this X == NA column comes from, I have no idea. I do not see
  #it in the text file. Maybe R bug?
  data$X = NULL
  
  data$label = NULL
  data$numOtus = NULL
  row.names(data) = data$Group
  data$Group = NULL
  return (data)
}

read.mothur.cons.taxonomy <- function(file_name) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  row.names(data) = data$OTU
  data$Taxa = laply(strsplit(as.character(data$Taxonomy),"\\([0-9]*\\);"),function(x) x[[pmatch("unclassified",x,nomatch=length(x)+1,dup=T)-1]])
  return (data)
}

read.mothur.otu.with.taxa <- function(otu.shared.file,cons.taxonomy.file) {
  otu.df = read.mothur.otu.shared(otu.shared.file)
  taxa.df = read.mothur.cons.taxonomy(cons.taxonomy.file)
  stopifnot(all(names(otu.df) == taxa.df$OTU))
  names(otu.df) = paste(taxa.df$Taxa,taxa.df$OTU,sep=".")
  ## Order columns by taxa name
  otu.df = otu.df[,order(names(otu.df))]
  return (otu.df)
}

read.mothur.otu.with.taxa.m_a <- function(...) {
  count = as.matrix(read.mothur.otu.with.taxa(...))
  return(list(count=count))
}

make.mothur.taxa.summary.clade.names <- function(taxa.summary) {
  taxon = taxa.summary$taxon
  names(taxon) = taxa.summary$rankID
  ind_unclass = which(taxon == "unclassified")
  ind_unclass.ini = ind_unclass
  while(length(ind_unclass)>0) {
    taxon[ind_unclass] = taxon[unlist(strsplit(names(taxon)[ind_unclass], "\\.[0-9]+$"))]
    ind_unclass = ind_unclass[taxon[ind_unclass]=="unclassified"]
  }
  ind_unclass.ini = ind_unclass.ini[taxon[ind_unclass.ini] != "unknown"]
  taxon = as.character(taxon)
  taxon[ind_unclass.ini] = paste("Unclassified",taxon[ind_unclass.ini],sep="_")
  return(taxon)
}

read.mothur.taxa.summary <- function(file_name) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  #where this X == NA column comes from, I have no idea. I do not see
  #it in the text file. Maybe R bug?
  data$X = NULL
  data$clade = as.factor(make.mothur.taxa.summary.clade.names(data))
  return (data)
}


multi.mothur.to.abund.m_a <- function(data,level) {
  data.level = data[data$taxlevel==level,]
  attr = c("taxlevel","rankID","taxon","daughterlevels","total","clade")
  x = split_count_df(data.level,col_ignore=attr)
  row.names(x$count) = x$attr$clade
  x$count = t(x$count)
  return (x)
}

multi.mothur.to.abund.df <- function(data,level) {
  x = multi.mothur.to.abund.m_a(data,level)
  return (as.data.frame(x$count))
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
##TODO: remove intermediate conversion into data.frame
merge.counts.with.meta <- function(x,y,suffixes=c("","meta")) {
  mrg = merge(x,y,by="row.names",suffixes=suffixes)
  #Row.names column is generated by merge() when by="row.names"
  #the assignment below also serves as assertion that count records were
  #not duplicated as a result of merge
  row.names(mrg) = mrg$Row.names
  mrg$Row.names = NULL
  all.names = colnames(mrg)
  attr.names = all.names[!(all.names %in% colnames(x))]
  return (split_count_df(mrg,col_ignore=attr.names))
}

aggregate.by.meta.data.m_a <- function(m_a,
                                       group_col,
                                       count_aggr=sum,
                                       attr_aggr=NULL,
                                       group_col_result_name="SampleID") {
  x = m_a
  
  groups = list()
  groups[[group_col_result_name]] = x$attr[,group_col]
  
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
  x$attr = droplevels(x$attr)
  if(!is.null(x$count)) {
    count_names = names(x$count)
    x$count = aggregate(x$count,groups,count_aggr)  
    row.names(x$count) = x$count[,group_col_result_name]
    x$count = drop_columns(x$count,c(group_col_result_name))
    return (merge.counts.with.meta(x$count,x$attr))
  }
  else {
    return (x)
  }
  
}


aggregate.by.meta.data <- function(meta_data,
                                   col_ignore=c(),
                                   ...) {
  m_a = split_count_df(meta_data,col_ignore=col_ignore)
  m_a = aggregate.by.meta.data.m_a(m_a,...)
  return (join_count_df(m_a))
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

plot.abund.meta <- function(m_a,
                            id.vars=c(),
                            value.name="abundance",
                            file_name=NULL,
                            ggp.comp=NULL,
                            facet_grid.margins=FALSE,
                            clades.order=NULL,
                            geom="bar",
                            n.top=20,
                            id.var.dodge=NULL,
                            flip.coords=T,
                            sqrt.scale=F) {
  
  if(is.null(id.var.dodge)) {
    id.vars.facet = id.vars
  }
  else {
    id.vars.facet = id.vars[id.vars != id.var.dodge]
  }
  if(flip.coords && sqrt.scale) {
    stop("SQRT coordinate transformation cannot be used in flipped horizontal and vertical coords")
  }
  
  #TODO: change code in this method to use m_a directly
  data = join_count_df(m_a)
  attr.names = names(m_a$attr)
  
  dat = melt.abund.meta(data,id.vars=id.vars,attr.names=attr.names,value.name=value.name)
  if (is.null(clades.order)) {
    dat$clade = order.factor.by.total(dat$clade,dat[[value.name]])
  }
  else {
    dat$clade = factor(dat$clade,levels=clades.order,ordered=T)
  }
  
  ##show only n.top
  clades = levels(dat$clade)
  clades = clades[1:min(length(clades),n.top)]
  dat = dat[dat$clade %in% clades,]
  
  if(geom == "violin") {
    ##violin will fail the entire facet if one clade has zero variance,
    ##so we perturb data a tiny bit
    val = dat[,value.name]
    sd_dodge = max(abs(val))*1e-6
    dat[,value.name] = val + rnorm(length(val),0,sd_dodge)
  }
  
  #print(ddply(dat, c("T1D","clade"), summarise, mean_abund = mean(abundance)))
  
  
  if(is.null(id.var.dodge)) {
    fill="clade"
    color="clade"
  }
  else {
    fill=id.var.dodge
    color=id.var.dodge
  }
  ## here is what is going on with scale transformations in ggplot2 (v.1.0.0):
  ## scale_y_sqrt() - transforms data points before anything else is done like
  ## stats calculation or range detection; scale ticks are distributed quadratically
  ## (unevenly); this one can be combined with coord_flip().
  ## A (serious) downside of this method is that it can change relative magnitude
  ## of clade mean values as computed and shown by stat_summary(fun.y=mean) because large
  ## outliers will be compressed stronger.
  ## coord_trans(y = "sqrt") transforms only final geometric representation (including
  ## all glyphs) and creates a non-linear tick marks as well; it cannot be combined with
  ## coord_flip().
  ## Because of the above, we switch between two representations: vertical with
  ## coord_trans(y = "sqrt") and horizontal (flipped) with original linear coords
  aes_s = aes_string(x="clade",y=value.name,
                     fill = fill,color = color)
  gp = ggplot(dat, aes_s)
  
  if(geom == "bar") {
    gp = gp + stat_summary(fun.y=mean, geom="bar", aes(width=0.5), 
                           position=position_dodge(width=0.9))
    #geom_obj = stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text")
    #geom_obj = geom_bar(stat=stat_summary(fun.y="mean"),width=0.4)
  }
  else if(geom == "violin") {
    gp = gp + geom_violin(scale= "width", trim=TRUE, adjust=1)
  }
  else if(geom == "boxplot") {
    gp = gp + geom_boxplot(fill=NA,na.value=NA,outlier.size = 0)
  }
  else {
    stop(paste("Unexpected parameter value: geom = ",geom))
  }
  
  #stat_summary(fun.data = mean_cl_boot, geom = "pointrange",color="black")+
  #coord_flip()+
  #geom_boxplot(color="black")+
  #geom_point(position = "jitter")+
  #stat_identity(geom="bar")+
  #geom_bar(stat="identity")
  #scale_fill_brewer(type = "seq", palette = 1)
  #labels facet with number of cases
  if(flip.coords) {
    gp = gp + coord_flip()
  }
  if(sqrt.scale) {
    gp = gp + coord_trans(y = "sqrt")
  }
  if(length(id.vars.facet) == 0) {
    wr = facet_null()
  }
  else if (length(id.vars.facet) == 1) {
    wr = facet_wrap(as.formula(paste("~",id.vars.facet[1],sep="")))
  }
  else {
    wr = facet_grid(as.formula(paste(id.vars.facet[2],id.vars.facet[1],sep="~")),
                    drop=T,margins=facet_grid.margins)
  }
  gp = gp + wr
  if(length(id.vars.facet) > 0) {
    clade.names = as.character(clades)
    #this will be used to label each facet with number of cases in it
    facet.cnt <- ddply(.data=data, id.vars.facet, function(x,clade.names) 
    { c(.n=nrow(x),
        .y=mean(colMeans(as.matrix(x[,clade.names]),na.rm=T),
                names=F,na.rm=T)) },
    clade.names)
    facet.cnt$.n = paste("n =", facet.cnt$.n)
    #facet.cnt$y = facet.cnt$V2
    
    if(!is.null(id.var.dodge)) {
      legend.position = "right"
    }
    else {
      legend.position = "none"
    }
    make.global(facet.cnt)
    gp = gp +
      geom_text(data=facet.cnt, aes(x=2.5, y=.y, label=.n, size=2), colour="black", inherit.aes=FALSE, parse=FALSE)+
      theme(legend.position = legend.position,
            axis.title=element_blank(),
            axis.text.y=element_text(color=c("black","black")))
  }
  else {
    gp = gp + 
      theme(
        axis.title=element_blank(),
        axis.text.y=element_text(color=c("black","black")))
    
  }
  
  gp = gp + 
    theme(
      plot.title = element_text(size = rel(2)),
      axis.text.x = element_text(size = rel(0.85),angle=90),
      axis.text.y = element_text(size = rel(0.85)))
  
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

read.table.m_a <- function(file.base) {
  fn.count = paste(file.base,"count.tsv",sep=".")
  count = as.matrix(read.table(fn.count,
                               sep="\t",
                               header=T,
                               stringsAsFactors=T))
  fn.attr = paste(file.base,"attr.tsv",sep=".")
  attr = read.table(fn.attr,
                    sep="\t",
                    header=T,
                    stringsAsFactors=T)
  rownames(count) = attr$SampleID
  return (list(count=count,attr=attr))
}

write.table.m_a <- function(m_a,file.base,row.names=F) {
  fn.count = paste(file.base,"count.tsv",sep=".")
  write.table(m_a$count,
              fn.count,
              sep="\t",
              row.names = row.names)
  fn.attr = paste(file.base,"attr.tsv",sep=".")
  write.table(m_a$attr,
              fn.attr,
              sep="\t",
              row.names = row.names)
  return (list(fn.count=fn.count,fn.attr=fn.attr))
}

write.table.file.report.m_a = function(m_a,name.base,descr=NULL,row.names=F) {
  ## if we write row.names, Excel shifts header row to the left when loading
  file.base = report$make.file.name(name.base)
  files = write.table.m_a(m_a=m_a,
                          file.base=file.base,
                          row.names = row.names)
  if (!is.null(descr)) {
    report$add.descr(paste("Wrote counts and metadata for",
                           descr,
                           "to files",
                           paste(files,collapse=",")))
  }
  return (files) 
}


export.taxa.meta <- function(m_a,
                             label,
                             descr="",
                             row.proportions=T,
                             row.names=F) {
  write.table.file.report.m_a(m_a=m_a,
                              name.base=paste("samples.raw",label,sep="."),
                              descr=paste("raw counts",descr),
                              row.names=row.names)  
  
  if(row.proportions) {
    
    write.table.file.report.m_a(m_a=norm.prop.m_a(m_a),
                                name.base=paste("samples.proportions",label,sep="."),
                                descr=paste("proportions counts",descr),
                                row.names=row.names)
  }
}

## Generate index into x that selects equal number
## of elements reps for each level of x. If resp==0,
## reps will be set to the smallest level count.
## Adopted from balanced.specaccum {BiodiversityR}
## If there is a need to also balance within a specific
## strata (e.g. by phentype but preserving as many
## family connections as possible), then we will have
## to use methods from package `sampling`.
balanced.sample <- function(x, grouped = TRUE, reps = 0) {
  x = factor(x)
  n <- length(x)
  levs <- levels(x)
  minimum <- min(summary(x))
  if (reps > 0) {
    alllevs <- summary(x)
    goodlevs <- alllevs > (reps - 1)
    levs <- names(alllevs[goodlevs])
    minimum <- reps
  }
  nl <- length(levs)
  seq2 <- array(nl * minimum)
  seq1 <- sample(n)
  strat <- sample(nl)
  count <- 0
  for (i in 1:nl) {
    for (j in 1:n) {
      if (x[seq1[j]] == levs[strat[i]]) {
        count <- count + 1
        if (count > i * minimum) {
          count <- count - 1
        }
        seq2[count] <- seq1[j]
      }
    }
  }
  if (grouped == FALSE) {
    seq3 <- sample(seq2)
    seq2 <- seq3
  }
  return(seq2)
}

mgsat.richness.counts <- function(m_a,n.rar.rep=400) {
  
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  #S.ACE & all se.* give NaN often
  x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
              .final=function(x) (x/n.rar.rep)) %dopar% 
{estimateR(rrarefy(m_a$count,n.rar))[c("S.obs","S.chao1"),]}
return(list(e=t(x)))
}

## This uses incidence data and therefore should be applicable to both raw count data as 
## well as to proportions
## or other type of measurement where non-zero value means "present"
mgsat.richness.samples <- function(m_a,group.attr=NULL,n.rar.rep=400) {
  
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  if(is.null(group.attr)) {
    pool = factor(rep("All",nrow(m_a$count)))
    do.stratify = F
  }
  else {
    pool = m_a$attr[,group.attr]
    do.stratify = T
  }
  count = m_a$count
  ##somehow just supplying .combine="+" generates an error,
  ##but both the function below or skipping .combine and
  ##applying Reduce("+",...) on the returned list work fine
  plus <-function(x,y) (x+y)
  x = foreach(seq(n.rar.rep),.packages=c("vegan"),
              .combine=plus,
              .export=c("balanced.sample"),
              .final=function(x) (x/n.rar.rep)) %dopar% 
{
  if(do.stratify) {
    strat.ind = balanced.sample(pool)
    count = count[strat.ind,]
    pool=pool[strat.ind]
  }
  specpool(rrarefy(count,n.rar),pool=pool)
}

se.ind = grep(".*[.]se",names(x))
e = x[,-se.ind]
se = x[,se.ind]
names(se) = sub("[.]se","",names(se))
return(list(e=e,se=se))
}

## This returns Hill numbers
mgsat.diversity.alpha.counts <- function(m_a,n.rar.rep=400,is.raw.count.data=T) {
  
  require(vegan)
  n.rar = min(rowSums(m_a$count))
  
  f.div = function(m) { cbind(N1=exp(diversity(m,index="shan")),
                              N2 = diversity(m,index="invsimpson")) }
  
  if(is.raw.count.data) {
    x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
                .final=function(x) (x/n.rar.rep)) %dopar% 
{f.div(rrarefy(m_a$count,n.rar))}
  }
else {
  x = f.div(m_a$count)
}

if(is.raw.count.data) {
  #this is already unbiased (determenistic wrt rarefication)
  div.unb.simpson = rarefy(m_a$count,2)-1
  #div.fisher.alpha = fisher.alpha(m_a$count)
  x = cbind(x,div.unb.simpson=div.unb.simpson)
}

return(list(e=x))
}

mgsat.diversity.beta.dist <- function(m_a,n.rar.rep=400,method="-1") {
  
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
              .final=function(x) (x/n.rar.rep)) %dopar% 
{
  betadiver(rrarefy(m_a$count,n.rar),method=method)
}
return(list(e=x))
}

mgsat.diversity.beta <- function(m_a,n.rar.rep=400,method="-1",
                                 group.attr=NULL,
                                 betadisper.task=list(),
                                 adonis.task=NULL) {
  
  require(vegan)
  
  res = new_mgsatres()
  
  beta.dist = mgsat.diversity.beta.dist(m_a,n.rar.rep=n.rar.rep,method=method)$e
  
  method.help = paste(grep(sprintf('\"%s\"',method),capture.output(betadiver(help=T)),value=T),
                      ", where number of shared species in two sites is a, 
                      and the numbers of species unique to each site are b and c.",sep="")
  
  report$add.descr(sprintf("Computed beta-diversity matrix using function betadiver {vegan}
                   with method %s",
                           method.help))
  
  if(!is.null(group.attr)) {
    betadisp = do.call(betadisper,
                       c(
                         list(beta.dist,group=m_a$attr[,group.attr]),
                         betadisper.task
                       )
    )
    report$add.descr(sprintf("Results of function betadisper {vegan}
                       for the analysis of multivariate homogeneity of group dispersions.
                       This is applied to sample beta diversity matrix to analyze it with
                       respect to a grouping variable %s. Arguments for the call are: %s",
                             group.attr,
                             arg.list.as.str(betadisper.task)))
    anova.betadisp = anova(betadisp)
    report$add(anova.betadisp)
    res$anova.betadisp = anova.betadisp
    report$add(plot(betadisp),caption=sprintf("Results of betadisper {vegan}. Distances from samples 
               to the group
               centroids are shown in the first two principal coordinates.
               Groups are defined by the variable %s.
               Sample beta-diversity matrix was generated with method %s",
               group.attr,method.help))
  }
  
  if(!is.null(adonis.task)) {
    adonis.task$data.descr = sprintf("Beta-diversity dissimilarity matrix created with method %s",
                                     method.help)
    tryCatchAndWarn({
      m_a.bd = m_a
      m_a.bd$count = beta.dist
      res$adonis = do.call(test.counts.adonis.report,
                           c(
                             list(m_a=m_a.bd),
                             adonis.task
                           )
      )
    })
  }
  return(res)
}


mgsat.divrich.accum.plots <- function(m_a,is.raw.count.data=T) {
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  y = poolaccum(rrarefy(m_a$count,n.rar))
  report$add(plot(y),caption=sprintf("Accumulation curves for extrapolated richness indices 
        for random ordering of samples (function poolaccum of package vegan;
             estimation is based on incidence data). Samples were rarefied
             to the the minimum sample size (%s).",n.rar))
  
  if(is.raw.count.data) {
    y = estaccumR(rrarefy(m_a$count,n.rar))
    report$add(plot(y),caption=sprintf("Accumulation curves for extrapolated richness indices 
        for random ordering of samples (function estaccumR of package vegan;
             estimation is based on abundance data). Samples were rarefied
             to the the minimum sample size (%s).",n.rar))
  }
  
  y = specaccum(m_a$count,method="exact")
  report$add(plot(y, ci.type="polygon", ci.col="yellow",xlab="Size",ylab="Species"),
             caption=sprintf("Accumulation curve for expected number of species (features)
             for a given number of samples (function specaccum of package vegan,
             using method 'exact'. Samples were rarefied
             to the the minimum sample size (%s).",n.rar))
}

mgsat.divrich.counts.glm.test <- function(m_a.divrich,
                                          divrich.names=NULL,                                          
                                          formula.rhs,
                                          glm.args=list()) {
  if(is.null(divrich.names)) {
    divrich.names = colnames(m_a.divrich$count)
  }
  
  res = new_mgsatres()
  
  for(divrich.name in divrich.names) {
    
    family = "gaussian"
    form.str = sprintf("%s~%s",divrich.name,formula.rhs)
    do.model = F
    
    if(divrich.name=="N1") {
      descr = "Hill number N~1~ (equals exp(Shannon index)"
      do.model = T
    }
    else if(divrich.name=="N2") {
      descr = "Hill number N~2~ (equals Inverted Simpson index)"
      do.model = T
    }
    else if(divrich.name %in% c("S.obs","S.chao1","S.ACE")) {
      descr = paste("Richness estimate",divrich.name)
      do.model = T
    }
    
    if(do.model) {
      ##To see how well normal fits:
      ##library(fitdistrplus)
      ##descdist(x,boot=1000)
      ##gof = gofstat(fitdist(x,"norm"))
      ##gof$kstest #conservative, will properly reject only at high sample count
      mod = do.call(glm,
                    c(
                      list(
                        formula(form.str),
                        family=family,
                        data=cbind(as.data.frame(m_a.divrich$count),m_a.divrich$attr)
                      ),
                      glm.args
                    )
      )
      
      report$add(summary(mod),
                 caption=sprintf("Association of abundance based %s with sample metadata.
                                 GLM with family %s and formula %s",
                                 descr,family,form.str)
      )
      
      res[[divrich.name]] = mod
      
    }  
  }
  
  return(res)
  
}

mgsat.plot.richness.samples <- function(rich) {
  
  var.names = c("chao","jack1","boot")
  var.names.se = var.names
  e = rich$e[,var.names]
  e$Group = rownames.ordered(e)
  
  se = rich$se[,var.names.se]
  se$Group = rownames.ordered(se)
  
  e.m = melt(e,"Group",var.names,variable.name="Index",value.name="e")
  se.m = melt(se,"Group",var.names.se,variable.name="Index",value.name="se")
  
  data = join(e.m,se.m,by=c("Group","Index"))
  dodge <- position_dodge(width=0.9)  
  pl = ggplot(data, aes(x = Index, y = e, fill = Group)) +  
    geom_bar(position = dodge,stat="identity",width=0.8) + 
    geom_errorbar(position = dodge,aes(ymin=e-se, ymax=e+se,width=0.5)) +
    xlab("Index name") +
    ylab("Index estimate and standard error")
  
  return(pl)
}


mgsat.divrich.report <- function(m_a,
                                 n.rar.rep=400,
                                 is.raw.count.data=T,
                                 group.attr=NULL,
                                 counts.glm.task=NULL,
                                 beta.task=NULL,
                                 plot.profiles.task=list(),
                                 do.plot.profiles=T) {
  
  report.section = report$add.header("Abundance and diversity estimates",section.action="push")
  group.descr = ""
  group.descr.short = ""
  if(!is.null(group.descr)) {
    group.descr = sprintf(" Incidence-based estimates are computed on sample pools split by
                          metadata attribute %s, and in each repetition, samples are also
                          stratified to balance the number of samples at each level
                          of the grouping variable.", group.attr)
    group.descr.short = sprintf(" for samples grouped by %s",group.attr)
  }
  report$add.descr(sprintf("Counts are rarefied to the lowest library size, abundance-based and
                   incidence-based alpha diversity indices and richness estimates are computed.
                   This is repeated multiple times (n=%s), and the results are averaged.
                   Beta diversity matrix is also computed by averaging over multiple 
                   rarefications.%s",
                           n.rar.rep,group.descr))
  report$add.package.citation("vegan")
  
  res = new_mgsatres()

  res$rich.samples = mgsat.richness.samples(m_a,group.attr=group.attr,n.rar.rep=n.rar.rep)
  report$add(mgsat.plot.richness.samples(res$rich.samples),
             caption=sprintf("Incidence based rihcness estimates and corresponding standard errors%s",
                             group.descr.short)
  )

  if(is.raw.count.data) {
    res$rich.counts = mgsat.richness.counts(m_a,n.rar.rep=n.rar.rep)
  }
  
  res$div.counts = mgsat.diversity.alpha.counts(m_a,n.rar.rep=n.rar.rep,is.raw.count.data=is.raw.count.data)
  
  divrich.counts = res$div.counts$e
  if(is.raw.count.data) {
    divrich.counts = cbind(divrich.counts,res$rich.counts$e)
  }
  
  m_a.dr=list(count=as.matrix(divrich.counts),attr=m_a$attr)
  
  
  if(!is.null(counts.glm.task)) {
    
    res$glm.res = do.call(mgsat.divrich.counts.glm.test,
                          c(list(m_a.dr),
                            counts.glm.task
                          )
    )
    
  }
  
  if(do.plot.profiles) {
    do.call(plot.profiles,
            c(list(m_a=m_a.dr,
                   feature.descr=sprintf("Abundance-based diversity indices (Hill numbers) and richness estimates",
                                         group.descr.short),
                   value.name="index"),
              plot.profiles.task
            )
    )
  }
  
  report$push.section(report.section)
  mgsat.divrich.accum.plots(m_a,is.raw.count.data=is.raw.count.data)
  report$pop.section()
    
  if(!is.null(beta.task)) {
    res$beta = do.call(mgsat.diversity.beta,
                       c(list(m_a,
                              n.rar.rep=n.rar.rep,
                              group.attr=group.attr),
                         beta.task
                       )
    )
  }
  
  report$pop.section()
  
  return(res)
}



plot.profiles <- function(m_a,
                          feature.order=NULL,
                          id.vars.list=list(c()),
                          clade.meta.x.vars=c(),
                          do.profile=T,
                          do.clade.meta=T,
                          value.name="abundance",
                          show.profile.task=list(
                            geoms=c("bar","violin","boxplot"),
                            dodged=T,
                            faceted=T
                          ),
                          show.clade.meta.task=list(),
                          feature.descr="abundance",
                          sqrt.scale=F) {
  
  report.section = report$add.header(sprintf("Plots of %s in multiple representations",feature.descr),section.action="push")
  report$add.descr("Plots are shown with relation to various combinations of meta 
                   data variables and in different graphical representations. Lots of plots here.")
  
  report$add.header("Iterating over all combinations of grouping variables")
  report$push.section(report.section)
  
  if(is.null(feature.order)) {
    feature.order = list(list(ord=NULL,ord_descr="original"))
  }
  
  for (id.vars in id.vars.list) {
    
    report$add.header(paste("Grouping variables",paste0(id.vars,collapse=",")))
    
    
    report$add.header(sprintf("Iterating over %s profile sorting order",feature.descr))
    report$push.section(report.section)
    
    for(pl.par in feature.order) {
      report$add.header(sprintf("%s profile sorting order: %s",feature.descr,pl.par$ord_descr))
      if(do.clade.meta) {
        clade.names.meta=if(is.null(pl.par$ord)) colnames(m_a$count) else pl.par$ord
        clade.names.meta = clade.names.meta[1:min(length(clade.names.meta),10)]
        
        report$add.header("Iterating over meta data variables")
        report$push.section(report.section)
        
        for(x.var in clade.meta.x.vars) {
          
          group.var = id.vars[id.vars != x.var]
          if(length(group.var)>0) {
            group.var = group.var[1]
          }
          else {
            group.var = NULL
          }
          tryCatchAndWarn({
            do.call(show.clade.meta,
                    c(
                      list(m_a=m_a,
                           clade.names=clade.names.meta,
                           x.var=x.var,
                           group.var=group.var,
                           value.name=value.name,
                           vars.descr=feature.descr),
                      show.clade.meta.task
                    )
            )
          })
          
        }
        
        report$pop.section()
      }
      
      if(do.profile) {
        
        report$add.header("Iterating over dodged vs faceted bars")
        report$add.descr("The same data are shown in multiple combinations of graphical representations. 
                         This is the same data, but each plot highlights slightly different aspects of it.
                         It is not likely that you will need every plot - pick only what you need.")
        report$push.section(report.section)
        
        within(show.profile.task, {
          if(!(dodged || faceted)) {
            faceted = T
          }
        })
        id.vars.dodge = list()
        if(show.profile.task$faceted) {
          id.vars.dodge[["faceted"]] = list(dodge=NULL,descr="faceted")
        }
        if(show.profile.task$dodged) {
          id.vars.dodge[["dodged"]] = list(dodge=id.vars[1],descr="dodged")
        }
        
        for(id.var.dodge in id.vars.dodge) {
          
          report$add.header(paste(id.var.dodge$descr,"bars. Iterating over orientation and, optionally, scaling"))
          report$push.section(report.section)
          
          for(other.params in list(
            list(flip.coords=T,
                 sqrt.scale=F,
                 descr="in flipped orientation, not scaled"),
            list(flip.coords=F,
                 sqrt.scale=sqrt.scale,
                 descr=paste("in original orientation", if(sqrt.scale) ", SQRT scaled" else "", sep=""))
          )) {
            
            report$add.header(paste(sprintf("%s are ",feature.descr), other.params$descr, ". Iterating over bar geometry",sep=""))
            report$push.section(report.section)
            
            for(geom in show.profile.task$geoms) {
              tryCatchAndWarn({
                pl.hist = plot.abund.meta(m_a=m_a,
                                          id.vars=id.vars,
                                          clades.order=pl.par$ord,
                                          geom=geom,
                                          file_name=NULL,
                                          id.var.dodge=id.var.dodge$dodge,
                                          flip.coords=other.params$flip.coords,
                                          sqrt.scale=other.params$sqrt.scale,
                                          value.name=value.name
                )
                #env=as.environment(as.list(environment(), all.names=TRUE))
                #print(names(as.list(env)))
                #print(evals("pl.hist",env=env))
                report$add(pl.hist,
                           caption=paste(sprintf("%s grouped by ",feature.descr),
                                         paste(id.vars,collapse=","),
                                         if(!is.null(pl.par$ord)) 
                                         {sprintf("in %s sorting order",pl.par$ord_descr)} 
                                         else {""})
                )
              })
            }
            report$pop.section()
          }
          report$pop.section()        
        }
        report$pop.section()
      }
    }
    report$pop.section()
  }
  report$pop.section()
  report$pop.section()
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
  data_all = count.filter(data_all,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)  
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
  
  
  abund_matr_norm = norm.prop(count_matr_from_df(taxa.meta.data,col_ignore=taxa.meta.attr.names))
  
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
  moth.taxa <- read.mothur.taxa.summary("4b92c86afff1e8ed5443fd52f540eb60.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count.filter(taxa.lev.all,col_ignore=c(),min_max_frac=0.001,min_row_sum=50,other_cnt="other")
  meta = load.meta.choc("CHOC_MiSEQ_June2014.txt")
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
    taxa.meta.data = count.filter(taxa.meta.data,
                                  col_ignore=taxa.meta.attr.names,
                                  min_max_frac=0.001,min_row_sum=50,other_cnt="other")
    
    if(do.power) {
      power.choc(taxa.meta.data,taxa.meta.attr.names)
    }
    res.tests = NULL
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


gen.tasks.choc.old <- function(taxa.meta) {
  
  glmer.descr.tpl = "The random effect terms in the model formula describe that: 
  - intercept varies among families 
  and among individuals within families (nested random effects);
  - when there are multiple observations (samples) per
  individual, we add a random effect for each observation to account
  for the overdispersion;
  The fixed effect is %s. The binomial family
  is used to build a set of univariate models, with each
  model describing the observed counts of one clade.
  P-values are estimated from the model under a null hypothesis
  of zero coefficients and a two-sided alternative. 
  Benjamini & Hochberg (1995) method is used 
  for multiple testing correction, and the significant clades
  are reported."
  
  task = within(
    list(),
{
  do.std.plots = T
  do.clade.meta=T
  do.profile=T
  
  do.heatmap = T
  do.std.tests = T
  
  do.stability=T
  do.tests=T
  do.genesel=T
  do.glmnet=T
  do.glmer=T
  do.adonis=T
  
  stability.steps = 600
  n.adonis.perm = 4000  
  
  resp.attr = "NULL" 
  stability.model.family = "binomial"
  
  genesel.resp.attr = resp.attr
  
  stability.transform.counts="ihs"
  
  adonis.tasks = list()
  
  glmer.task = list()
  
  plot.group = list()            
  
  clade.meta.x.vars=c()
  
  heatmap.task = list()
  
  
}
  )

task1 = within(
  task,
{
  
  descr = "All samples, no aggregation"
  
  do.stability=T
  do.tests=F
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = resp.attr
  
  plot.group = list(
    c("TherapyStatus","Sample.type"),
    c("Sample.type","visit")
  )            
  
  clade.meta.x.vars=c("visit","age")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type.1","visit","age","Antibiotic"),
    attr.row.labels="SampleID"
  )
  
}
)

task2 = within(
  task,
{
  
  descr = "Samples before therapy aggregated by SubjectID"
  
  do.clade.meta=F
  
  taxa.meta.aggr = aggregate.by.meta.data(subset(taxa.meta$data,TherapyStatus=="before.chemo"),
                                          "SubjectID",
                                          taxa.meta$attr.names)
  resp.attr = "Sample.type" 
  stability.model.family = "binomial"
  
  genesel.resp.attr = resp.attr
  
  adonis.tasks = list(
    list(formula_rhs="Sample.type",
         strata=NULL,
         descr="Association with the patient/sibling status unpaired"),
    list(formula_rhs="Sample.type",
         strata="FamilyID",
         descr="Association with the patient/sibling status paired")
  )
  
  glmer.task = list(
    formula.rhs = "Sample.type + (1|FamilyID/SubjectID)",
    linfct=c("Sample.typesibling = 0"),
    descr=sprintf(glmer.descr.tpl,"Sample.type")
  )
  
  plot.group = list(
    c("Sample.type")
    #c("Sample.type","TherapyStatus"),
    #c("Sample.type","visit")
    #c("SubjectID","TherapyStatus"),
    #c("Sample.type.1")
    #c("SampleID.1")
  )            
  
  clade.meta.x.vars=c("visit")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type","Antibiotic"),
    attr.row.labels="SampleID"
  )
  
}
)


task3 = within(
  task,
{
  
  descr = "Patients only, no aggregation"
  
  do.std.plots = F
  do.clade.meta=T
  do.profile=T
  
  do.heatmap = F
  do.std.tests = T
  
  do.stability=T
  do.tests=T
  do.genesel=T
  do.glmnet=T
  do.glmer=T
  do.adonis=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="patient")
  
  resp.attr = "visit" 
  stability.model.family = "gaussian"
  
  genesel.resp.attr = "TherapyStatus" 
  
  adonis.tasks = list(
    list(formula_rhs="visit",
         strata=NULL,
         descr="Association with the visit number unpaired"),
    list(formula_rhs="visit",
         strata="SubjectID",
         descr="Association with the visit number paired by subject")
  )
  
  glmer.task = list(
    formula.rhs = "visit + (1|SubjectID) + (1|SampleID)",
    linfct=c("visit = 0"),
    descr=sprintf(glmer.descr.tpl,"visit")
  )
  
}
)

task4 = within(
  task,
{
  
  descr = "All siblings and patients with two or more visits, no aggregation"
  
  do.stability=F
  do.tests=F
  do.clade.meta=F
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="sibling" | visit.max>=2)
  
  resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = resp.attr
  
  plot.group = list(
    c("Sample.type","visit")
  )            
  
  clade.meta.x.vars=c("visit","age")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type.1","visit","age","Antibiotic"),
    attr.row.labels="SampleID"
  )
  
}
)

task5 = within(
  task,
{
  
  descr = "Patients only, no aggregation. Checking antibiotic use."
  
  do.stability=F
  do.tests=F
  do.clade.meta=T
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="patient")
  
  resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = resp.attr
  
  plot.group = list(
    c("Antibiotic","visit")
  )            
  
  clade.meta.x.vars=c("visit")
  
  heatmap.task = list(
    attr.annot.names=c("visit","Fever","Antibiotic"),
    attr.row.labels="SampleID"
  )
  
}
)


return (list(task1,task2,task3,task4,task5))

}


proc.choc <- function() {
  #taxa.levels = c(2,3,4,5,6)
  taxa.levels = c(2)
  
  do.summary.meta = T
  
  report$add.descr("Set of analysis routines is applied
                   in nested loops over combinations of defined subsets of samples (if any) and 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out its place in the report hierarchy.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture.")
  
  report.section = report$get.section()
  
  if(do.summary.meta) {
    
    report$add.header("Summary of metadata variables")
    
    taxa.meta = read.choc(taxa.level=2)
    
    xtabs.formulas = list("~Sample.type+TherapyStatus","~Sample.type+visit","~FamilyID","~Sample.type.1","~SubjectID")
    for(xtabs.formula in xtabs.formulas) {
      fact.xtabs = xtabs(as.formula(xtabs.formula),data=taxa.meta$data,drop.unused.levels=T)
      report$add.table(fact.xtabs,show.row.names=T,caption=paste("Sample cross tabulation",xtabs.formula))
      report$add.printed(summary(fact.xtabs))
    }
    
    with(taxa.meta$data,{
      report$add.printed(summary(aov(age~Sample.type)),
                         caption="ANOVA for age and sample type")
      report$add(qplot(Sample.type,age,geom="violin"),
                 caption="Violin plot for age and sample type")
    })
    
    with(taxa.meta$data,{
      report$add.printed(cor.test(age,
                                  visit,
                                  method="spearman"),
                         caption="Spearman RHO for age and visit")
      
    })
    with(taxa.meta$data[taxa.meta$data$Sample.type=="patient",],{
      report$add.printed(cor.test(age,
                                  visit,
                                  method="spearman"),
                         caption="Spearman RHO for age and visit, patients only")
      
    })
    
    #summary(glht(lm(age~visit,data=meta[meta$Sample.type=="patient",]),linfct="visit=0"))
    #summary(glht(lmer(age~visit+(visit|Sample.type),data=meta),linfct="visit=0"))
    report$add(ggplot(taxa.meta$data,aes(x=visit,y=age,color=Sample.type))+
                 geom_point()+
                 stat_smooth(method="loess", se = T,degree=1,size=1),
               caption="Plot for age and visit with Loess trend line")
    
  }
  
  
  report$add.header("Iterating over taxonomic levels")
  report$push.section(report.section)
  
  for (taxa.level in taxa.levels) {
    taxa.meta = read.choc(taxa.level)
    
    label = paste("16s","l",taxa.level,sep=".",collapse=".")
    report$add.header(paste("Taxonomic level:",taxa.level))
    report$push.section(report.section)
    
    make.global(taxa.meta)
    #taxa.meta$data = taxa.meta$data[taxa.meta$data$Sample.type.1 != "sibling",]
    
    report$add.p(paste("Number of samples:",nrow(taxa.meta$data)))
    
    report$add.header("Iterating over subsets of data")
    report$push.section(report.section)
    
    tasks = gen.tasks.choc(taxa.meta)
    
    for (task in tasks) {
      
      report$add.header(paste("Subset:",task$descr))
      report$push.section(report.section)
      
      make.global(task)
      
      with(task,{
        
        report$add.p(paste("After aggregating/subsetting samples:",nrow(taxa.meta.aggr$data)))
        
        taxa.meta.aggr$data = count.filter(taxa.meta.aggr$data,
                                           col_ignore=taxa.meta.aggr$attr.names,
                                           min_max=10)
        #taxa.meta.data = count.filter(taxa.meta.data,col_ignore=taxa.meta.attr.names,
        #                              min_mean=0.002,
        #                              min_max_frac=0.1,min_max=10,
        #                              min_row_sum=500,other_cnt="other")
        
        
        report$add.p(paste("After count filtering,",
                           (ncol(taxa.meta.aggr$data)-length(taxa.meta.aggr$attr.names)),"clades left."))
        
        
        res.tests = NULL
        if (do.std.tests) {
          res.tests = tryCatchAndWarn(
            test.counts.project(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,
                                label=label,
                                do.stability=do.stability,
                                do.tests=do.tests,
                                do.genesel=do.genesel,
                                do.glmnet=do.glmnet,
                                do.glmer=do.glmer,
                                do.adonis=do.adonis,
                                genesel.resp.attr = genesel.resp.attr,
                                resp.attr = resp.attr,
                                stability.model.family = stability.model.family,
                                stability.transform.counts=stability.transform.counts,
                                adonis.tasks = adonis.tasks,
                                glmer.task = glmer.task,
                                stability.steps = stability.steps,
                                n.adonis.perm = n.adonis.perm
                                
            )
          )
        }
        if( do.std.plots || do.heatmap ) {
          report$add.header("Plots")
        }
        if (do.std.plots) {
          
          tryCatchAndWarn({
            std.plots(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,
                      id.vars.list=plot.group,
                      label=label,
                      res.tests=res.tests,
                      do.clade.meta=do.clade.meta,
                      do.profile=do.profile,
                      #clade.meta.x.vars=c("Specimen.Collection.Date","age")                  
                      clade.meta.x.vars=clade.meta.x.vars
            )
          })
        }
        if (do.heatmap) {
          tryCatchAndWarn({
            heatmap.counts(taxa.meta.aggr$data,
                           taxa.meta.aggr$attr.names,
                           attr.annot.names=heatmap.task$attr.annot.names,
                           attr.row.labels=heatmap.task$attr.row.labels,
                           caption="Heatmap of abundance profile",
                           stand.show="max",
                           trans.show=sqrt
            )
          })
        }
      })
      report$pop.section()
    }
    report$pop.section()
    report$pop.section()
  }
  report$pop.section()
}

## Default method for loading metadata file
## This will work OK if your file follows certain conventions:
## 1. Tab delimited with a header row
## 2. All strings can be treated as factors
## 3. All fields that look like numbers are OK to treat like numbers.
##    Specifically, make sure that all ID fields (SampleID, SubjectID etc)
##    start with a letter, not a number.
## 4. There has to be a unique key field SampleID, that matches SampleID in abundance tables
##
## If your file cannot follow the conventions above, then write your own method and pass it to 
## read.data.project(). The method should create SampleID factor key field and also set row.names to
## the values of that field.

load.meta.default <- function(file.name) {
  meta = read.delim(file_name,header=T,stringsAsFactors=T)
  meta$SampleID = as.factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  return (meta)
}


read.data.project.yap <- function(taxa.summary.file,
                                  otu.shared.file,
                                  cons.taxonomy.file,
                                  meta.file,
                                  load.meta.method,
                                  load.meta.options=list(),
                                  count.filter.options=NULL,
                                  taxa.level=3) {
  if (taxa.level == "otu") {
    taxa.lev.all = read.mothur.otu.with.taxa.m_a(otu.shared.file=otu.shared.file,cons.taxonomy.file=cons.taxonomy.file)
    count.file = otu.shared.file
  }
  else {
    moth.taxa <- read.mothur.taxa.summary(taxa.summary.file)
    taxa.lev.all = multi.mothur.to.abund.m_a(moth.taxa,taxa.level)
    count.file = taxa.summary.file
  }  
  report$add.p(sprintf("Loaded %i records for %i clades from count file %s for taxonomic level %s",
                       nrow(taxa.lev.all$count),ncol(taxa.lev.all$count),count.file,taxa.level))
  
  if(!is.null(count.filter.options)) {
    report$add.p(paste("Filtering initial records with arguments",arg.list.as.str(count.filter.options)))
    report$add.p("Note that many community richness estimators will not work correctly 
                 if provided with abundance-filtered counts")
    taxa.lev = do.call(count.filter.m_a,c(list(taxa.lev.all),count.filter.options))
    report$add.p(sprintf("After filtering, left %i records for %i clades for taxonomic level %s",
                         nrow(taxa.lev$count),ncol(taxa.lev$count),taxa.level))
  }
  else {
    ##this should only order features by mean abundance but not drop anything
    taxa.lev = count.filter.m_a(taxa.lev.all)
  }
  meta = do.call(load.meta.method,c(file.name=meta.file,load.meta.options))  
  m_a = merge.counts.with.meta(taxa.lev$count,meta)
  report$add.p(sprintf("After merging with metadata, %i records left",
                       nrow(m_a$count)))
  return (m_a)
}


summary.meta.method.default <- function(taxa.meta) {
  
  report$add.header("Summary of metadata variables")
  
  m_a = split_count_df(taxa.meta$data,col_ignore=taxa.meta$attr.names)
  
  report$add.printed(summary(m_a$attr),caption="Summary of metadata variables")
  
}


mgsat.16s.task.template = within(list(), {
  
  main.meta.var = "Group"
  
  read.data.method=read.data.project.yap  
  
  read.data.task = list(
    taxa.summary.file=NULL,
    otu.shared.file=NULL,
    cons.taxonomy.file=NULL,
    meta.file=NULL,
    load.meta.method=load.meta.default,
    load.meta.options=list()
  )
  
  taxa.levels = c("otu",2,3,4,5,6)
  
  get.taxa.meta.aggr<-function(m_a) { return (m_a) }
  
  count.filter.sample.options=list(min_row_sum=400,max_row_sum=100000)
  
  summary.meta.method=summary.meta.method.default
  
  do.summary.meta = F
  
  do.plots = T
  
  do.tests = T
  
  summary.meta.task = list(
    meta.x.vars = c(),
    group.var = c(main.meta.var)
  )
  
  test.counts.task = within(list(), {
    
    do.deseq2 = T
    do.genesel=T
    do.stabsel=T
    do.glmer=T
    do.adonis=T
    do.divrich=c("otu",6)
    do.plot.profiles.abund=T
    do.heatmap.abund=T
    do.select.samples=c()
    
    alpha = 0.05
    
    divrich.task = within(list(),{
      n.rar.rep=400
      is.raw.count.data=T
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })
      beta.task = within(list(),{
        method="-1"
        betadisper.task=list()
        adonis.task=NULL #will be replaced by global adonis.task
      })
      do.plot.profiles = T
    })
    
    count.filter.feature.options=list(min_mean_frac=0.0005)
    
    norm.count.task = within(list(), {
      method = "norm.ihs.prop"
      method.args = list()
      drop.features=list("other")
      ##method="norm.rlog"
      ##method.args=list(dds=NA) #signals to pull Deseq2 object
    })
    
    deseq2.task = within(list(), {
      formula.rhs=main.meta.var
      test.task=list(
        fitType="local"
      )
    })
    
    stabsel.task = list(
      resp.attr = main.meta.var,      
      args.fitfun = list(
        family="binomial",
        standardize=T                                     
      ),
      args.stabsel = list(
        PFER=0.05,
        sampling.type="SS",
        assumption="r-concave"
      )
    )
    
    
    genesel.task = list(
      resp.attr = main.meta.var
    )
    
    select.samples.task = list (
      n.species = 20,
      n.samples = 20
    )
    
    adonis.task = within(list(), {
      
      tasks=list(list(
        formula.rhs=main.meta.var,
        strata=NULL,
        descr=paste("Association with",main.meta.var)
      ))
      
      n.perm=4000
      dist.metr="bray"
      col.trans="range"
    })
    
    glmer.task = within(list(), {
      
      tasks = list(list(
        descr.extra = "",
        formula.rhs = paste(main.meta.var,"(1|SampleID)",sep="+"),
        linfct=c("YearsSinceDiagnosis = 0")
      ))
      
    })
    
    plot.profiles.task = within(list(), {
      id.vars.list = list(c(main.meta.var))
      clade.meta.x.vars=c()
      do.profile=T
      do.clade.meta=F
      show.profile.task=list(
        geoms=c("bar","violin","boxplot"),
        dodged=T,
        faceted=T
      )
      show.clade.meta.task=list()      
    })
    
    plot.profiles.abund.task = within(list(), {
      
      norm.count.task = within(list(), {
        method = "norm.prop"
        drop.features=c("other")
      })
      
    })
    
    heatmap.abund.task = within(list(), {
      attr.annot.names=c(main.meta.var)
      attr.row.labels=NULL
      stand.clust="standardize"
      dist.metr="euclidean"
      caption="Heatmap of abundance profile"
      stand.show="range"
    })
    
  })
  
})


start.cluster.project <- function() {
  cl<-makeCluster(getOption("mc.cores", 2L),type = "SOCK") #number of CPU cores
  registerDoSNOW(cl)
  return(cl)
}

stop.cluster.project <- function(cl) {
  stopCluster(cl)
}

## create new MGSAT result object
new_mgsatres <- function(...) {
  x = list(...)
  class(x) <- append(class(x),"mgsatres")
  return(x)
}

## extract ranking of features created by a specific method
## Value: named list with at least name `ranked` that is an ordered 
## vector of feature names.
get.feature.ranking <- function(x, ...) {
  UseMethod('get.feature.ranking', x)
}

get.feature.ranking.default <- function(x.f) {
  stop("Not defined for arbitrary objects")
}

get.feature.ranking.stabsel <- function(x.f,only.names=T) {
  ranked = x.f$max[order(x.f$max,decreasing = T)]
  if(only.names) {
    ranked = names(ranked)
  }
  return(list(ranked=ranked))
}

get.feature.ranking.mgsatres <- function(x.f,method="stabsel") {
  meth.res = x.f[[method]]
  if(method=="stabsel") {
    res = get.feature.ranking(meth.res)
  }
  else {
    stop(paste("I do not know what to do for method",method))
  }
  return (res)
}

count.filter.report <- function(m_a,count.filter.options,descr) {
  m_a = do.call(count.filter.m_a,
                c(list(m_a),
                  count.filter.options)
  )
  
  report$add.p(paste("After",descr,"count filtering with arguments",
                     arg.list.as.str(count.filter.options),
                     ",",
                     ncol(m_a$count),
                     "features and",
                     nrow(m_a$count),
                     "samples left."))
  return(m_a)
}

proc.project <- function(
  task.generator.method
) {
  
  report$add.descr("Set of analysis routines is applied
                   in nested loops over combinations of defined subsets of samples (if any) and 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out its place in the report hierarchy.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture.")
  
  report.section = report$get.section()
  
  report$add.header("Iterating over subsets of data")
  report$push.section(report.section) #1 {
  
  tasks = task.generator.method()
  
  cl = start.cluster.project()
  
  res.tasks = lapply(tasks,function(task) {
    
    report$add.header(paste("Subset:",task$descr))
    report$push.section(report.section) #2 {
    
    res.task = new_mgsatres(task=task)
    
    res.task$res.taxa.levels = with(task,{
      
      report$add.header("Iterating over taxonomic levels")
      report$push.section(report.section) #3 {
      
      res.taxa.levels = list()
      
      for (taxa.level in taxa.levels) {
        
        res.level = new_mgsatres(taxa.level=taxa.level)
        
        label = paste("16s","l",taxa.level,sep=".",collapse=".")
        report$add.header(paste("Taxonomic level:",taxa.level))
        report$push.section(report.section) #4 {
        
        m_a = do.call(read.data.method,
                      c(
                        list(taxa.level=taxa.level),
                        read.data.task
                      )
        )
        
        m_a = get.taxa.meta.aggr(m_a)
        
        report$add.p(paste("After aggregating/subsetting, sample count is:",nrow(m_a$count)))
        
        m_a = count.filter.report(m_a,
                                  count.filter.options=count.filter.sample.options,
                                  "sample")
        
        if(do.summary.meta) {
          summary.meta.method(m_a)
          do.call(report.sample.count.summary,c(
            list(m_a),
            summary.meta.task
          )
          )
          ##only do summary once
          do.summary.meta = F
        }
        
        export.taxa.meta(m_a,
                         label=label,
                         descr=descr,
                         row.proportions=T,
                         row.names=F)
        
        res.tests = NULL
        
        if (do.tests) {
          
          ## modify a copy
          test.counts.task.call = test.counts.task
          test.counts.task.call$do.select.samples = (taxa.level %in% test.counts.task.call$do.select.samples)
          test.counts.task.call$do.divrich = (taxa.level %in% test.counts.task.call$do.divrich)
          
          res.tests = tryCatchAndWarn(
            do.call(test.counts.project,
                    c(
                      list(
                        m_a=m_a,
                        label=label
                      ),
                      test.counts.task.call
                    )
            )
          )          
        }
        
        res.level$res.tests = res.tests
        
        report$pop.section() #4 }
        res.taxa.levels[[as.character(taxa.level)]] = res.level
      }
      report$pop.section() #3 }
      res.taxa.levels
    })
    
    report$pop.section() #2 }
    
    res.task
  })
  
  report$pop.section() #1 }
  
  stop.cluster.project(cl)
  
  res = new_mgsatres(res.tasks=res.tasks)
  
  return (res)
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

make.sample.summaries <- function(m_a.abs) {
  x = data.frame(
    count.sum = rowSums(m_a.abs$count)
  )
  return (list(count=x,attr=m_a.abs$attr))
}


show.sample.summaries.meta <- function(m_a,
                                       x.var,
                                       group.var,
                                       summ.names=NULL,
                                       value.name="Y",
                                       trans="ident",
                                       vars.descr="sample summary properties") {
  if(is.null(summ.names)) {
    summ.names = names(m_a$count)
  }
  show.clade.meta(m_a=m_a,
                  clade.names=summ.names,
                  x.var=x.var,
                  group.var=group.var,
                  value.name=value.name,
                  trans=trans,
                  vars.descr=vars.descr) 
  
}


show.clade.meta <- function(m_a,
                            clade.names,
                            x.var,
                            group.var,
                            value.name="abundance",
                            trans="boxcox",
                            vars.descr="abundances") {
  
  
  count = m_a$count[,clade.names,drop=F]
  
  count = switch(trans,
                 boxcox=norm.boxcox(count),
                 ihs=ihs(count,1),
                 ident=count,
                 binary=(count > 0))  
  
  if(trans != "ident") {
    trans.msg = paste("(transformed with",trans,")")
  }
  else {
    trans.msg = ""
  }
  
  if(is.null(group.var)) {
    group.var.msg = ""
  }
  else {
    group.var.msg = paste("split by",group.var)
  }
  
  report.section = report$add.header(paste("Plots of",
                                           vars.descr,
                                           trans.msg,
                                           "as a function of",
                                           x.var,
                                           group.var.msg),
                                     section.action="push")
  
  id.vars = c(x.var,group.var)
  dat = cbind(m_a$attr[,id.vars,drop=F],count)
  dat = melt.abund.meta(dat,id.vars=id.vars,attr.names=id.vars,value.name=value.name)
  smooth_method = "loess" #"lm"
  for(clade.name in clade.names) {
    pl = ggplot(dat[dat$clade==clade.name,], aes_string(x=x.var, y=value.name,color=group.var)) +
      geom_point() +
      #geom_line(alpha=0.3, linetype=3) + 
      #geom_smooth(aes(group=group,color=group), method='lm', formula=y~x+I(x^2)+I(x^3)) + 
      stat_smooth(method=smooth_method, se = T,degree=1,size=1)
    #scale_x_date() +
    #labs(title=title)+
    #facet_wrap(~clade,scales="free")
    report$add(pl,
               caption=paste("Value",
                             trans.msg,
                             "of",
                             clade.name,
                             "as a function of",
                             x.var,
                             group.var.msg)
    )
  }
}


##Code adapted from Jyoti Shankar. Find best alpha through cross-validation as
##per the recipe from cv.glmnet help page.
##Extra argument q=dfmax+1 to make this function compatible with 
##the parameter list for stabs::stabsel (see how this is used in
##stabsel.report).
cv.glmnet.alpha <- function(y, x, family, q=NULL, seed=NULL, standardize=T,...) {
  # This seed is taken from LASSO's example. 
  if (!is.null(seed) ) { set.seed(seed) }
  # Setting the number of folds to the number of samples (leave one out)
  #is not recommended by cv.glmnet help page
  numfolds <- min(15,dim(x)[1])
  # Grid for alpha crossvalidation
  alphas <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999,1.0)
  # In a k fold crossvalidation, fold ID gives the iteration in which the sample would be in the test set.
  
  foldid <- sample(rep(seq(numfolds),length=dim(x)[1]))
  # Go through alpha grid
  # Run crossvalidation for lambda.
  # Each model for each alpha is run by a parallel core and put into a list called lassomodels
  if(!is.null(q)) {
    dfmax = q - 1
  }
  else {
    dfmax = ncol(x) + 1
  }
  lassomodels <- foreach(i = c(1:length(alphas)),.packages=c("glmnet")) %dopar% {
    #re-import required for windows parallelism with doSNOW
    library(glmnet)
    # set.seed(seed)
    # the function finds the best lambda for a given alpha
    # within each model there is cross-validation happening for lambda for each alpha.
    # lambda1 = lambda*alpha 
    # lambda2 = lambda*(1-alpha)
    model <- try({cv.glmnet(x=x, y=y, family=family,
                            nfolds=numfolds, 
                            type.measure="deviance", 
                            foldid=foldid,
                            standardize=standardize, 
                            alpha=alphas[i],
                            dfmax=dfmax,
                            ...)})
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
  
  if (best_alpha_index != 0) {
    # print the lassomodel at the best_alpha_index
    lasso_model <- lassomodels[[best_alpha_index]]
    alpha <- alphas[best_alpha_index]
    # Use lambda which gives the lowest cross validated error
    lambda <- lasso_model$lambda.min
    out <- list(param=list(alpha=alpha),
                glmnet.model=lasso_model, 
                lambda=lambda
    )
  }
  else {
    out = list()
  }
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
                            do.select.samples=F,
                            stability.transform.counts="ihs") {
  
  report.section = report$add.header("Data analysis",section.action="push")
  n.batch.levels = sum(table(data$Batch) > 0)
  resp.attr = "T1D"
  stability.model.family = "binomial"
  n.adonis.perm = 4000
  
  #data = data[data$Batch %in% c(1,2,3),]
  
  ##only families with 1 sibling and 1 patient
  #tbl = with(data,table(Family,T1D))
  #fams = unique(row.names(tbl[(tbl[,"T1D"]>0 & tbl[,"Control"]>0),]))
  #data = data[data$Family %in% fams,]
  
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  
  data.norm = norm.prop(data,col_ignore=attr.names)
  
  #make.global(tbl)
  #make.global(fams)
  #make.global(data.norm)
  
  all.clades = get.clade.names(data,attr.names)
  m_a = split_count_df(data.norm,col_ignore=attr.names)
  
  make.global(m_a)
  #DEBUG
  #stopifnot(F)
  
  m = m_a$count
  
  
  
  res = list()  
  
  #m = (m_a.abs$count > 0)
  #storage.mode(m) = "integer"
  
  count = switch(stability.transform.counts,
                 boxcox=norm.boxcox(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m_a.abs$count > 0))
  
  if (do.stability) {
    
    if(do.genesel) {
      ## We need to use row-normalized counts here, because we perform Wilcox test
      ## inside, and could false significance due to different depth of sequencing
      res.stab_sel_genesel = stab_sel_genesel(m_a$count,m_a$attr[,resp.attr])
      res$stability.genesel = res.stab_sel_genesel
      report$add.header("GeneSelector stability ranking")
      report$add.package.citation("GeneSelector")
      report$add.descr("Univariate test (Wilcoxon) is applied to each clade on random
                       subsamples of the data. Consensus ranking is found with a
                       Monte Carlo procedure. The top clades of the consensus ranking
                       are returned, along with the p-values computed on the full
                       original dataset (with multiple testing correction).")
      report$add.table(res.stab_sel_genesel$stab_feat,
                       caption=
                         paste("GeneSelector stability ranking for response ",resp.attr)
      )
      report$add.package.citation("vegan")
      report$add(
        plot.features.mds(m_a,species.sel=(colnames(m_a$count) %in% 
                                             res.stab_sel_genesel$stab_feat$name),
                          sample.group=m_a$attr[,resp.attr]),
        caption="metaMDS plot of abundance profile. 'x' marks top selected clades, 'o' marks samples"
      )
    }
    if(do.glmnet) {
      standardize.glm = T
      if(standardize.glm) {
        count = decostand(count,method="standardize",MARGIN=2)
      }
      
      cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
      registerDoSNOW(cl)  
      cv.res = cv.glmnet.alpha(m_a$attr[,resp.attr],count,family=stability.model.family,standardize=F)
      stopCluster(cl)
      
      penalty.alpha = cv.res$alpha
      #alpha = 0.8
      stab.res.c060 = stabpath(m_a$attr[,resp.attr],count,weakness=0.8,
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
      #                     resp.attr,
      #                     ")"
      #                     )
      #)
      report$add.header(paste(
        "Glmnet stability path analysis for response (",
        resp.attr,
        ")"
      )
      )
      report$add.package.citation("c060")
      report$add.descr(paste("This multivariate feature selection method builds glmnet models 
                             for a given response variable at 
                             varying strengths of L1 norm
                             regularization parameter and on multiple random subsamples at each
                             level of the parameter. The features (e.g. taxonomic clades)
                             are ranked according to their overall probability to be included in the
                             model over the range of regularization parameter values.
                             Input frequency profiles are transformed with",
                             stability.transform.counts,
                             "then standardized across samples."))
      report$add.descr(paste("Alpha penalty parameter was chosen in cross-validation as",
                             penalty.alpha))
      report$add.printed(stab.feat.c060$stable,
                         caption="Features that passed FWER in stability analysis:")
      
      if(!is.null(pl.stab)) {
        report$add(pl.stab,caption="Stability path. Probability of each variable 
                   to be included in the model as a function of L1 regularization
                   strength. Paths for top ranked variables are colored. The variables
                   (if any) that passed family wide error rate (FWER) cutoff are
                   plotted as solid lines.")
      }
      
      res$glmnet.stability=stab.feat.c060
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
        #formula_rhs = paste(formula_rhs,"+TimestampMonth",sep="")
        linfct=c("T1DT1D = 0")
      }
      res$glmer = test.counts.glmer(m_a.abs,alpha=alpha,
                                    formula_rhs=formula_rhs,
                                    linfct=linfct)
      report$add.header("Binomial mixed model analysis")
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
      report$add.header("PermANOVA (adonis) analysis of taxonomic profiles")
      report$add.package.citation("vegan")
      report$add.descr(
        paste(
          "Non-parametric multivariate test for association between
          taxonomic profiles and meta-data variables. Profile is normalized
          to proportions across clades and then to range across samples.
          Distance metric is",
          adonis.dist
        )
      )
      
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
  }
  
  if( do.select.samples && do.stability && do.glmnet ) {
    
    tryCatchAndWarn({    
      species.sel=res$glmnet.stability.genesel$stab_feat$name
      species.sel=names(take_first(res$glmnet.stability.glmnet$mean.order.p,20))
      select.samples(data.norm,attr.names,
                     species.sel=species.sel,
                     sample.group.name=resp.attr)
    })
    
  }
  
  return (res)
}


test.counts.choc.Nov_2013_gran_proposal <- function(data,attr.names,label,alpha=0.05,
                                                    do.tests=T,do.stability=T,
                                                    stability.transform.counts="ihs") {
  
  n.adonis.perm = 400
  res = list()
  m_a.abs = split_count_df(data,col_ignore=attr.names)
  
  data.norm = norm.prop(data,col_ignore=attr.names)
  
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
                 boxcox=norm.boxcox(m),
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
    res$glmnet.stability=stab.feat.c060
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


deseq2.report <- function(m_a,
                          formula.rhs,
                          test.task=list(),
                          result.tasks=list(list())
) {
  report.section = report$add.header("DESeq2 tests and data normalization",section.action="push")
  report$add.package.citation("DESeq2")
  dds = as.dds.m_a(m_a=m_a,
                   formula.rhs=formula.rhs,
                   force.lib.size=T)
  dds = do.call(DESeq,
                c(list(object=dds),
                  test.task)
  )
  res.all = foreach(result.task=result.tasks) %do% {
    res = do.call(results,
                  c(list(object=dds),
                    result.task)
    )
    res = res[order(res$padj),]  
    res.descr = paste(capture.output(head(res,0))[1:2],collapse=";")
    res.df = cbind(feature = rownames(res),as.data.frame(res))
    report$add.table(as.data.frame(res.df),
                     caption=paste("DESeq2 results for task:",
                                   formula.rhs,
                                   arg.list.as.str(result.task),
                                   res.descr,
                                   sep=";"),
                     show.row.names=F)
    res
  }
  return (list(dds=dds,results=res.all))
}

glmnet.stabpath.c060.report <- function(m_a,
                                        resp.attr,
                                        family="gaussian",
                                        steps=600,
                                        weakness=0.8,
                                        standardize.count=T,
                                        transform.count="ident",
                                        pred.attr=NULL,
                                        fwer.alpha=0.05) {
  
  report$add.header(paste(
    "Glmnet stability path analysis for response (",
    resp.attr,
    ")"
  )
  )
  
  m = m_a$count
  count = switch(transform.count,
                 boxcox=norm.boxcox(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m > 0),
                 clr=norm.clr(m))
  
  
  attr = m_a$attr
  if(standardize.count) {
    count = decostand(count,method="standardize",MARGIN=2)
  }
  
  cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
  registerDoSNOW(cl)  
  cv.res = cv.glmnet.alpha(attr[,resp.attr],count,family=family,standardize=F)
  stopCluster(cl)
  
  penalty.alpha = cv.res$alpha
  #alpha = 0.8
  stab.res.c060 = stabpath(attr[,resp.attr],count,weakness=weakness,
                           family=family,steps=steps,
                           alpha=penalty.alpha,standardize=F)
  #stab.res.c060 = stabpath(m_a$attr$A1C[m_a$attr$T1D=="T1D"],
  #                               count[m_a$attr$T1D=="T1D",],
  #                               weakness=0.9,
  #                               family="gaussian",steps=600,
  #                               alpha=penalty.alpha,standardize=F)
  
  fwer = fwer.alpha
  pi_thr = 0.6
  stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
  #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
  stab.feat.c060$penalty.alpha = penalty.alpha
  
  pl.stab = tryCatchAndWarn({
    #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
    p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
    #ggsave(stab.path.file)
    p
  })
  
  #report$add.printed(stab.res.c060$fit,
  #                   caption=paste(
  #                     "Glmnet stability path analysis for response (",
  #                     resp.attr,
  #                     ")"
  #                     )
  #)
  
  report$add.printed(summary(m_a$attr[,resp.attr]),
                     caption="Summary of response variable")
  
  report$add.package.citation("c060")
  report$add.descr("This multivariate feature selection method implements 
                  stability selection procedure by Meinshausen and Buehlmann (2010) 
                  The features (e.g. taxonomic clades)
                   are ranked according to their probability to be selected
                   by models built on multiple random subsamples of the input dataset.")
  
  report$add.descr(paste("Alpha penalty parameter was chosen in cross-validation as",
                         penalty.alpha))
  
  report$add.vector(stab.feat.c060$stable,
                    caption="Features that passed FWER in stability analysis")
  
  if(!is.null(pl.stab)) {
    report$add(pl.stab,caption="Stability path. Probability of each variable 
               to be included in the model as a function of L1 regularization
               strength. Paths for top ranked variables are colored. The variables
               (if any) that passed family wide error rate (FWER) cutoff are
               plotted as solid lines.")
  }
  
  return (stab.feat.c060)
}


stabsel.report <- function(m_a,
                           resp.attr,
                           fitfun="glmnet.lasso",
                           parfitfun="cv.glmnet.alpha",
                           args.fitfun = list(
                             family="gaussian",
                             standardize=T                                        
                           ),
                           args.stabsel = list(
                             PFER=0.05,
                             sampling.type="SS",
                             assumption="r-concave"
                           )
) {
  
  require(stabs)
  report$add.header(paste(
    "Stability selection analysis for response (",
    resp.attr,
    ")"
  )
  )
  
  resp = m_a$attr[,resp.attr]
  count = m_a$count
  
  #MB's paper
  q = ceiling(sqrt(0.8*nrow(count)))
  
  param.fit = do.call(parfitfun,
                      c(list(x=count,y=resp,q=q),
                        args.fitfun
                      )
  )$param
  
  args.fitfun = c(args.fitfun,param.fit)
  args.stabsel$q = q
  stab.res = do.call(stabsel,c(
    list(x=count,y=resp,
         fitfun=fitfun,
         args.fitfun=args.fitfun),
    args.stabsel
  )
  )
  
  report$add.printed(summary(resp),
                     caption=paste("Summary of response variable",resp.attr))
  
  report$add.package.citation("stabs")
  report$add.descr("This multivariate feature selection method implements 
                  stability selection procedure by Meinshausen and Buehlmann (2010) 
                  and the improved error bounds by Shah and Samworth (2013). 
                  The features (e.g. taxonomic clades)
                   are ranked according to their probability to be selected
                   by models built on multiple random subsamples of the input dataset.")
  
  report$add.descr(paste("Base selection method parameters that were chosen based on
                         cross-validation are:",
                         arg.list.as.str(param.fit)))
  
  report$add.descr(paste("All base selection method parameters are:",
                         arg.list.as.str(args.fitfun)))  
  
  report$add.descr(paste("Stability selection method parameters are:",
                         arg.list.as.str(args.stabsel)))
  
  #report$add.printed(stab.res)
  
  cutoff.descr = sprintf("Probability cutoff=%s corresponds to per family error rate PFER=%s",
                         stab.res$cutoff,
                         stab.res$PFER)  
  
  report$add.vector(get.feature.ranking(stab.res,only.names=F)$ranked,
                    caption=paste("Selection probability for the variables.",
                                  cutoff.descr)
  )
  
  report$add(plot(stab.res, main = NULL, type = "maxsel",
                  ymargin = 16, np = 20),
             caption=paste("Selection probability for the top ranked variables.",
                           cutoff.descr,
                           "(vertical line).")
  )
  
  return (stab.res)
}


genesel.stability.report <- function(m_a,resp.attr) {
  
  ## m_a$count passed here should be normalized for library size, because we perform Wilcox test
  ## inside, and could find false significance due to different depth of sequencing
  
  report$add.header("GeneSelector stability ranking")
  report$add.package.citation("GeneSelector")
  report$add.descr("Univariate test (Wilcoxon rank-sum) is applied to each feature (clade, gene)
                   on random
                   subsamples of the data. Consensus ranking is found with a
                   Monte Carlo procedure. The top clades of the consensus ranking
                   are returned, along with the p-values, statistic and effect size 
                  computed on the full
                   original dataset. P-values are reported with and without the 
                   multiple testing correction Benjamini & Hochberg. The effect sizes
                   for Wilcoxon rank-sum test are reported as: common-language effect
                   size (proportion of pairs where observations from the first group
                   are larger than observations from the second group; no effect
                   corresponds to 0.5); rank-biserial
                   correlation (common language effect size minus its complement; no
                   effect corresponds to 0; range is [-1;1]) and
                   absolute value of r (as defined in Cohen, J. (1988). Statistical power 
                   analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Erlbaum.)
                   ")
  
  
  res.stab_sel_genesel = stab_sel_genesel(m_a$count,m_a$attr[,resp.attr])
  report$add.printed(summary(m_a$attr[,resp.attr]),
                     caption=paste("Summary of response variable",resp.attr))
  
  report$add.table(res.stab_sel_genesel$stab_feat,
                   caption=
                     paste("GeneSelector stability ranking for response ",resp.attr))
  
  report$add.package.citation("vegan")
  m_a.mds = m_a
  m_a.mds$count = m_a.mds$count[,res.stab_sel_genesel$stab_feat$name]
  report$add(
    plot.features.mds(m_a.mds,sample.group=m_a$attr[,resp.attr]),
    caption="metaMDS plot of clades selected by GeneSelector. 'x' marks clades, 'o' marks samples"
  )
  
  return (res.stab_sel_genesel)
  
}

test.counts.glmer.report <- function(m_a,
                                     tasks,
                                     alpha=0.05) {
  descr.tpl = "
  Mixed model analysis of count data.
  The binomial family
  is used to build a set of univariate models, with each
  model describing the counts for one clade.
  We add a random effect for each sample to account
  for the overdispersion;
  P-values are estimated from the model under a null hypothesis
  of zero coefficients and a two-sided alternative. 
  Benjamini & Hochberg (1995) method is used 
  for multiple testing correction, and the significant clades
  are reported.
  %s
  "
  
  report$add.header("Binomial mixed model analysis")
  report$add.package.citation("lme4")
  res = lapply(tasks,function(task) {
    with(task,{
      res.glmer = test.counts.glmer(m_a,alpha=alpha,
                                    formula_rhs=formula.rhs,
                                    linfct=linfct)
      
      descr = sprintf(descr.tpl,descr.extra)
      report$add.descr(descr)
      report.counts.glmer(report,res.glmer)
      res.glmer
    })
  })
  return(res)
}

test.counts.adonis.report <- function(m_a,
                                      tasks,
                                      n.perm=4000,
                                      dist.metr="bray",
                                      col.trans="range",
                                      data.descr="proportions of counts") {
  ##Negative values break bray-curtis and jaccard distances; we standardize to "range" to reduce
  ##the influence of very abundant species:
  count = m_a$count
  is.dist = inherits(count,"dist")
  if(!is.dist && !is.null(col.trans) && col.trans != "ident") {
    count = decostand(count,method=col.trans,MARGIN=2)
    col.trans.descr = sprintf(" Profile columns are normalized with %s method of decostand function.",col.trans)
  }
  else {
    col.trans.descr = ""
  }
  
  if(!is.dist) {
    dist.metr.descr = sprintf(" Dissimilarity index is %s.",dist.metr)
  }
  else {
    dist.metr.descr = " Using supplied distance matrix."
  }
  
  #print(ad.res)
  report$add.header(paste("PermANOVA (adonis) analysis of ",data.descr))
  report$add.package.citation("vegan")
  report$add.descr(sprintf("Non-parametric multivariate test for association between
                           %s and meta-data variables.%s%s",
                           data.descr,
                           col.trans.descr,
                           dist.metr.descr))
  
  res = lapply(tasks,function(task) {
    with(task,{
      formula_str = paste("count",formula.rhs,sep="~")
      ad.res = adonis(
        as.formula(formula_str),
        data=m_a$attr,
        strata=if(!is.null(strata)) m_a$attr[,strata] else NULL,
        permutations=n.perm,
        method=dist.metr)
      
      report$add.printed(ad.res,caption=paste(descr,
                                              "with formula",
                                              pandoc.escape.special(formula_str)))
      report$add.table(ad.res$aov.tab,
                       show.row.names=T,
                       caption=paste(descr,"AOV Table"))
      ad.res
    })
  })
  
  return (res)
  
}

select.samples.report <- function(m_a,
                                  feature.ranking,
                                  resp.attr,
                                  n.species,
                                  n.samples) {
  
  species.sel=take_first(feature.ranking,
                         n.species)
  select.samples(m_a=m_a,
                 species.sel=species.sel,
                 sample.group.name=resp.attr,
                 n.select=n.samples)
  
}


plot.profiles.abund <- function(m_a,
                                label,
                                plot.profiles.task,                                
                                res.tests=NULL,
                                norm.count.task=NULL) {
  
  m_a.norm <- norm.count.report(m_a,
                                res.tests=res.tests,
                                descr="abundance plots",
                                norm.count.task=norm.count.task)
  
  tryCatchAndWarn({
    
    feature.order = list(list(ord=NULL,ord_descr="average abundance",sfx="ab"))
    if (!is.null(names(res.tests)) && !is.null(res.tests$stabsel)) {
      feature.order[[2]] = list(ord=get.feature.ranking(res.tests,"stabsel")$ranked,
                                ord_descr="Stability selection",
                                sfx="stabsel")
    }
    
    do.call(plot.profiles,
            c(list(m_a=m_a.norm,
                   feature.order=feature.order),
              plot.profiles.task
            )
    )
    
  })
  
}

norm.count.report <- function(m_a,norm.count.task=NULL,res.tests=NULL,descr=NULL) {
  
  if(is.null(descr)) {
    descr = ""
  }
  else {
    descr = paste("for",descr)
  }
  
  if(!is.null(norm.count.task)) {
    
    descr = paste("Count normalization method",descr,":",arg.list.as.str(norm.count.task))
    
    ## if dds is required but not defined (set to NA)
    if(!is.null(norm.count.task$method.args$dds) && 
         is.na(norm.count.task$method.args$dds)) {
      norm.count.task$method.args$dds = res.tests$deseq2$dds
    }
    
    m_a.norm <- do.call(norm.count.m_a,
                        c(list(m_a),
                          norm.count.task)
    )
    
  }
  else {
    
    descr = paste("Counts are not normalized",descr)
    m_a.norm = m_a
    
  }
  
  report$add.descr(descr)
  
  return(m_a.norm)
}

test.counts.project <- function(m_a,
                                label,
                                do.genesel=T,
                                do.stabsel=T,
                                do.glmer=T,
                                do.adonis=T,
                                do.select.samples=F,
                                do.deseq2=T,
                                do.divrich=T,
                                do.plot.profiles.abund=T,
                                do.heatmap.abund=T,
                                count.filter.feature.options=NULL,
                                norm.count.task=NULL,
                                stabsel.task=NULL,
                                genesel.task=NULL,
                                adonis.task=NULL,
                                glmer.task=NULL,
                                select.samples.task=NULL,
                                deseq2.task=NULL,
                                divrich.task=NULL,
                                plot.profiles.task=NULL,
                                plot.profiles.abund.task=NULL,
                                heatmap.abund.task=NULL,
                                alpha=0.05,
                                do.return.data=T
) {
  
  report.section = report$add.header("Data analysis",section.action="push")
  
  res = new_mgsatres()
  
  if(do.divrich) {
    if(is.null(divrich.task$beta.task$adonis.task)) {
      if(!is.null(divrich.task$beta.task)) 
      { divrich.task$beta.task$adonis.task = adonis.task }
      if(!do.adonis) {
        divrich.task$beta.task$adonis.task = NULL
      }
    }
    tryCatchAndWarn({ 
      res$divrich <- do.call(mgsat.divrich.report,
                             c(list(m_a,
                                    plot.profiles.task=plot.profiles.task),
                               divrich.task)
      )
    })
  }
  
  if(!is.null(count.filter.feature.options)) {
    m_a = count.filter.report(m_a,
                              count.filter.options=count.filter.feature.options,
                              "feature")
  }
  
  if(do.deseq2) {
    tryCatchAndWarn({ 
      res$deseq2 = do.call(deseq2.report,c(list(m_a=m_a),deseq2.task))
    })
  }
  
  ## this is done after an optional call to deseq2 in case the norm.method
  ## wants deseq2 normalization
  
  m_a.norm <- norm.count.report(m_a,
                                res.tests=res,
                                descr="data analysis",
                                norm.count.task)
  
  make.global(m_a.norm)
  make.global(m_a)
  
  if(do.genesel) {
    tryCatchAndWarn({ 
      res$genesel = do.call(genesel.stability.report,c(list(m_a=m_a.norm),genesel.task))
    })
  }
  
  if(do.stabsel) {
    tryCatchAndWarn({ 
      res$stabsel = do.call(stabsel.report,c(list(m_a=m_a.norm),
                                             stabsel.task
      )
      )
    })
  }
  
  if(do.glmer) {  
    tryCatchAndWarn({ 
      res$glmer = do.call(test.counts.glmer.report,
                          c(
                            list(m_a=m_a,
                                 alpha=alpha),
                            glmer.task
                          )
      )
    })
  }
  
  if(do.adonis) {
    tryCatchAndWarn({     
      res$adonis = do.call(test.counts.adonis.report,
                           c(
                             list(m_a=m_a.norm),
                             adonis.task
                           )
      )
    })
  }
  
  if( do.select.samples && do.stabsel ) {
    
    tryCatchAndWarn({ 
      do.call(select.samples.report,
              c(
                list(m_a=m_a.norm,
                     feature.ranking=get.feature.ranking(res)$ranked
                ),
                select.samples.task
              )
      )
    })
    
  }
  
  if( do.plot.profiles.abund ) {
    
    tryCatchAndWarn({ 
      do.call(plot.profiles.abund,
              c(
                list(m_a=m_a,
                     label=label,
                     res.tests=res,
                     plot.profiles.task=plot.profiles.task),
                plot.profiles.abund.task
              )
      )
    })
  }
  
  if (do.heatmap.abund) {
    
    tryCatchAndWarn({
      do.call(heatmap.counts,
              c(list(m_a=m_a.norm),
                heatmap.abund.task
              )
      )
      
    })
  }
  
  
  if(do.return.data) {
    res$m_a = m_a
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

library(Heatplus)
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
heatmap.counts <- function(m_a,attr.annot.names,
                           attr.row.labels=NULL,
                           caption="Heatmap",
                           max.species.show=30, 
                           stand.clust=NULL,
                           trans.clust=NULL,
                           dist.metr="euclidean",
                           stand.show="range",
                           trans.show=NULL,
                           attr.order=NULL,
                           agglo.fun.order=sum,
                           cluster.row.cuth=2) {
  
  library(RColorBrewer)
  library(Heatplus)
  library(vegan)
  
  ##permute samples to make sure that our dendrogram
  ##clustering is not influenced by the original order
  perm.ind = sample(nrow(m_a$count))
  
  count.src = m_a$count[perm.ind,]
  count = count.src
  
  if(!is.null(trans.clust)) {
    count = do.call(trans.clust,list(count))
  }
  if(!is.null(stand.clust)) {
    count = decostand(count,method=stand.clust,MARGIN=2)
  }
  
  attr = m_a$attr[perm.ind,]
  data.dist.samp <- vegdist(count, method = dist.metr)
  row.clus <- hclust(data.dist.samp, "ward.D2")
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
  count.sub = count[,1:min(ncol(count),max.species.show)]
  # you have to transpose the dataset to get the taxa as rows
  data.dist.taxa <- vegdist(t(count.sub), method = dist.metr)
  col.clus <- hclust(data.dist.taxa, "ward.D2")
  
  ## go back to un-normalized
  count.sub = count.src[,1:min(ncol(count),max.species.show)]
  if(!is.null(trans.show)) {
    count.sub = do.call(trans.show,list(count.sub))
  }
  if(!is.null(stand.show)) {
    count.sub = decostand(count.sub,method=stand.show,MARGIN=2)
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
  report$add(plot(pl.heat),caption=caption)
}

heatmap.t1d <- function(meta.data,attr.names,caption="Heatmap") {
  report.section = report$add.header("Heatmap of abundance profile",section.action="push")  
  caption.ini = caption
  meta.data.ini = meta.data
  for(i in 1:6) {
    meta.data = meta.data.ini
    if(i>1) {
      meta.data$TimestampMonth = sample(meta.data$TimestampMonth)
      caption = paste(caption.ini,"randomization ",i,"of TimestampMonth")
    }
    heatmap.counts(meta.data,
                   attr.names,
                   attr.annot.names=c("T1D","TimestampMonth"),
                   attr.row.labels=NULL,
                   caption=caption,
                   stand.show="max",
                   trans.show=sqrt,
                   attr.order="TimestampMonth",
                   agglo.fun.order=mean)
    if(i==1) {
      report.section = report$add.header("Heatmaps with randomly permuted meta-data variable")
      report$add.descr("Compare these with the heatmap of non-randomized data")
    }
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

load.meta.t1d.sebastian_format <- function(file_name,batch=NULL,aggr_var=NULL) {
  
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
  meta_qa = meta[meta$Sample.QA=="PASS",]
  
  isBioRepeatFirst = duplicated(meta_qa$SubjectID) & ! duplicated(meta_qa$BioSampleID)
  
  BioSampleIDRep = meta_qa$BioSampleID[isBioRepeatFirst]
  
  ## mark all failed samples or bio repeats following at least one good sample
  meta$isBioRepeatOrFailed = (meta$BioSampleID %in% BioSampleIDRep) | 
    ! (meta$Sample.QA=="PASS")
  
  SampleIDRep = meta_qa$SampleID[duplicated(meta_qa$SubjectID)]
  
  ## mark all failed samples or repeats following at least one good sample
  meta$isRepeatOrFailed = (meta$SampleID %in% SampleIDRep) | 
    ! (meta$Sample.QA=="PASS")
  
  
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
              ! (meta$Sample.QA=="PASS")) #should not see otherwise
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
      (meta$Sample.QA=="PASS") & 
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
    
    meta = aggregate.by.meta.data(meta,
                                  aggr_var,
                                  names(meta))  
  }
  
  return (meta)
}



read.t1d <- function(taxa.level="3",batch=NULL,aggr_var_meta = NULL) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("43aa6.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary.txt")
  ##Sebastian's run:
  #moth.taxa <- read.mothur.taxa.summary("ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("acc5a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  ##Andrey's run
  #taxa.file = "dec1a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  ##1st MiSeq run
  #taxa.file = "147936b4aacf1a450a0dfe2d97c3dd5a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  ##Verification MiSeq run
  #taxa.file = "2014-05_final_53ea024a8bdf736a4dc90f4ae6d02d6a.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  #otu.shared.file="otu/6ce0d3fda8e528db61a72ac69cf5592c.files_x1.sorted.0.03.shared"
  #cons.taxonomy.file="otu/8478e43c54e283847f7de06e9e14267f.files_x1.sorted.0.03.cons.taxonomy"
  #meta.file="aliq_id_to_metadata_for_yap_inp_t1d_v1_v3_2014-04-11_fixed_diagnosis_dates.tsv"
  ##All samples up to 2014-09-22
  taxa.file = "a93eeedeef2878be17d30f27b1b0de1c.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
  otu.shared.file="d9f44114ac6369cc66978002228bc5d3.files_x1.sorted.0.03.shared"
  cons.taxonomy.file="4272870500ad07a7270f9772581fe7fe.files_x1.sorted.0.03.cons.taxonomy"
  meta.file="aliq_id_to_metadata_for_T1D_YAP_run_20140922.tsv"
  if (taxa.level == "otu") {
    taxa.lev.all = read.mothur.otu.with.taxa(otu.shared.file=otu.shared.file,cons.taxonomy.file=cons.taxonomy.file)
  }
  else {
    moth.taxa <- read.mothur.taxa.summary(taxa.file)
    ##Vishal's run
    #moth.taxa <- read.mothur.taxa.summary("vishal_all_e61bb.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
    #moth.taxa <- read.mothur.taxa.summary("vishal_aa_isdef_3b4ea.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
    taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
    make.global(taxa.lev.all)
  }
  #taxa.lev = count.filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
  taxa.lev = taxa.lev.all
  make.global(taxa.lev)
  meta = load.meta.t1d(meta.file,batch=batch,aggr_var = aggr_var_meta)
  make.global(meta)
  #meta = load.meta.t1d("annotation20130819.csv")
  #meta = load.meta.t1d("annotation20130819.csv","all_batches_cleaned_yap_input.csv")
  report$add.p(sprintf("Loaded %i records from taxonomic assignment file %s for taxonomic level %s",
                       nrow(taxa.lev),taxa.file,taxa.level))
  taxa.meta = merge.counts.with.meta(taxa.lev,meta)
  report$add.p(sprintf("After merging with metadata, %i records left",
                       nrow(taxa.meta)))
  return(taxa.meta)
}

plot.features.mds <- function(m_a,species.sel=NULL,sample.group=NULL,show.samples=T,show.species=T) {
  ##https://stat.ethz.ch/pipermail/r-help/2008-April/159351.html
  ##http://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
  
  m = m_a$count
  ##default distance is bray
  m = decostand(m,method="range",MARGIN=2)
  sol = metaMDS(m,autotransform = F,trymax=40)
  if(show.samples) {
    site.sc <- scores(sol, display = "sites")
    plot(site.sc)
    unique_sample.group = unique(sample.group)
    points(sol,display="sites",col=sample.group)
    if(!is.null(unique_sample.group)) {
      legend(-0.5,1,unique_sample.group,col=1:length(sample.group),pch=1)  
    }
  }
  if(show.species) {
    species.sc <- scores(sol, display = "species")
    if(!is.null(species.sel)) {
      species.sel.mask=(colnames(m_a$count) %in% species.sel)
    }
    else {
      species.sel.mask=rep(T,length(colnames(m_a$count)))
    }
    points(sol,display="species",pch="x",select=species.sel.mask)
    text(sol,display="species", cex=0.7, col="blue",select=species.sel.mask)
  }
  #points(site.sc,col=m_a$attr$T1D)
  #points(species.sc)
  #points(species.sc["Streptococcus_0.1.11.1.2.6.2",1:2,drop=F],pch="+")
  
}

cut.top.predictions<-function(scores,labels,sample.ids,n.cut,filter.by.label=F) {
  lab.levels = levels(labels)
  stopifnot(length(lab.levels)==2)
  n = length(scores)
  stopifnot(length(labels)==n && length(sample.ids)==n)
  
  n.cut = min(n.cut,round(n/2))
  
  res = list()
  
  if (filter.by.label) {
    
    lab = lab.levels[1]
    ord = order(scores)
    cum.cnt = cumsum(labels[ord]==lab)
    i.ord = first(cum.cnt>=n.cut)
    stopifnot(length(i.ord)>0) #index found
    cut = scores[ord[i.ord]]
    ids = sample.ids[(scores <= cut) & (labels == lab)]
    res[[lab]] = list(cut=cut,ids=ids)
    
    lab = lab.levels[2]
    ord = order(-scores)
    cum.cnt = cumsum(labels[ord]==lab)
    i.ord = first(cum.cnt>=n.cut)
    stopifnot(length(i.ord)>0) #index found
    cut = scores[ord[i.ord]]
    ids = sample.ids[(scores >= cut) & (labels == lab)]
    res[[lab]] = list(cut=cut,ids=ids)
    
  }
  else {
    ord = order(scores)
    cut = scores[ord[n.cut]]
    ids = sample.ids[ord[1:n.cut]]
    res[[lab.levels[1]]] = list(cut=cut,ids=ids)
    
    ord = order(-scores)
    cut = scores[ord[n.cut]]
    ids = sample.ids[ord[1:n.cut]]
    res[[lab.levels[2]]] = list(cut=cut,ids=ids)
  }
  return (res)
  
}

show.pred.perf <- function(pred,pred.score,labels) {
  library(ROCR)
  report.section = report$add.header("Prediction performance measures",section.action="push")  
  report$add.package.citation("ROCR")
  #print(table(labels,pred.score>0))
  pred.perf = prediction(pred.score,labels)
  report$add.table(table(labels,pred),caption="Confusion table")
  # Plot ROC curve
  perf <- performance(pred.perf, measure = "tpr", x.measure = "fpr")
  report$add(plot(perf),caption="ROC curve")
  # Plot precision/recall curve
  perf <- performance(pred.perf, measure = "prec", x.measure = "rec")
  report$add(plot(perf),caption="Precision/recall")
  perf <- performance(pred.perf, measure = "acc")
  report$add(plot(perf),caption="Accuracy")
}

select.samples <- function(m_a,
                           species.sel,
                           sample.group.name,
                           n.select=20,
                           filter.by.label=T,
                           selection.file=NULL) {
  library(kernlab)
  library(caret)
  
  report.section = report$add.header("Selecting most different samples with regard to a phenotype",section.action="push")
  report$add.package.citation("kernlab")
  report$add.package.citation("caret")
  report$add.descr(paste("This procedure selects those samples that are most different with regard to a grouping variable"
                         ,sample.group.name,". This can be used to select a subset for WGS or 
                         transcriptomics sequencing based on 16S profiles. Support Vector Machine is built using previously selected clades,
                         and",n.select,"samples corresponding to each of the two levels
                         of the grouping variable are picked. Samples are picked as 
                         predicted correctly by the linear SVM after applying to the frequency 
                         profiles the 
                         inverse hyperbolic sign transform and normalizing across columns, 
                         and selecting samples that are maximally
                         distant from the SVM separating plane. The nuisanse parameter
                         of the SVM is picked through a grid search maximizing
                         prediction accuracy in training with resampling.
                         Accuracy of the final model is reported both on the
                         training set and with cross-validation. Several
                         performance charts are shown for the training set.
                         After selecting samples, the same metrics are shown
                         for the selected samples only. Additionally, the
                         abundance profiles are plotted for the selected samples.
                         Typically, you might expect a fairly poor
                         accuracy for the full sample set, and very good - for the
                         selected samples assuming that the number of selected
                         samples is a small fraction of the total. For reference, a purely
                         random binary classifier would result in a 50% accuracy."))
  
  sample.group = m_a$attr[,sample.group.name]
  
  ## nasty "feature" - when species.sel is a factor,
  ## the subscripting of count seems to use the integer level of a factor
  species.sel = as.character(species.sel)
  
  m = m_a$count[,species.sel]
  
  #If building kernel from a distance:
  #S = np.exp(-D * gamma), where one heuristic for choosing gamma is 1 / num_features
  #or
  #S = 1. / (D / np.max(D))
  #m = norm.boxcox(m)
  
  m = ihs(m)
  
  m = decostand(m,method="standardize",MARGIN=2)
  
  report$add.vector(colnames(m),name="Clade",caption="Using these clades for sample selection.")
  report$add.header("Models trained on the full dataset and their performance")  
  #make.global(species.sel)
  #make.global(m)
  #make.global(sample.group)
  bootControl <- trainControl(number = 40)
  #set.seed(2)
  scaled = F
  sigma = sigest(m,scaled=scaled)[2]
  ##print(paste("sigma",sigma))
  ##using prob.mod=T screws up the model - the plane is no longer at 0 decision
  ##threshold. Maybe bias is lost?
  if(T) {
    cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
    registerDoSNOW(cl)      
    #class.weights=1/table(sample.group),
    #tuneGrid = expand.grid(sigma=sigma,C = 4**(-2:7)),
    mod.fit <- train(m, sample.group,
                     method = "svmLinear", #"svmRadial" 
                     prob.model = F,
                     cross=3,
                     class.weights=1/table(sample.group),
                     tuneGrid = expand.grid(C = 4**(-2:7)),
                     trControl = bootControl, scaled = scaled,
                     metric="Accuracy")  #"Kappa"
    stopCluster(cl)
    report$add.printed(mod.fit,caption="Results of SVM parameter fitting")
    mod = mod.fit$finalModel
  }
  else {
    for (C in c(1)) {
      mod = ksvm(x=m, y=sample.group, kernel = "vanilladot",
                 kpar = "automatic", C = C, cross = 6, prob.model = F,
                 class.weights=1/table(sample.group),scaled=scaled)
      #print(cross(mod))
      print(mod)
    }
  }
  
  report$add.printed(mod,caption="Best model trained on the full input set")  
  
  #mod.pred.prob = predict(mod, m, type = "probabilities")
  #make.global(mod.pred.prob)
  #print(sum(mod.pred.prob[,"Control"]>0.5))
  mod.pred.score = predict(mod,m,type="decision")
  mod.pred = predict(mod,m)
  #make.global(mod.pred.score)
  
  report$push.section(report.section)
  show.pred.perf(mod.pred,mod.pred.score,sample.group)  
  report$pop.section()
  
  #}
  report$add.header("Performance of the models and plots for selected samples")  
  sample.ids = rownames(m)
  cut.res = cut.top.predictions(mod.pred.score,
                                sample.group,
                                sample.ids,
                                n.cut=n.select,
                                filter.by.label=filter.by.label)
  ids.sel = c(laply(cut.res,function(x) x[["ids"]]))
  mask.sel = sample.ids %in% ids.sel
  ids.sel = sample.ids[mask.sel]
  pred.score.sel = mod.pred.score[mask.sel]
  pred.sel = mod.pred[mask.sel]
  sample.group.sel = sample.group[mask.sel]
  samples.sel = data.frame(SampleID=ids.sel,Group=sample.group.sel,Score=pred.score.sel)
  samples.sel = samples.sel[order(samples.sel$Score),]
  report$add.table(samples.sel,
                   caption="Selected samples, sorted by score. You can pick subsets at the opposite score extremes.")
  samples.all = data.frame(SampleID=sample.ids,Group=sample.group,Score=mod.pred.score)  
  samples.all = samples.all[order(samples.all$Score),]
  fn = report$make.file.name("sample.selection.tsv")
  write.table(samples.all,fn,sep="\t",row.names = FALSE)
  report$add.descr(paste(
    "All sample IDs with decision score and grouping variable are saved into file ",
    fn,"that you can use to fine tune the selection."))
  report$push.section(report.section)
  show.pred.perf(pred.sel,pred.score.sel,sample.group.sel)
  
  report$add.header("Abundance Plots")
  
  m_a.sel = m_a
  
  m_a.sel$attr = m_a.sel$attr[mask.sel,]
  
  m_a.sel$count = m_a.sel$count[mask.sel,species.sel]
  
  pl.hist = plot.abund.meta(m_a=m_a.sel,
                            id.vars=c(sample.group.name),
                            geom="bar",
                            id.var.dodge=sample.group.name
  )
  #env=as.environment(as.list(environment(), all.names=TRUE))
  #print(names(as.list(env)))
  #print(evals("pl.hist",env=env))
  report$add(pl.hist,
             caption=paste("Abundance profile of samples maximally different with regard to ",
                           sample.group.name,"(only clades that were used for selection are shown)")
  )
  
  m_a.sel$count = m_a$count[mask.sel,]
  
  pl.hist = plot.abund.meta(m_a=m_a.sel,
                            id.vars=c(sample.group.name),
                            geom="bar",
                            id.var.dodge=sample.group.name
  )
  #env=as.environment(as.list(environment(), all.names=TRUE))
  #print(names(as.list(env)))
  #print(evals("pl.hist",env=env))
  report$add(pl.hist,
             caption=paste("Abundance profile of samples maximally different with regard to ",
                           sample.group.name,"(the most abundant clades are shown, even those not used for selection)")
  )
  report$pop.section()
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
## as a function of attr metadata var (e.g. patient-control classification)
## (Binomial regression model. Add a sample specific random effect
## to account for overdispersion.)
## Value: p value of testing against a null hypothesis that 
## model coefficients are zero (ignoring intercept; two-sided)
test.counts.glmer.col <- function(taxa.count,
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
                      data=attr,
                      family="binomial",
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
    
    print("Warnings or errors in glmer, returning NA")
    
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
test.counts.glmer <- function(m_a,
                              formula_rhs,
                              linfct,
                              alpha=0.05,
                              p_adjust_method="BH") {
  ##http://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html
  library(lme4)
  library(multcomp)
  count = m_a$count
  count_rowsum = rowSums(count)
  #attr = m_a$attr[,c("T1D","FamilyID","SampleID","Batch")]
  attr = m_a$attr
  # Not working properly with .parallel=T in the presence
  # of warnings or errors - getting all NAs no matter what
  # I try for error handling in tryCatch
  ##TODO: try foreach
  p_vals = aaply(count,
                 2,
                 test.counts.glmer.col,
                 attr,
                 count_rowsum,
                 formula_rhs=formula_rhs,
                 linfct=linfct,
                 .inform=F,
                 .parallel=F,
                 .paropts=list(.packages=c("lme4","multcomp")))
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

report.counts.glmer <- function(report,x,...) {
  report$add.p("Test results from fitting mixed effects binomial model for each clade")
  report$add.p("Right-hand side of formula:")
  report$add.p(x$formula_rhs)
  report$add.p("Specification of the linear hypotheses (see glht):")
  report$add.p(x$linfct)
  
  ## our default table style of "grid" somehow screws up row representation, so we use "rmarkdown" here
  report$add.vector(x$p_vals_signif,"p_value",
                    caption="p-values that pass significance cutoff before multiple testing correction.",
                    style="rmarkdown")
  report$add.vector(x$p_vals_adj_signif,"p_value",
                    caption="p-values that pass significance cutoff after BH multiple testing correction.",
                    style="rmarkdown")
  report$add.vector(x$p_vals_adj[1:min(20,length(x$p_vals_adj))],"p_value",
                    caption="Top 20 p-values after BH multiple testing correction.",
                    style="rmarkdown")
  failed.names = names(x$p_vals[is.na(x$p_vals)])
  if(length(failed.names)>0) {
    report$add(list(failed.names),
               caption="Clades for which the model could not be built for any reason")
  }
}

## Wilcoxon rank-sum effect size r taken from Andy Fields
## two-sided
wilcox.eff.size.r <-function(pval, n){
  z<- qnorm(pval/2)
  r<- z/ sqrt(n)
  return(abs(r))
}

wilcox.eff.size <- function(stat,pval,group) {
  group = factor(group)
  taby <- table(group)
  stopifnot(nlevels(group) == 2)
  npairs = prod(taby)
  nsamp = sum(taby)
  common.lang = stat/npairs
  rbs = 2*common.lang - 1
  r = wilcox.eff.size.r(pval,nsamp)
  return (data.frame(common.lang.eff.size=common.lang,
                     rank.biserial.corr.eff.size=rbs,
                     r.eff.size=r)
  )
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
                             maxrank=20,
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
  
  y = factor(y)
  
  ranking_methods = c(RankingWilcoxon,RankingLimma,RankingFC)
  ranking_method = RankingWilcoxon
  n_feat = nrow(x)
  n_samp = ncol(x)
  #contrary to the help page, does not work when y is a factor - need to convert
  fold_matr = GenerateFoldMatrix(y = as.numeric(y), k=trunc(n_samp*(1-fold_ratio)),replicates=replicates)
  stopifnot(type=="unpaired") #effect size not implemented
  rnk = ranking_method(x,y,type=type,pvalues=T)
  rnk.vals = toplist(rnk,n_feat)
  rnk.vals[rnk.vals$index,] = rnk.vals
  rnk.vals$pval.adjusted = p.adjust(rnk.vals$pval,method="BH")
  rnk.vals = cbind(data.frame(name=rownames(x)),rnk.vals)
  
  rnk.vals = cbind(rnk.vals,
                   wilcox.eff.size(stat=rnk.vals$statistic,
                                   pval=rnk.vals$pval,
                                   group=y))
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
  #pvals are always NA somehow, and we do not need the 'index' in
  #the output, so not adding 'selected' to the output
  index = selected$index
  rnk.vals$index = NULL
  rnk.vals = rnk.vals[index,]
  return (list(stab_feat=rnk.vals,gsel=gsel))
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


pair.counts <- function(m_a,pair.attr) {
  m_m = melt(m_a$count,varnames=c("SampleID","clade"),value.name="cnt")
  attr_sub = m_a$attr[,c("SampleID",pair.attr)]
  m_m = join(m_m,attr_sub,by="SampleID",type="inner",match="first")
  
}

report.sample.count.summary <- function(m_a,meta.x.vars=c(),group.var=NULL) {
  report.section = report$get.section()
  
  m_a.summ=make.sample.summaries(m_a)
  
  report$add.vector(c(summary(m_a.summ$count[,"count.sum"])),caption="Summary of counts per sample")
  
  if(!(is.null(meta.x.vars) & is.null(group.var))) {
    
    report$add.header("Iterating over meta data variables")
    report$push.section(report.section)
    
    for(x.var in meta.x.vars) {
      
      show.sample.summaries.meta(m_a=m_a.summ,
                                 x.var=x.var,
                                 group.var=group.var)
      
    }
    
    report$pop.section()
    
  }
  
}

proc.t1d <- function() {
  #taxa.levels = c("otu")
  #taxa.levels = c("otu",2,3,4,5,6)
  taxa.levels = c(6)
  #batches = list(c(1,2,3),c(1),c(2),c(3),c(2,3),c(1,3))
  batches = list(c(4))
  do.std.plots = T
  do.tests = T
  do.heatmap = T
  do.clade.meta = T
  do.sample.summary.meta = T
  do.profile = T
  do.summary.meta = T
  
  report$add.descr("Set of analysis routines is applied
                   in nested loops over combinations of batches (if present) and 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out the batch combination and taxonomic level.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture.")
  
  report.section = report$get.section()
  report$add.header("Iterating over batch combinations")
  report$push.section(report.section)
  for (batch in batches) {
    
    report$add.header(paste(c("Batch combination:",batch),collapse=" "))
    label_batch = paste("16s","b",paste(batch,collapse="-"),sep=".")
    report$add.header("Iterating over taxonomic levels")
    report$push.section(report.section)
    for (taxa.level in taxa.levels) {
      taxa.meta = read.t1d(taxa.level,batch=batch,aggr_var_meta = "AliquotID")
      report$add.header(paste("Taxonomic level:",taxa.level))
      report$push.section(report.section)
      report$add.header("Data preparation")
      report$add.p(paste("Before filtering for QAed and redundant samples:",nrow(taxa.meta$data)))
      
      label = paste(label_batch,"l",taxa.level,sep=".",collapse=".")
      
      aggrBySubject = F
      if(!aggrBySubject) {
        taxa.meta$data = taxa.meta$data[taxa.meta$data$Sample.QA=="PASS",]
        taxa.meta.aggr = taxa.meta
        report$add.p(paste("After filtering for QAed samples:",nrow(taxa.meta$data)))      
        
        taxa.meta$data = taxa.meta$data[!taxa.meta$data$isRepeatOrFailed,] 
        taxa.meta.aggr = taxa.meta
        report$add.p(paste("After filtering for repeated samples per subject:",nrow(taxa.meta$data)))
        
      }
      else {
        taxa.meta$data = taxa.meta$data[!taxa.meta$data$isAnnualRepeatOrFailed,]
        report$add.p(paste("After filtering for QAed and non-annual repeat samples:",nrow(taxa.meta$data)))
        
        #aggr_var = "AliquotID"
        aggr_var = "SubjectID"
        
        taxa.meta.aggr = aggregate.by.meta.data(taxa.meta$data,
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
      
      taxa.meta.aggr$data = count.filter(taxa.meta.aggr$data,
                                         col_ignore=taxa.meta.aggr$attr.names,
                                         min_max=10,
                                         max_row_sum=100000)
      #taxa.meta.data = count.filter(taxa.meta.data,col_ignore=taxa.meta.attr.names,
      #                              min_mean=0.002,
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
      
      taxa.meta = taxa.meta.aggr
      
      if(do.summary.meta) {
        
        #switch it off for other iterations
        do.summary.meta = F
        
        report$add.header("Summary of metadata variables after filtering samples")
        
        m_a = split_count_df(taxa.meta$data,col_ignore=taxa.meta$attr.names)
        
        report$add.printed(summary(m_a$attr),caption="Summary of metadata variables")
        
        xtabs.formulas = list("~T1D","~T1D+FamilyID","~FamilyID","~SubjectID")
        for(xtabs.formula in xtabs.formulas) {
          fact.xtabs = xtabs(as.formula(xtabs.formula),data=taxa.meta$data,drop.unused.levels=T)
          report$add.table(fact.xtabs,show.row.names=T,caption=paste("Sample cross tabulation",xtabs.formula))
          report$add.printed(summary(fact.xtabs))
        }
        
        with(taxa.meta$data,{
          
          report$add.printed(summary(aov(age~T1D)),
                             caption="ANOVA for age and cohort")
          report$add(qplot(T1D,age,geom="violin"),
                     caption="Violin plot for age and cohort")
          
          #summary(glht(lm(age~visit,data=meta[meta$Sample.type=="patient",]),linfct="visit=0"))
          #summary(glht(lmer(age~visit+(visit|Sample.type),data=meta),linfct="visit=0"))
          report$add(ggplot(taxa.meta$data,aes(x=age,y=YearsSinceDiagnosis,color=T1D))+
                       geom_point()+
                       stat_smooth(method="loess", se = T,degree=1,size=1),
                     caption="Plot for age and time since diagnosis with Loess trend line")
          
        })
        
        with(taxa.meta$data,{
          report$add.printed(cor.test(YearsSinceDiagnosis,
                                      Timestamp,
                                      method="spearman"),
                             caption="Spearman RHO for metadata variables")
          report$add.printed(cor.test(YearsSinceDiagnosis,
                                      age,
                                      method="spearman"),
                             caption="Spearman RHO for metadata variables")          
        })
        
        report.sample.count.summary(taxa.meta)        
        
      }
      
      res.tests = NULL      
      if (do.tests) {
        
        if(taxa.level==6) {
          do.select.samples=T
        }
        else {
          do.select.samples=F
        }
        
        res.tests = tryCatchAndWarn({
          test.counts.t1d(taxa.meta.aggr$data,taxa.meta.aggr$attr.names,
                          label=label,
                          stability.transform.counts="ihs",
                          do.stability=T,
                          do.tests=T,
                          do.genesel=T,do.glmnet=T,
                          do.glmer=T,do.adonis=T,
                          do.select.samples=do.select.samples
          )}
        )
        
      }
      if (do.std.plots) {
        report$add.header("Plots")
        if(do.heatmap) {
          tryCatchAndWarn({
            heatmap.t1d(taxa.meta.aggr$data,taxa.meta.aggr$attr.names)
          })
        }
        tryCatchAndWarn({
          plot.group = list(c("T1D"))
          if(n.batches>1) {
            plot.group[[2]] = c("T1D","Batch")
          }
          std.plots(taxa.meta.aggr$data,
                    taxa.meta.aggr$attr.names,
                    id.vars.list=
                      plot.group,
                    label=label,
                    res.tests=res.tests,
                    do.clade.meta=do.clade.meta,
                    do.sample.summary.meta=do.sample.summary.meta,
                    do.profile=do.profile,
                    clade.meta.x.vars=c("YearsSinceDiagnosis","TimestampDate","age")
          )
        })
      }
      report$pop.section()
    }
    report$pop.section()
  }
  report$pop.section()
}

read.t1d.mg <-function(annot.type,level) {
  annot.dir = "../BATCH_01_02_META"
  mgrast.dir = paste(annot.dir,"BATCH_01-02_METAGENOMICS_MGRAST",sep="/")
  if (annot.type == "humann") {
    counts = read.humann.summary(paste(annot.dir,"humann/output/04b-keg-mpt-cop-nul-nve-nve.txt",sep="/"))
    counts = count.filter(counts,col_ignore=c(),min_max_frac=0,min_max=0,min_row_sum=0,other_cnt="other")      
  }
  else if (annot.type %in% c("cog","kegg","subsys")) {
    make.global(level)
    counts = read.mgrast.summary(paste(mgrast.dir,paste(annot.type,"tsv",sep="."),sep="/"),
                                 file_name.id.map=paste(mgrast.dir,"mgrast_to_samp_id.tsv",sep="/"))
    make.global(counts)
    print(level)
    counts = mgrast.to.abund.df(counts,level)
    counts = count.filter(counts,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=0,other_cnt="other")
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
        
        res.tests = tryCatchAndWarn({test.counts.t1d(count.meta.data,count.meta.attr.names,
                                                     label=label,
                                                     stability.transform.counts="ihs",
                                                     do.stability=T,
                                                     do.tests=T)})
      }
      if (do.std.plots) {
        tryCatchAndWarn({
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
        })
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
  
  data.norm = norm.prop(data,col_ignore=attr.names)
  
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
                 boxcox=norm.boxcox(m),
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
    res$glmnet.stability=stab.feat.c060
  }
  
  if (do.tests) {
    ##Negative values break bray-curtis and jaccarda= distances; we standardize to "range" to reduce
    ##the influence of very abundant species
    count = decostand(m_a$count,method="range",MARGIN=2)
    adonis.dist = "bray"
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
  taxa.lev = count.filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
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
      res.tests = tryCatchAndWarn({
        test.counts.mr_oralc(taxa.meta.data,taxa.meta.attr.names,
                             label=label,
                             stability.transform.counts="ihs",
                             do.stability=T,
                             do.tests=T)
      })
    }
    if (do.std.plots) {
      tryCatchAndWarn({
        std.plots(taxa.meta.data,taxa.meta.attr.names,id.vars.list=
                    list(
                      c("Subject","sample.type"),
                      c("sample.type"),
                      c("Subject")
                    ),
                  label=label,
                  res.tests=res.tests
        )
      })
    }
  }
}


count.summary <- function(count,fun,group) {
  group = as.data.frame(group)
  group.summ = ddply(cbind(as.data.frame(count),group),names(group),colwise(fun))
  return (group.summ)
}

group.mean.ratio <- function(count,group) {
  group.mean = count.summary(count,mean,group)
  stopifnot(dim(group.mean)[1] == 2)
  x = (group.mean[1,-1] / group.mean[2,-1])
  row.names(x) = paste(group.mean[,1],collapse=".")
  return (x)
}

read.pieper.t1d <- function(file.name="aim3/September 28 analysis_T1D.txt") {
  
  file_name = file.name
  #data =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  data =read.delim(file_name, header=T,sep="\t",stringsAsFactors=FALSE)
  #row.names(meta) = meta$ID
  id.ind = 13
  feat.attr.ind.last = 15
  feat.attr = data[,c(1:feat.attr.ind.last)]
  row.names(feat.attr) = data[,id.ind]
  row.names(data) = data[,id.ind]
  data = data[,-c(1:feat.attr.ind.last)]
  data = t(data)
  meta = as.data.frame(t(as.data.frame(strsplit(row.names(data),"[._]"))))
  row.names(meta) = row.names(data)
  names(meta) = c("group","FamilyID","TechReplRelID")
  meta$SubjectID = paste(meta$group,meta$FamilyID,sep=".")
  meta$SampleID = row.names(meta)
  attr.names = names(meta)
  meta.data = join_count_df(list(count=data,attr=meta))
  return(list(data=meta.data,attr.names=attr.names,feat.attr=feat.attr))
}

power.pieper.t1d <- function(
  n = 50,
  alpha.sim = 0.05,
  alpha.orig = 0.05,
  min.mean = 100,
  mult.adj = "fdrtool",
  data.file="aim3/September 28 analysis_T1D.txt",
  R = 4000,
  feat.sig.names=NULL,
  targeted=F
) {
  if(F) {
    taxa.meta = read.pieper.t1d(file.name=data.file)
    aggr_var = "SubjectID"
    taxa.meta = aggregate.by.meta.data(taxa.meta$data,
                                       aggr_var,
                                       taxa.meta$attr.names)
    
    make.global(taxa.meta)
  }
  taxa.meta.attr.names = taxa.meta$attr.names    
  dim.data.orig = dim(count_matr_from_df(taxa.meta$data,taxa.meta.attr.names))
  taxa.meta.data.raw = count.filter(taxa.meta$data,col_ignore=taxa.meta.attr.names,
                                    min_max_frac=0,min_max=0,min_mean=min.mean,min_row_sum=0,other_cnt=NULL)
  print(dim(taxa.meta.data.raw))
  
  taxa.meta.data = norm.meta.data(taxa.meta.data.raw,col_ignore=taxa.meta.attr.names,norm.func=ihs)
  make.global(taxa.meta.data)
  
  #clade.names = get.clade.names(taxa.meta.data,taxa.meta.attr.names)
  #print(clade.names)
  #pvals = wilcox.test.multi(data=taxa.meta.data,resp.vars=clade.names,group.var="group",subset=NULL)
  m = count_matr_from_df(taxa.meta.data,taxa.meta.attr.names)
  
  if(!is.null(feat.sig.names)) {
    ind.sig = which(colnames(m) %in% feat.sig.names)
  }
  else {
    #ind.sig = which(pvals.adj<=alpha.orig)
    ind.sig = which(eff.raw <= 0.66 | eff.raw >= 1.5)
  }
  
  if(targeted) {
    m = m[,ind.sig]
    ind.sig=1:length(ind.sig)
  }
  
  print(length(ind.sig))
  #return(0)  
  
  dim.data.filt = dim(m)
  #mtp.res = MTP(X=t(m),Y=taxa.meta.data$group,get.adjp=F)
  tr.m = t(m)
  group = taxa.meta.data$group
  make.global(m)
  make.global(group)
  m.raw = count_matr_from_df(taxa.meta.data.raw,taxa.meta.attr.names)
  eff.raw = group.mean.ratio(m.raw,group)
  
  pvals = wilcox.test.multi.fast(tr.m,group)  
  #make.global(pvals)
  
  if(mult.adj=="fdrtool") {
    pvals.adj = fdrtool(pvals,statistic="pvalue",plot=F,verbose=F)$qval
  }
  else {
    pvals.adj = p.adjust(pvals,method=mult.adj)
  }
  
  #make.global(pvals.adj)
  
  group.mean = count_matr_from_df(count.summary(m,mean,group),c("group"))
  group.sd = count_matr_from_df(count.summary(m,sd,group),c("group"))
  group.mean.sig = group.mean[,ind.sig]
  #print(group.mean.sig)
  group.sd.sig = group.sd[,ind.sig]
  #print(group.sd.sig)
  #group.n = c(table(group))[rownames(group.mean.sig)]
  group.n = count(group)
  cohens.d = foreach(i.var=seq(ncol(group.mean.sig)),.combine=c) %do% {
    cohens.d.from.mom(mean.gr=group.mean.sig[,i.var],var.gr=group.sd.sig[,i.var]**2,n.gr=group.n)
  }
  #print(cohens.d)
  pvals.boot = boot(data=m,
                    statistic=booted.wilcox.test.multi.fast,
                    R=R,
                    n=n,
                    strata=group,
                    group=group,
                    test.method=wilcox.test.multi.fast)$t
  ## some numerical instability creates pvals that are slightly larger than 1
  pvals.boot[pvals.boot>1] = 1  
  make.global(pvals.boot)
  if(mult.adj=="fdrtool") {
    pvals.boot.adj = aaply(pvals.boot,1,function(x) fdrtool(x,statistic="pvalue",plot=F,verbose=F)$qval)
  }
  else {
    pvals.boot.adj = aaply(pvals.boot,1,p.adjust,method=mult.adj)
  }
  make.global(pvals.boot.adj)
  power.sig = colMeans(pvals.boot.adj[,ind.sig] <= alpha.sim)
  names(power.sig) = names(pvals.adj[ind.sig])
  make.global(power.sig)
  #print(power.sig)
  #print(mean(power.sig))
  return(list(cohens.d=cohens.d,
              mean.cohens.d=mean(cohens.d),
              group.mean.sig=group.mean.sig,
              group.sd.sig=group.sd.sig,
              group.n=group.n,
              power.sig=power.sig,
              mean.power.sig=mean(power.sig),
              pvals.adj.sig=pvals.adj[ind.sig],
              dim.data.orig=dim.data.orig,
              dim.data.filt=dim.data.filt,
              alpha.orig=alpha.orig,
              alpha.sim=alpha.sim))
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
  power.res = power.pieper.t1d(
    n = 50,
    data.file="aim3/September 28 analysis_T1D.txt"
  )
  print("Power results Aim 3:")
  print(power.res)
  mom = pediatric.cancer.2013.aim2()
  print("Power results Aim 2:")
  print(paste("Cohen's d Lactoferrin:",mom$cohens.d))
  print(paste("Mean Cohen's d for significant proteins from Aim 3:",power.res$mean.cohens.d))
}

power.pieper.prostate.cancer.2014<-function() {
  power.res = power.pieper.t1d(
    n = 170,
    ## this is under the proposal directory
    data.file = "data/T1D_proteome/Original Collapesed APEX (All information).AT.tsv",
    min.mean = 1200,
    alpha.sim = 0.05,
    alpha.orig = 0.05,
    R = 400
  )
  print("Power results:")
  print(power.res)
}

power.madupu.kidney_diabetes<-function() {
  prot.ids = as.data.frame(
    matrix(
      c(
        "MAN2B1","O00754",
        "GP5","P40197",
        "FUCA2","Q9BTY2",
        "ATP1A1","P05023",
        "CDH5","P33151",
        "ACE2","Q9BYF1"
      ),
      ncol=2,
      byrow=T
    )
  )
  names(prot.ids) = c("name","uniref100")
  make.global(prot.ids)
  power.res = power.pieper.t1d(
    n = 300,
    ## this is under the proposal directory
    data.file = "data/T1D_proteome/Original Collapesed APEX (All information).AT.tsv",
    min.mean = 0,
    alpha.sim = 0.05,
    alpha.orig = 0.05,
    R = 400,
    mult.adj="bonferroni",
    feat.sig.names=prot.ids$uniref100,
    targeted=T
  )
  print("Power results:")
  print(power.res)  
}
