library(hdi)
library(reshape2)

cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
registerDoSNOW(cl)  

stab_sel = function(do.perm=F) {
  y = m_a.norm$attr$T1D
  if(do.perm) {
    y = sample(y)
  }
  stability(x=m_a.norm$count,
            y=y,
            EV=1,
            args.model.selector=list(standardize=T,family="binomial")
  ) 
}

freq = stab_sel(do.perm=F)$freq
freq.perm = foreach(i=1:300,
        .combine=cbind,
        .packages=c("hdi")) %dopar% {
  stab_sel(do.perm=T)$freq
}

freq.rank = freq[order(freq,decreasing = T)]
freq.perm.rank = t(aaply(freq.perm,2,function(x) x[order(x,decreasing = T)]))
rownames(freq.perm.rank) = NULL

freq.perm.rank.mlt = melt(freq.perm.rank)
freq.perm.rank.mlt$series = "perm"

freq.rank.mlt = data.frame(
  Var1=seq(length(freq.rank)),
  X1 = "samp",
  value=freq.rank,
  series="samp")

freq.mlt = rbind(freq.rank.mlt,freq.perm.rank.mlt)

ggplot(freq.mlt, 
       aes(Var1,value)) + 
  scale_x_continuous(limits=c(0, 10)) + 
  geom_line(aes(color = series))

rank.p.val = rowMeans(freq.rank >= freq.perm.rank)
names(rank.p.val) = names(freq.rank)

## how many features in permutations, on average, have as good prob of 
##selection as the best feature in the sample
print(mean(colSums(freq.perm.rank >= 0.79)))

call_hdi<-function() {
  ##it only understands boolean vars for bionomial models.
  ##can also use ridge.proj here
  fit.proj = lasso.proj(x=m_a.norm$count,y=m_a.norm$attr$T1D=="T1D",standardize=T,family="binomial",multiplecorr.method="WY")
  head(fit.proj$pval.corr[order(fit.proj$pval.corr)])
  head(confint(fit.proj,level=0.95)[order(fit.proj$pval.corr),])
  fit.proj = multi.split(x=m_a.norm$count,y=m_a.norm$attr$T1D=="T1D",
                         classical.fit=glm.pval,args.model.
                         selector=list(standardize=T,family="binomial"),
                         args.classical.fit=list(family="binomial"))
}

call_stabs<-function() {
  require(stabs)
  fit.proj = stabsel(x=m_a.norm$count,y=m_a.norm$attr$T1D=="T1D",fitfun=glmnet.lasso,
                     args.fitfun=list(alpha=0.05,standardize=T,family="binomial"),
                     cutoff=0.9,PFER=1,sampling.type="SS",assumption="r-concave")
}