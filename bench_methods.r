## Supporting methods for comparing multivariate abundance datasets with some notion of ground truth

drop.zero.columns <- function(x) {
  x[,colSums(x)!=0]
}

g.test.pairwise <- function(x,p=NULL,p.adjust.method="BY") {
  if(!is.matrix(x)) {
    stop("x must be a matrix")
  }
  if(nrow(x)<2) {
    stop("x must have at least two rows")
  }
  method = NA
  n.row = nrow(x)
  if(is.null(p)) {
    p.val = matrix(nrow=n.row,ncol=n.row,dimnames=list(rownames(x),rownames(x)))
    cram.v = p.val
    ## we do diagonal as a sanity check - should ge pval 1 and effect size 0
    for(i in seq(n.row)) for(j in seq(i,n.row)) {
      t.res = g.test(drop.zero.columns(x[c(i,j),]))
      p.val[i,j] = p.val[j,i] = t.res$p.value
      cram.v[i,j] = cram.v[j,i] = cramers.v(t.res)
      method = t.res$method
    }
    ## we temporarily set diagonal as NA to exclude it from multiple testing adjustment because
    ## it is deterministically fixed
    diag(p.val) = NA
    p.val.adj = matrix(p.adjust(p.val,method=p.adjust.method),nrow=n.row,ncol=n.row)
    dimnames(p.val.adj) = dimnames(p.val)
    diag(p.val.adj) = 1
    diag(p.val) = 1
  }
  else {
    p.val = rep(NA,n.row)
    names(p.val) = rownames(x)
    cram.v = p.val
    for(i in seq(n.row)) {
      t.res = g.test(x[i,],p=p)
      make.global(t.res)
      p.val[i] = t.res$p.value
      cram.v[i] = cramers.v(t.res)
      p.val.adj = p.adjust(p.val,method=p.adjust.method)
      names(p.val.adj) = names(p.val)
      method = t.res$method
    }
  }
  return (structure(list(p.value=p.val,
                         p.value.adj=p.val.adj,
                         cramers.v=cram.v,
                         method=method,
                         p.adjust.method=p.adjust.method),
                    class="htest"))
}

report.test.pairwise <- function(x,descr) {
  report.section = report$get.section()
  descr = sprintf("Pairwise test of %s (%s)",descr,x$method)
  report$add.header(descr)
  report$push.section(report.section)
  
  call.rep.method = function(y,name=NULL,...) {
    
    if(is.matrix(y)) {
      report$add.table(y,...)
    }
    else {
      report$add.vector(y,name=name,...)
    }
  }
  
  call.rep.method(x$p.value,
                  name="p.value",
                  caption=sprintf("%s. Unadjusted p-value",descr),
                  show.row.names=T)
  call.rep.method(x$p.value.adj,
                  name="p.value.adj",
                  caption=sprintf("%s. Adjusted p-value (method %s).",descr,x$p.adjust.method),
                  show.row.names=T)
  call.rep.method(x$cramers.v,
                  name="cramers.v",
                  caption=sprintf("%s. Cramer's V",descr),show.row.names=T)
  report$pop.section()
}

test.multinom.counts <- function(m_a,
                                 ground.truth.sample.id,
                                 p.adjust.method="BY") {
  gt.ind = which(rownames(m_a$count) == ground.truth.sample.id)
  all.cnt = m_a$count
  cnt = all.cnt[-gt.ind,]
  ## test is undefined when column has no information
  cnt = drop.zero.columns(cnt)
  ## We use all non-zero columns to pairwise compare samples excluding the ground truth samples.
  ## Thus, this comparison includes Unexpected.Taxa column
  pair.indep = g.test.pairwise(cnt,p.adjust.method=p.adjust.method)
  ## Goodness of fit test va ground truth proportions does not make sense for zero probability cells, 
  ## so we have to drop Unexpected.Taxa column. Thus, we can get good fit even if a lot of sequences
  ## were put into unexpected taxa. Assigning some arbitrary small value to p[Unexpected.Taxa] would not help
  ## because it would mean that we know a prior on this type of misclassification. If it is too small,
  ## even a few misclassified sequences can cause significant test results.
  all.cnt = all.cnt[,all.cnt[gt.ind,]>0]
  gt.p = norm.prop(all.cnt[gt.ind,])
  cnt = all.cnt[-gt.ind,]
  good.of.fit = g.test.pairwise(cnt,p=gt.p,p.adjust.method=p.adjust.method)
  res = list(pair.indep=pair.indep,good.of.fit=good.of.fit)
  report.test.pairwise(res$good.of.fit,"Goodness of fit, zero probability categories are dropped")
  report.test.pairwise(res$pair.indep,"Independence")
  return(res)
}

multinom.ci.matrix <- function(x,conf.level=0.95,...) {
  library(DescTools)
  aaply(x,1,function(y) MultinomCI(y,conf.level=conf.level,...))
}


plot.dissim.to.gt <- function(dist.m,caption) {
  pl.data = dist.m[-1,1,drop=F]
  pl.data = data.frame(SampleID=factor(rownames(pl.data)),Dissimilarity=pl.data[,1])
  min.x = min(pl.data$Dissimilarity)
  pl.data$min.x = min.x
  max.x = max(pl.data$Dissimilarity)
  pl.data$max.x = max.x
  pl = ggplot(pl.data, aes(x = Dissimilarity, y = SampleID)) +  
    geom_point(size=rel(4)) +
    geom_segment(aes(xend = Dissimilarity,yend=SampleID,y=SampleID,
                     x=min.x-(max.x-min.x)*0.1)) +
    #geom_bar(stat="identity",width=0.8) +
    scale_x_continuous(limits=c(
      min.x-(max.x-min.x)*0.1,
      max.x)
    ) +
    #coord_flip() +
    theme(axis.title=element_blank(),
          axis.text.y=element_text(color=c("black","black")),
          plot.title = element_text(size = rel(2)),
          axis.text.x = element_text(size = rel(2)), #,angle=90
          axis.text.y = element_text(size = rel(1)))
  
  report$add(pl,caption=paste(caption,"relative to ground truth"))      
  return(pl)
}

discrete.distro.dissim.report <- function(m_a,m_a.norm,reference.rowid,distro.type="multinomial",
                                          vegdist.method = "manhattan",
                                          decostand.method = "range",                                       
                                          do.plot.profiles=F,plot.profiles.task=NULL,sub.report=T) {
  
  report.section = report$add.header("Distance measures between estimated distributions",section.action="push",sub=sub.report)
  
  if(distro.type=="multinomial") {
    ## This adds a small offset to the matrix for zero elements, so that
    ## the distance with self is defined (it has to compute x*log(x+y)).
    dist.m = as.matrix(dist.js(m_a.norm$count))
    caption="Jensen-Shannon dissimilarity"
    report$add.table(dist.m,show.row.names=T,
                     caption=caption)
    plot.dissim.to.gt(dist.m=dist.m,caption=caption)
    ## This dissimilarity is easier to understand in the continuous case,
    ## where it is 1 minus an integral of a product of densities of two variables.
    ## In the discreet case, the rational is that sqrt(p) is a unit vector in L2 norm,
    ## and the distance is defined as Euclidian distance between sqrt(p).
    dist.m = as.matrix(vegdist(sqrt(m_a.norm$count),method = "euclidean")/sqrt(2))
    caption="Hellinger dissimilarity"
    report$add.table(dist.m,show.row.names=T,
                     caption=caption)
    plot.dissim.to.gt(dist.m=dist.m,caption=caption)
    res.mult = test.multinom.counts(m_a,reference.rowid)
  }
  
  m_a.norm$count = decostand(m_a.norm$count,
                             method=decostand.method,
                             MARGIN=2)
  
  if(decostand.method == "standardize") {
    ## if all values are equal, we get NaNs; replace with 0.
    m_a.norm$count[is.nan(m_a.norm$count)] = 0.
  }
  
  decostand.descr = paste("Taxa standardized to equal weight using method",decostand.method)
  
  if(do.plot.profiles) {
    plot.profiles.task = within(plot.profiles.task, {
      show.profile.task=within(show.profile.task, {
        geoms=c("bar")
      })
    })
    
    do.call(plot.profiles,
            c(list(m_a=m_a.norm,
                   ci=NULL,
                   feature.order=list(list(ord=colnames(m_a.norm$count),
                                           ord_descr="Mean abundance")),
                   feature.descr=paste(decostand.descr,aggr.descr,sep=". ")),
              plot.profiles.task
            )
    )
  }
  dist.m = as.matrix(vegdist(m_a.norm$count,method=vegdist.method))
  caption=sprintf("%s distance. %s.",vegdist.method,decostand.descr)
  report$add.table(dist.m,show.row.names=T,
                   caption=caption)
  plot.dissim.to.gt(dist.m=dist.m,caption=caption)
  report$pop.section()  
}
