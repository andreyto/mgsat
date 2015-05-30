##For different runs of pipelines , we could use
##McNemar's test because the reads are the same - we are putting them
##into bins in each run. Then, we need to get the initial number of read pairs,
##and create a (giant) bin for all rejected reads. But then, this would
##be a useless test, because the difference will be always significant due to
##a different ratio of reads rejected at each run.

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

load.ground.truth.m_a <- function(abund.file,
                                   aggr.type=c("genus.abund","genus.otus","otu"),
                                   pseudo.sum=50000) {
  x = read.delim(abund.file,header=T,sep="\t")
  x = x[,c("Organism","Genus","Profile","WGS.Proportion.16S.Gene","ProfileID")]
  names(x)[4] = "Prop"
  x$Prop = round(x$Prop * pseudo.sum)
  x$Organism = sanitize.taxa.names(x$Organism)
  x$Genus = sanitize.taxa.names(x$Genus)
  if(aggr.type %in% c("genus.abund","genus.otus")) {
    aggr = ddply(x,c("Profile","ProfileID","Genus"),summarise,Abund=sum(Prop),Count=sum(Prop>0))
    if(aggr.type=="genus.abund") {
      value.var = "Abund"
    }
    else if(aggr.type=="genus.otus") {
      value.var = "Count"
    }
    count = acast(aggr,ProfileID~Genus,value.var=value.var)
  }
  else if(aggr.type == "otu") {
    count = acast(x,ProfileID~Organism,value.var="Prop")
  }
  attr = unique(x[,c("ProfileID","Profile")])
  rownames(attr) = attr$ProfileID
  attr = attr[rownames(count),]
  m_a = list(count=count,attr=attr)
  m_a
}


load.wl.cdhit.m_a <- function(data.file,meta.file=NULL,aggr.type=c("genus.abund","genus.otus","otu")) {
  x = read.delim(data.file,header=T,sep="\t",stringsAsFactors=F)
  ## extract all fields up to total field plus genus
  x = x[,c(names(x)[seq(match("total",names(x))-1)],"family","genus")]
  
  ## make taxa names compatible with Mothur output
  mask.rec = (x$genus == "Incertae Sedis")
  x$genus[mask.rec] = paste(x$family[mask.rec],"incertae sedis",sep="_")
  mask.rec = grepl("uncultured.*",x$genus,ignore.case=T) | (!grepl("\\S+",x$genus))
  x$genus[mask.rec] = paste("Unclassified",x$family[mask.rec],sep="_")
  
  x$genus = gsub("\\s","_",x$genus)
  x$genus = gsub("-","_",x$genus) # Escherichia-Shigella
  x$genus[x$genus == "Clostridium"] = "Clostridium_sensu_stricto"
  
  x$OTU = paste(x$OTU,x$genus,sep=".")
  rownames(x) = x$OTU
  x$OTU = NULL
  x$family = NULL
  
  if(aggr.type=="otu") {
    x$genus = NULL
    x = t(x)
  }
  else {
    aggr.func = switch(aggr.type,
                       genus.abund = sum,
                       genus.otus = function(y) sum(y>0)
    )
    x = aggregate(x[,!("genus" == names(x))],list(genus = x$genus),aggr.func)
    
    rownames(x) = x$genus
    x$genus = NULL
    x = t(x)
  }
  
  meta = load.meta.default(meta.file)
  return (merge.counts.with.meta(x,meta))
}

benchmark.load.data <- function(otu.count.filter.options,
                                taxa.level,
                                count.basis,
                                col.match,
                                aggr.type,
                                aggr.descr
) {
  
  report.section = report$add.header(sprintf("Loading annotation files for taxonomic aggregation level %s",
                                             aggr.type),
                                     section.action="push", sub=T)
  
  meta.file.samples = "refdata/bench_bei_meta.txt"
  
  res = with(mgsat.16s.task.template,{
    
    runs.files = list(
#             list(
#               RunID="YAP.comm.V13",
#               run.descr="YAP currently deployed in common area; V13",
#               read.data.task=within(read.data.task, {
#                 taxa.summary.file = "yap/bei/v13/common/7913675baaa4dd38c5ad2c19b10bdda5.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
#                 otu.shared.file="yap/bei/v13/common/5f44040ee0ae9f533d83e48deded4bac.files_x1.sorted.0.03.shared"
#                 cons.taxonomy.file="yap/bei/v13/common/dce008a0990dc67e8148b3a31fb0bb02.files_x1.sorted.0.03.cons.taxonomy"
#                 taxa.summary.file.otu = "yap/bei/v13/common/0211bc387719c7e0a4dc9d4ea0eccee9.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
#                 meta.file=meta.file.samples
#               })
#             )
      #       ,  
      #   list(
      #     RunID="YAP.base",
      #     read.data.task=within(read.data.task, {
      #       taxa.summary.file = "yap/bei/v13/baseline/bd568ff2897e717e1ab6c948f9509344.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      #       otu.shared.file="yap/bei/v13/baseline/21c401833a6ba4d011974fa46f05e476.files_x1.sorted.0.03.shared"
      #       cons.taxonomy.file="yap/bei/v13/baseline/ff4e455254bb6759e890e170151a4e37.files_x1.sorted.0.03.cons.taxonomy"
      #       taxa.summary.file.otu = "yap/bei/v13/baseline/1baf5634efbb3aa711408089eda1faf5.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"      
      #       meta.file=meta.file.samples
      #     })
      #   ),
      #   list(
      #     RunID="YAP.mc.ch1",
      #     read.data.task=within(read.data.task, {
      #       taxa.summary.file = "yap/bei/v13/make.contigs.chim_full_ref/5526056cfd0ce087ba07cb3377c90b4f.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      #       otu.shared.file="yap/bei/v13/make.contigs.chim_full_ref/60d6c77f675da1ab6bada2c013e8a8aa.files_x1.sorted.0.03.shared"
      #       cons.taxonomy.file="yap/bei/v13/make.contigs.chim_full_ref/3ca693e09ea5fb1d1e769bb5e7cef101.files_x1.sorted.0.03.cons.taxonomy"
      #       taxa.summary.file.otu = "yap/bei/v13/make.contigs.chim_full_ref/83324aef0f7627d8195b55efa56e5fb7.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
      #       meta.file=meta.file.samples
      #     })
      #   ),  
      #       list(
      #         RunID="YAP.upd",
      #         run.descr="YAP updated to use Mothur make.contigs and Mothur reference DB v.10",    
      #         read.data.task=within(read.data.task, {
      #           taxa.summary.file = "yap/bei/v13/make.contigs.chim_gold/d4db63db7e566943c9ba8277f0d6b253.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      #           otu.shared.file="yap/bei/v13/make.contigs.chim_gold/57782e4f219e8f2212aa5e3ba3f91452.files_x1.sorted.0.03.shared"
      #           cons.taxonomy.file="yap/bei/v13/make.contigs.chim_gold/4021347ead1b4c49a7723f682c9362a1.files_x1.sorted.0.03.cons.taxonomy"
      #           taxa.summary.file.otu = "yap/bei/v13/make.contigs.chim_gold/e7cdd0d2622b10ea7d2a58666eb713f4.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
      #           meta.file=meta.file.samples
      #         })
      #       ),    
      #       list(
      #         RunID="YAP.upd.prefilt",
      #         run.descr="YAP updated to use Mothur make.contigs, Mothur reference DB v.10 and drop singletons after pre-clustering step",    
      #         read.data.task=within(read.data.task, {
      #           taxa.summary.file = "yap/bei/v13/min_precluster_size//3b8c4c926594400416e38756cd299875.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      #           otu.shared.file="yap/bei/v13/min_precluster_size//f6e65ff8146987bcfacb533a2e2a95c9.files_x1.sorted.0.03.shared"
      #           cons.taxonomy.file="yap/bei/v13/min_precluster_size//3b8c4c926594400416e38756cd299875.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
      #           taxa.summary.file.otu = "yap/bei/v13/min_precluster_size//4979c38da8610fd77c1ca50da226a859.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
      #           meta.file=meta.file.samples
      #         })
      #       ),    
      #       list(
      #         RunID="YAP.upd.split",
      #         run.descr="YAP updated to use Mothur make.contigs, Mothur reference DB v.10, three-stage split pre-clustering step and drop singletons after pre-clustering step",    
      #         read.data.task=within(read.data.task, {
      #           taxa.summary.file = "yap/bei/v13/min_precluster_size_split/9b08b971d51b7ed71c3dc1e968778957.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      #           otu.shared.file="yap/bei/v13/min_precluster_size_split/51fc031dfaeff9acc580d26676cd2ce4.files_x1.sorted.0.03.shared"
      #           cons.taxonomy.file="yap/bei/v13/min_precluster_size_split/9b08b971d51b7ed71c3dc1e968778957.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
      #           taxa.summary.file.otu = "yap/bei/v13/min_precluster_size_split/1ba957695051f8733a71572a5cadc143.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
      #           meta.file=meta.file.samples
      #         })
      #       )
      #       ,          
#       list(
#         RunID="YAP.rc.2015-04-09",
#         run.descr="YAP RC 2015-04-09",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = "yap/bei/v13/rc.2015-04-09/d95c583431096967cf6279fdb7feb3cd.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
#           otu.shared.file="yap/bei/v13/rc.2015-04-09/d41c472173ef00e1ce21205c22111779.files_x1.sorted.0.03.shared"
#           cons.taxonomy.file="yap/bei/v13/rc.2015-04-09/d95c583431096967cf6279fdb7feb3cd.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
#           taxa.summary.file.otu = "yap/bei/v13/rc.2015-04-09/a67ca553bd7a28470116ff39c303c03d.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
#           meta.file=meta.file.samples
#         })
#       )
#       ,                
#       list(
#         RunID="YAP.clean_01",
#         run.descr="YAP cleaning 01",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v13/clean_01/ee3bf4399aa1b11b4a666b7fd3779b4b.files_x1.sorted.0.03.shared"
#           cons.taxonomy.file="yap/bei/v13/clean_01/b414a7c2dba8a90d1d6c6ab49823d0ac.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )
      #      ,                
      #      list(
      #        RunID="YAP.clean_02",
      #        run.descr="YAP cleaning 02",    
      #        read.data.task=within(read.data.task, {
      #          taxa.summary.file = NA
      #          otu.shared.file="yap/bei/v13/clean_02/b57962317b029feaea759068656731f0.files_x1.sorted.0.03.shared"
      #            cons.taxonomy.file="yap/bei/v13/clean_02/c815c2f7b84bb80f5d11ef854a9f8ab7.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
      #          taxa.summary.file.otu = NA
      #          meta.file=meta.file.samples
      #        })
      #      )     
#       ,
#       list(
#         RunID="YAP.clean_03",
#         run.descr="YAP cleaning 03",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v13/clean_03/68e401b15d999da15ae0deee7abf20cc.files_x1.sorted.0.03.shared"
#           cons.taxonomy.file="yap/bei/v13/clean_03/1002f15a8e73fa19f39c5e69e296c23c.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )     
#       ,
#       list(
#         RunID="YAP.clean_05",
#         run.descr="YAP cleaning 05",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v13/clean_05/7ab5391481764492bbe645c317078d37.files_x1.sorted.0.03.shared"
#           cons.taxonomy.file="yap/bei/v13/clean_05/862aaa213d961e10074a8f39a8fbea12.files_x1.sorted.0.03.cons.taxonomy.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )
#      ,
      list(
        RunID="YAP.rc.2015-04-18.V13",
        run.descr="YAP RC 2015-04-18 V13",    
        read.data.task=within(read.data.task, {
          taxa.summary.file = NA
          otu.shared.file="yap/bei/v13/rc.2015-04-18/*.shared"
          cons.taxonomy.file="yap/bei/v13/rc.2015-04-18/*.seq.taxonomy"
          taxa.summary.file.otu = NA
          meta.file=meta.file.samples
        })
      )
#       ,
#       list(
#         RunID="YAP.rc.2015-04-18.V4.Bill",
#         run.descr="Bill's dataset; YAP RC 2015-04-18 V4",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v4/rc.2015-04-18/*.shared"
#           cons.taxonomy.file="yap/bei/v4/rc.2015-04-18/*.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )           
,
list(
  RunID="YAP.rc.2015-05-30.V4.Bill",
  run.descr="Bill's dataset; YAP RC 2015-05-30 V4",    
  read.data.task=within(read.data.task, {
    taxa.summary.file = NA
    otu.shared.file="yap/bei/v4/rc.2015-05-30/*.shared"
    cons.taxonomy.file="yap/bei/v4/rc.2015-05-30/*.seq.taxonomy"
    taxa.summary.file.otu = NA
    meta.file=meta.file.samples
  })
)           
#       ,
#       list(
#         RunID="YAP.rc.2015-04-18.no_prefilt.V4.Bill",
#         run.descr="Bill's dataset; YAP RC 2015-04-18 singletons ARE NOT dropped after pre-clustering; V4",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v4/rc.2015-04-18.no_prefilt/*.shared"
#           cons.taxonomy.file="yap/bei/v4/rc.2015-04-18.no_prefilt/*.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )           
#       ,
#       list(
#         RunID="YAP.rc.2015-04-18.V4.Bill.no_pcr_ref",
#         run.descr="Bill's dataset; YAP RC 2015-04-18 V4 No reference PCR",    
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = NA
#           otu.shared.file="yap/bei/v4/rc.2015-04-18.no_pcr_ref/*.shared"
#           cons.taxonomy.file="yap/bei/v4/rc.2015-04-18.no_pcr_ref/*.seq.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )           
      ,
      list(
        RunID="YAP.rc.2015-05-30.V4.default.Alex",
        run.descr="Alex's dataset one pair of samples; YAP RC 2015-05-30 V4 Default",    
        read.data.task=within(read.data.task, {
          taxa.summary.file = NA
          otu.shared.file="alex/v4/yap/4_4/default/*.shared"
          cons.taxonomy.file="alex/v4/yap/4_4/default/*.seq.taxonomy"
          taxa.summary.file.otu = NA
          meta.file=meta.file.samples
        })
      )
,
list(
  RunID="YAP.rc.2015-05-30.V4.default.pooled.Alex",
  run.descr="Alex's dataset all samples pooled by community type; YAP RC 2015-05-30 V4 Default",    
  read.data.task=within(read.data.task, {
    taxa.summary.file = NA
    otu.shared.file="alex/v4/yap/pooled/default/*.shared"
    cons.taxonomy.file="alex/v4/yap/pooled/default/*.seq.taxonomy"
    taxa.summary.file.otu = NA
    meta.file=meta.file.samples
  })
)           
# ,
# list(
#   RunID="YAP.rc.2015-05-30.V4.Alex.no_prefilt",
#   run.descr="Alex's dataset; YAP RC 2015-05-30 V4 No singleton removal",    
#   read.data.task=within(read.data.task, {
#     taxa.summary.file = NA
#     otu.shared.file="alex/v4/yap/4_4/no_sing_remove/*.shared"
#     cons.taxonomy.file="alex/v4/yap/4_4/no_sing_remove/*.seq.taxonomy"
#     taxa.summary.file.otu = NA
#     meta.file=meta.file.samples
#   })
# )           
,     
      list(
        RunID="Mothur.01.V13",
        run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP by Andrey (reference-based chimera removal) V13",
        read.data.task=within(read.data.task, {
          taxa.summary.file = "sop/v13/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
          otu.shared.file="sop/v13/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
          cons.taxonomy.file="sop/v13/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
          taxa.summary.file.otu = "bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
          meta.file=meta.file.samples
        })
      )
#       #       ,
#       #      list(
#       #        RunID="Mothur.prefilt.V13",
#       #        run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP but singletons dropped after pre.cluster()",
#       #        read.data.task=within(read.data.task, {
#       #          taxa.summary.file = "bei/v13.prefilt/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.tax.summary"
#       #          otu.shared.file="bei/v13.prefilt/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared"
#       #          cons.taxonomy.file="bei/v13.prefilt/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy"
#       #          taxa.summary.file.otu = NA
#       #          meta.file=meta.file.samples
#       #        })
#       #      )
#       ,
#       list(
#         RunID="Mothur.auto.V13",
#         run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons dropped after pre.cluster(), automated most parameter selection V13",
#         read.data.task=within(read.data.task, {
#           taxa.summary.file = "sop/v13/v13.oligos/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.tax.summary"
#           otu.shared.file="sop/v13/v13.oligos/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.shared"
#           cons.taxonomy.file="sop/v13/v13.oligos/input.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.an.unique_list.0.03.cons.taxonomy"
#           taxa.summary.file.otu = NA
#           meta.file=meta.file.samples
#         })
#       )
# #       ,
# #       list(
# #         RunID="Mothur.auto.V4.Bill",
# #         run.descr="Bill's dataset; Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons dropped after pre.cluster(), automated most parameter selection V4",
# #         read.data.task=within(read.data.task, {
# #           taxa.summary.file = NA
# #           otu.shared.file="sop/v4/auto/*.shared"
# #           cons.taxonomy.file="sop/v4/auto/*.0.03.cons.taxonomy"
# #           taxa.summary.file.otu = NA
# #           meta.file=meta.file.samples
# #         })
# #       )
,
list(
  RunID="Mothur.no_prefilt.V4.Bill",
  run.descr="Bill's dataset; Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons ARE NOT dropped after pre.cluster(), hand-tuned parameter selection; V4",
  read.data.task=within(read.data.task, {
    taxa.summary.file = NA
    otu.shared.file="sop/v4/manual.no_prefilt/*.shared"
    cons.taxonomy.file="sop/v4/manual.no_prefilt/*.0.03.cons.taxonomy"
    taxa.summary.file.otu = NA
    meta.file=meta.file.samples
  })
)
,
list(
  RunID="Mothur.no_prefilt.V4.Alex",
  run.descr="Alex's dataset one pair of samples; Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons ARE NOT dropped after pre.cluster(), hand-tuned parameter selection; V4",
  read.data.task=within(read.data.task, {
    taxa.summary.file = NA
    otu.shared.file="alex/v4/sop/4_4/default/*.shared"
    cons.taxonomy.file="alex/v4/sop/4_4/default/*.0.03.cons.taxonomy"
    taxa.summary.file.otu = NA
    meta.file=meta.file.samples
  })
)
,
list(
  RunID="Mothur.no_prefilt.V4.pooled.Alex",
  run.descr="Alex's dataset samples pooled by community profile; Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons ARE NOT dropped after pre.cluster(), hand-tuned parameter selection; V4",
  read.data.task=within(read.data.task, {
    taxa.summary.file = NA
    otu.shared.file="alex/v4/sop/pooled/default/*.shared"
    cons.taxonomy.file="alex/v4/sop/pooled/default/*.0.03.cons.taxonomy"
    taxa.summary.file.otu = NA
    meta.file=meta.file.samples
  })
)
# ,
# list(
#   RunID="Mothur.prefilt.V4.Bill",
#   run.descr="Bill's dataset; Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP, singletons ARE dropped after pre.cluster(), hand-tuned parameter selection; V4",
#   read.data.task=within(read.data.task, {
#     taxa.summary.file = NA
#     otu.shared.file="sop/v4/manual/*.shared"
#     cons.taxonomy.file="sop/v4/manual/*.0.03.cons.taxonomy"
#     taxa.summary.file.otu = NA
#     meta.file=meta.file.samples
#   })
# )           

      #     ,     
      #   list(
      #   RunID="Mothur.Sarah.01",
      #   read.data.task=within(read.data.task, {
      #     taxa.summary.file = "sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
      #     otu.shared.file="sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
      #     cons.taxonomy.file="sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
      #     taxa.summary.file.otu = "sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
      #     meta.file=meta.file.samples
      #   })
      #  )
      #       ,
      #       list(
      #         RunID="Mothur.Sarah.02",
      #         run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP by Sarah (dataset-based chimera removal)",
      #         read.data.task=within(read.data.task, {
      #           taxa.summary.file = "sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
      #           otu.shared.file="sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
      #           cons.taxonomy.file="sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
      #           taxa.summary.file.otu = "sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
      #           meta.file=meta.file.samples
      #         })
      #       )
      
    )
    
    ground.truth.WGS.file = "refdata/sarah.2015-03-08/derived/bei_abund_v.5.txt"
    ground.truth.run.id = "Ground.Truth.WGS"
    ground.truth.run.descr = "Mock samples sequenced with WGS, reads mapped to reference to compute coverage
which is multiplied by 16S copy number and normalized to proportions to get expected ground truth organism
abundance in the biological sample. Proportions multiplied by a large constant to imitate 16S read counts
so this dataset can be used in diversity and abundance estimates together with actual 16S annotation runs."
    
    runs.data = lapply(runs.files,
                       function(run.files) {
                         with(run.files,{
                           read.data.task$count.basis = count.basis
                           read.data.task$otu.count.filter.options = otu.count.filter.options
                           report$add.header(paste("Loading inputs for run:",run.descr))
                           m_a = do.call(read.data.method,
                                         c(
                                           list(taxa.level=taxa.level),
                                           read.data.task
                                         )
                           )
                           m_a$attr$RunID = RunID
                           m_a$attr$run.descr = run.descr
                           m_a
                         })
                       }
    )
    
    m_a = cbind.m_a(runs.data,batch.attr="RunID",col.match=col.match)
    
    m_a.gt = load.ground.truth.m_a(ground.truth.WGS.file,aggr.type=aggr.type)
    
    m_a.gt$attr$SampleID = m_a.gt$attr$ProfileID
    m_a.gt$attr$RunID = ground.truth.run.id
    m_a.gt$attr$run.descr = ground.truth.run.descr
    m_a.gt$attr$IdSfx = m_a.gt$attr$RunID
    
#     m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.00005.txt",
#                                meta.file=meta.file.samples,
#                                aggr.type=aggr.type)
#     m_a.wl$attr$RunID = "WL.00005.V13"
#     m_a.wl$attr$run.descr = "Weizhong's CD-HIT-OTU pipeline that discarded OTUs below relative abundance of 0.00005 V13"
#     m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
#     m_a.wl.1 = m_a.wl
#     m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.0001.txt",
#                                meta.file=meta.file.samples,
#                                aggr.type=aggr.type)
#     m_a.wl$attr$RunID = "WL.0001.V13"
#     m_a.wl$attr$run.descr = "Weizhong's CD-HIT-OTU pipeline that discarded OTUs below relative abundance of 0.0001 V13"
#     m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
#     m_a.wl.2 = m_a.wl
    
    m_a$attr$IdSfx = ""
    #m_a.norm = norm.count.m_a(m_a,method=norm.method.basic)
    #make.global(m_a.norm)
    m_a.abs = cbind.m_a(list(m_a.gt,m_a),batch.attr="IdSfx",col.match=col.match)
    m_a.abs$attr$IdSfx = NULL
    
    list(m_a.abs=m_a.abs,ground.truth.run.id=ground.truth.run.id)
  })
  
  report$pop.section()
  return (res)
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

benchmark.abund <- function(m_a.abs,
                            aggr.type,
                            aggr.descr,
                            drop.taxa=c(),
                            true.taxa.only,
                            other_cnt,
                            norm.count.task,
                            constrasts.stage.min.feature.prop,
                            do.contrasts,
                            ground.truth.run.id,
                            contrasts.method,
                            plot.profiles.task,
                            do.plot.profiles=T) {
  
  report.section = report$get.section()
  
  #make.global(name="dbgenv")
  
  true.taxa = colMeans(m_a.abs$count[m_a.abs$attr$RunID==ground.truth.run.id,,drop=F])>0
  true.taxa = names(true.taxa)[true.taxa]
  
  m_a.abs = subset.m_a(m_a.abs,select.count=!(colnames(m_a.abs$count) %in% drop.taxa))
  
  if(true.taxa.only) {
    m_a.abs = count.filter.m_a(m_a.abs,drop.zero=T,
                               keep.names=true.taxa,
                               other_cnt=other_cnt)
  }
  
  make.global(m_a.abs)
  m_a.norm = norm.count.report(m_a.abs,
                               res.tests=NULL,
                               descr="Normalization",
                               norm.count.task)
  are.props = (norm.count.task$method == "norm.prop")
  if(do.contrasts) {
    m_a.prop = norm.prop.m_a(m_a.abs)
    m_a.norm = subset.m_a(m_a.norm,select.count=
                            colMeans(m_a.prop$count[
                              m_a.prop$attr$Profile=="EVEN" &
                                m_a.prop$attr$RunID==ground.truth.run.id,,drop=F])>0 &
                            colMeans(m_a.prop$count[
                              m_a.prop$attr$Profile=="STAG" &
                                m_a.prop$attr$RunID==ground.truth.run.id,,drop=F])>=
                            constrasts.stage.min.feature.prop
    )
    contrasts = c(STAG=1,EVEN=-1)
    s.c = sample.contrasts(m_a.norm, group.attr = "Profile", block.attr = "RunID")
    make.global(s.c)
    if(contrasts.method=="diff") {
      m_a.norm = s.c$m_a.contr
    }
    else {
      m_a.norm = contrasts.groups.log.fold.change(s.c$m_a.groups)
    }
    are.props = F
    #levels.last.first = names(s.c$contrasts)
    #m_a.c = s.c$m_a.contr
    #x = m_a.c$count
    #y.relev = rep(1,nrow(x))
  }
  ## CLR transform can create all-zero columns from non-zero columns
  #m_a.norm = subset.m_a(m_a.norm,select.count=colSums(abs(m_a.norm$count))>0)
  
  make.global(m_a.norm)
  gt.rowid = rownames(m_a.norm$attr)[m_a.norm$attr$RunID==ground.truth.run.id]
  
  if(do.plot.profiles) {
    if(!are.props) {
      geoms.profiles=c("bar")
      ci.prop=NULL
      plot.descr = aggr.descr
    }
    else {
      geoms.profiles=c("bar_stacked","bar")
      ci.prop=multinom.ci.matrix(m_a.abs$count)
      ci.prop[gt.rowid,,c("lwr.ci","upr.ci")] = NA
      plot.descr = paste(aggr.descr,
                         "Error bars show simultaneous confidence 
                       intervals for multinomial proportions
                       computed with MultinomCI {DescTools}.
                       They represent uncertainty in estimating
                       the proportions due to random sampling of
                       reads.",sep=". ")
    }
    
    plot.profiles.task = within(plot.profiles.task, {
      show.profile.task=within(show.profile.task, {
        geoms=geoms.profiles
      })
    })
    
    do.call(plot.profiles,
            c(list(m_a=m_a.norm,
                   ci=ci.prop,
                   feature.order=list(list(ord=colnames(m_a.norm$count),
                                           ord_descr="Mean abundance")),
                   feature.descr=plot.descr),
              plot.profiles.task
            )
    )
  }
  report$add.header("Distance measures between samples")
  report$push.section(report.section)
  
  if(are.props) {
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
    dist.m = as.matrix(vegdist(sqrt(m_a.norm$count),method = "euclidian")/sqrt(2))
    caption="Hellinger dissimilarity"
    report$add.table(dist.m,show.row.names=T,
                     caption=caption)
    plot.dissim.to.gt(dist.m=dist.m,caption=caption)
    res.mult = test.multinom.counts(m_a.abs,gt.rowid)
  }
  
  if(norm.count.task$method == "ident" || are.props) {
    vegdist.method = "manhattan"
    decostand.method = "range"
  }
  else {
    vegdist.method = "euclidean"
    decostand.method = "standardize"
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

benchmark.project <- function() {
  
  cl = start.cluster.project()
  
  report$add.descr("16S annotation count matrices from different pipelines and parameter combinations are
                   compared against the expected ground truth.")
  
  drop.taxa = c("Nothing") #c("Clostridium_sensu_stricto") #c("Helicobacter") #c("Streptococcus")
  other_cnt = "Unexpected_Taxa"
  otu.count.filter.options = mgsat.16s.task.template$read.data.task$otu.count.filter.options
  otu.count.filter.options$min_max_frac = 0.0001
  
  do.plot.profiles = F
  
  true.taxa.only = T
  contrasts.methods = c("lfc","diff")
  constrasts.stage.min.feature.prop = 0
  do.contrasts = T
  
  norm.count.task.comp = within(mgsat.16s.task.template$test.counts.task$norm.count.task, {
    method = "norm.clr"
    ##offset 1 only makes sense if we create pseudocounts for ground truth
    method.args = list(offset=1,tol=0.0000001)
    drop.features=list("other")
  })
  
  norm.count.task.prop = within(mgsat.16s.task.template$test.counts.task$norm.count.task, {
    method = "norm.prop"
    method.args = list()
    drop.features=list("other")
  })  
  
  norm.count.task.ident = within(mgsat.16s.task.template$test.counts.task$norm.count.task, {
    method = "ident"
    method.args = list()
    drop.features=list("other")
  })    
  
  report.section = report$get.section()
  
  for(task.aggr in list(
    list(aggr.type="genus.abund",
         aggr.descr="Genus relative abundance",
         norm.count.task=norm.count.task.prop)
        ,
        list(aggr.type="genus.abund",
             aggr.descr="Genus relative abundance, Centered Log Ratio transform",
             norm.count.task=norm.count.task.comp)    
        ,
        list(aggr.type="genus.otus",
             aggr.descr="OTU counts per genus",
             norm.count.task=norm.count.task.ident)
        ,
        list(aggr.type="genus.otus",
             aggr.descr="OTU relative counts per genus",
             norm.count.task=norm.count.task.prop)
#         ,
#         list(aggr.type="otu",
#              aggr.descr="OTU abundance",
#              norm.count.task=norm.count.task.prop)
  )) {
    
    
    with(task.aggr,{
      
      report$add.header(paste("Annotation level:",aggr.descr),
                        report.section=report.section,sub=T)
      
      
      count.basis = "seq"
      taxa.level = 6
      col.match = T
      do.divrich = F
      do.abund = T
      do.summary.meta = F
      
      if(aggr.type == "genus.otus") {
        count.basis = "otu"
        taxa.level = 6
        do.divrich = F
        do.summary.meta = T
        do.contrasts = F
      }
      else if(aggr.type == "otu") {
        taxa.level = "otu"
        col.match = F
        do.summary.meta = T
        do.divrich = T
        do.abund = F
      }
      
      if(norm.count.task$method != "norm.prop") {
        contrasts.methods = c("diff")
      }
      
      
      test.counts.task = within(mgsat.16s.task.template$test.counts.task, {
        divrich.task = within(list(),{
          n.rar.rep=400
          is.raw.count.data=T
          group.attr = "SampleID"
          counts.glm.task = NULL
          beta.task = NULL
          counts.genesel.task = NULL
          do.plot.profiles = T
          do.incidence=F
          do.abundance=T
          do.rarefy=F
          do.accum=F
        })
        
        plot.profiles.task = within(plot.profiles.task, {
          id.vars.list = list(c("SampleID"))
          feature.meta.x.vars=NULL
          do.profile=T
          do.feature.meta=F
          show.profile.task=within(show.profile.task, {
            geoms=c("bar")
            sqrt.scale=T
            stat_summary.fun.y="identity"
          })
        })
        
      })
      
      load.res = benchmark.load.data(otu.count.filter.options=otu.count.filter.options,
                                     taxa.level=taxa.level,
                                     count.basis=count.basis,
                                     col.match=col.match,
                                     aggr.type=aggr.type,
                                     aggr.descr=aggr.descr)
      m_a.abs = load.res$m_a.abs
      ground.truth.run.id = load.res$ground.truth.run.id
      
      attr.rep = unique(m_a.abs$attr[,c("RunID","run.descr")])
      attr.rep$run.descr = gsub("\n"," ",attr.rep$run.descr)
      report$add.table(attr.rep,caption="Description of annotation runs")
      
      report$add.header("Data analysis")
      report$push.section(report.section)
      
      #if(aggr.type == "otu") {
      #  m_a.abs = subset.m_a(m_a.abs,subset=!(m_a.abs$attr$RunID==ground.truth.run.id))
      #}
      
      m_a.abs = count.filter.m_a(m_a.abs,drop.zero=T)
      
      make.global(m_a.abs)
      
      if(do.summary.meta) {
        
        summary.meta.task = within(mgsat.16s.task.template$summary.meta.task, {
          group.vars = NULL
          show.sample.totals=T
          show.sample.means=F
          sub.report=F
        })
        
        do.call(report.sample.count.summary,c(
          list(m_a.abs),
          summary.meta.task
        )
        )
      }
      
      
      if(do.divrich) {
        for(ProfileID in levels(factor(m_a.abs$attr$ProfileID))) {
          m_a.abs.prof = subset.m_a(m_a.abs,subset=(m_a.abs$attr$ProfileID==ProfileID))
          divrich.prof.descr = sprintf("Community profile %s",ProfileID)
          
          for(do.rarefy in c(F,T)) {
            divrich.task = test.counts.task$divrich.task
            divrich.task$do.rarefy = do.rarefy
            divrich.task$extra.header = sprintf("%s %s",divrich.prof.descr,ifelse(do.rarefy,
                                                                                  "with rarefication to common depth",
                                                                                  "without rarefication to common depth"))
            tryCatchAndWarn({ 
              do.call(mgsat.divrich.report,
                      c(list(m_a.abs.prof,
                             plot.profiles.task=test.counts.task$plot.profiles.task),
                        divrich.task)
              )
            })
          }
        }
      }
      
      if(do.abund) {
        abund.subtasks = list()
        for(ProfileID in levels(factor(m_a.abs$attr$ProfileID))) {
          abund.subtask = list()
          abund.subtask$m_a.abs = subset.m_a(m_a.abs,subset=(m_a.abs$attr$ProfileID==ProfileID))
          abund.subtask$prof.descr = sprintf("Community profile %s",ProfileID)
          abund.subtask$do.contrasts=F
          abund.subtask$contrasts.method=NULL
          abund.subtasks = c(abund.subtasks,list(abund.subtask))
        }
        if(do.contrasts) {
          for(contrasts.method in contrasts.methods) {
            abund.subtask = list(m_a.abs=m_a.abs,
                                 prof.descr=paste("Difference between community profiles.",
                                                  switch(contrasts.method,
                                                         diff="Values are subtracted.",
                                                         lfc="Log-fold-change of values."),
                                                  "Even community is used as base."),
                                 do.contrasts=T,
                                 contrasts.method=contrasts.method)
            abund.subtasks = c(abund.subtasks,list(abund.subtask))
          }
        }
        
        for(abund.subtask in abund.subtasks) {
          report$add.header(abund.subtask$prof.descr,
                            section.action="push",
                            report.section=report.section,sub=T)
          benchmark.abund(m_a.abs=abund.subtask$m_a.abs,
                          aggr.type=aggr.type,
                          aggr.descr=aggr.descr,
                          drop.taxa=drop.taxa,
                          true.taxa.only=true.taxa.only,
                          other_cnt=other_cnt,
                          norm.count.task=norm.count.task,
                          constrasts.stage.min.feature.prop=constrasts.stage.min.feature.prop,
                          ground.truth.run.id=ground.truth.run.id,
                          do.contrasts=abund.subtask$do.contrasts,
                          contrasts.method=abund.subtask$contrasts.method,
                          plot.profiles.task=test.counts.task$plot.profiles.task,
                          do.plot.profiles=do.plot.profiles)      
          report$pop.section()
        }
      }
      report$pop.section()
    })
    report$pop.section()
  }
  
  stop.cluster.project(cl)
  
}


## number of cores to use on multicore machines
options(mc.cores=4)
options(boot.ncpus=4)
## parallel backend
options(boot.parallel="snow")
library("BiocParallel")
register(SnowParam(4))


## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"

source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)

## Uncomment next line to install packages needed by MGSAT (!!!comment it out
## in all subsequent runs once the packages have been installed!!!).
## Note: you should also pre-install Pandoc program from http://johnmacfarlane.net/pandoc/
## or using your OS package manager (if running on Linux)

#install_required_packages()

## loads dependency packages (which already must be installed)
load_required_packages()

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

## extra dependencies
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

evalsOptions("width",800)
evalsOptions("height",600)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="atovtchi@jcvi.org",
                       title="Report on benchmarking of 16S annotation pipelines",
                       incremental.save=F)
report.section <- report$reset.section()

benchmark.project()

report$save()
