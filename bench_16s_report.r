
cbind.m_a <- function(m_a.list,batch.attr,col.match=T) {
  cols.map = unlist(lapply(m_a.list,function(x) colnames(x$count)))
  batch.id.rows = unlist(lapply(m_a.list,function(x) (x$attr[,batch.attr])))
  batch.id.cols = unlist(lapply(m_a.list,function(x) {
    b.attr = unique(x$attr[,batch.attr])
    stopifnot(length(b.attr)<=1)
    rep(b.attr,ncol(x$count))
    }))
  names(cols.map) = paste(cols.map,ifelse(batch.id.cols!="",".",""),batch.id.cols,sep="")
  if(col.match) {
    cols = unique(cols.map)
  }
  else {
    cols.map[] = names(cols.map)
    cols = cols.map
  }
  m_a = foreach(m_a=m_a.list,
          .final=function(m_a.list) {
            m_a = list()
            m_a$count = do.call(rbind,lapply(m_a.list,function(x) {x$count}))
            m_a$attr = do.call(rbind,lapply(m_a.list,function(x) {x$attr}))
            m_a
          }) %do% {
            batch.attr.val = unique(m_a$attr[,batch.attr])
            stopifnot(length(batch.attr.val) <= 1)
            if(batch.attr.val != "") {
              sep="."
            }
            else {
              sep=""
            }
            cols.keys = paste(colnames(m_a$count),batch.attr.val,sep=sep)
            rows = paste(rownames(m_a$count),batch.attr.val,sep=sep)
            count = matrix(0.,
                           nrow=nrow(m_a$count),
                           ncol=length(cols),
                           dimnames=list(rows,cols))
            count[,cols.map[cols.keys]] = m_a$count
            m_a$count = count
            rownames(m_a$attr) = rows
            m_a
  }
  m_a
}

load.ground.thruth.m_a <- function(abund.file,aggr.type=c("genus.abund","genus.otus","otu")) {
  x = read.delim(abund.file,header=T,sep="\t")
  x = x[,c("Organism","Genus","Profile","WGS.Proportion.16S.Gene","ProfileID")]
  names(x)[4] = "Prop"
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
  x = read.delim(data.file,header=T,sep="\t")
  ## extract all fields up to total field plus genus
  x = x[,c(names(x)[seq(match("total",names(x))-1)],"genus")]
  x$genus = gsub("\\s","\\.",x$genus)
  x$OTU = paste(x$OTU,x$genus,sep=".")
  rownames(x) = x$OTU
  x$OTU = NULL
  
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
    unk.mask = grepl("uncultured.*",x$genus,ignore.case=T) | (!grepl("\\S+",x$genus))
    x = x[!unk.mask,]
    rownames(x) = x$genus
    x$genus = NULL
    x = t(x)
  }
  
  meta = load.meta.default(meta.file)
  return (merge.counts.with.meta(x,meta))
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

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="atovtchi@jcvi.org",
                       title="Report on benchmarking of 16S annotation pipelines",
                       incremental.save=F)

meta.file.samples = "refdata/bench_bei_meta.txt"
ground.truth.WGS.file = "refdata/sarah.2015-03-08/derived/bei_abund_v.5.txt"

with(mgsat.16s.task.template,{
runs.files = list(
  list(
    RunID="YAP.01",
    read.data.task=within(read.data.task, {
      taxa.summary.file = "yap/bei/v13/baseline/bd568ff2897e717e1ab6c948f9509344.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      otu.shared.file="yap/bei/v13/baseline/21c401833a6ba4d011974fa46f05e476.files_x1.sorted.0.03.shared"
      cons.taxonomy.file="yap/bei/v13/baseline/ff4e455254bb6759e890e170151a4e37.files_x1.sorted.0.03.cons.taxonomy"
      meta.file=meta.file.samples
    })
  ),
  list(
    RunID="Mothur.01",
    read.data.task=within(read.data.task, {
      taxa.summary.file = "bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
      otu.shared.file="bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
      cons.taxonomy.file="bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
      meta.file=meta.file.samples
    })
  )
)

taxa.level = 6
aggr.type = "genus.abund"
col.match = T
ProfileID = "HM782D"
drop.taxa = c("Nothing") #c("Clostridium_sensu_stricto") #c("Helicobacter") #c("Streptococcus")
vegdist.method = "manhattan"

runs.data = lapply(runs.files,
                   function(run.files) {
                     with(run.files,{
                       print(run.files)
                       m_a = do.call(read.data.method,
                                                   c(
                                                     list(taxa.level=taxa.level),
                                                     read.data.task
                                                   )
                       )
                       m_a$attr$RunID = RunID
                       m_a
                     })
                   }
                   )

m_a = cbind.m_a(runs.data,batch.attr="RunID",col.match=col.match)
make.global(m_a)
m_a.gt = load.ground.thruth.m_a(ground.truth.WGS.file,aggr.type=aggr.type)
make.global(m_a.gt)
m_a.gt$attr$SampleID = m_a.gt$attr$ProfileID
m_a.gt$attr$RunID = "Ground.Truth.WGS"
m_a.gt$attr$IdSfx = m_a.gt$attr$RunID
make.global(m_a.gt)

m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.00005.txt",
                           meta.file=meta.file.samples,
                           aggr.type=aggr.type)
m_a.wl$attr$RunID = "WL.00005"
m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
m_a.wl.1 = m_a.wl
m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.0001.txt",
                       meta.file=meta.file.samples,
                       aggr.type=aggr.type)
m_a.wl$attr$RunID = "WL.0001"
m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
m_a.wl.2 = m_a.wl

m_a$attr$IdSfx = ""
m_a.prop = norm.count.m_a(m_a,method="norm.prop")
make.global(m_a.prop)
m_a.prop = cbind.m_a(list(m_a.prop,m_a.gt,m_a.wl.1,m_a.wl.2),batch.attr="IdSfx",col.match=col.match)
m_a.prop$attr$IdSfx = NULL
m_a.prop = subset.m_a(m_a.prop,subset=(m_a.prop$attr$ProfileID==ProfileID))
m_a.prop = count.filter.m_a(m_a.prop,drop.zero=T)
m_a.prop = subset.m_a(m_a.prop,select.count=!(colnames(m_a.prop$count) %in% drop.taxa))
#m_a.prop = subset.m_a(m_a.prop,select.count=seq(9))
m_a.prop = norm.count.m_a(m_a.prop,method="norm.prop")
#m_a.prop = norm.count.m_a(m_a.prop,method="norm.alr",method.args=list(offset=0,tol=0.005,divcomp=2))
#m_a.prop = subset.m_a(m_a.prop,select.count=colSums(abs(m_a.prop$count))>0)
make.global(m_a.prop)

plot.profiles.task = within(test.counts.task$plot.profiles.task, {
  id.vars.list = list(c("RunID"))
  clade.meta.x.vars=NULL
  do.profile=T
  do.clade.meta=F
  show.profile.task=within(show.profile.task, {
    geoms=c("bar_stacked","bar")
  })
})
do.call(plot.profiles,
        c(list(m_a=m_a.prop,
               feature.order=NULL),
          plot.profiles.task
        )
)

report$add.table(as.matrix(vegdist(m_a.prop$count,method = vegdist.method)),show.row.names=T,
                 caption=sprintf("%s distance",vegdist.method))
report$add.table(as.matrix(vegdist(sqrt(m_a.prop$count),method = "euclidian")/sqrt(2)),show.row.names=T,
                 caption="Hellinger distance")
report$add.table(as.matrix(dist.js(m_a.prop$count)),show.row.names=T,
                 caption="Jensen-Shannon distance")

})

report$save()
