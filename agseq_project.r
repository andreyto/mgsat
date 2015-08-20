
read.agseq.jincheng <- function(count.file,
                                min.count=0,
                                to.ranks=T) {
  require(data.table)
  
  data = fread(count.file, header=F, sep="\t", stringsAsFactors=FALSE, data.table=T)
  setnames(data,c("seq","count"))
  if(min.count>0) {
    data = data[count>=min.count,,drop=F]
  }
  data = data[order(data[,count],decreasing=T),]
  if(to.ranks) {
    new.seq = sprintf("F%08i",1:nrow(data))
    data$seq = new.seq
  }
  
  data
}

agseq.count.file.name.parser.jincheng <- function(count.file) {
  ##Hs638-48op-K-L_S6_L001_R2.CDR.H3.prot.count
  require(stringr)
  m = str_match(count.file,"([[:alnum:]]+)-([[:digit:]]+)([[:alpha:]]+)-(.*)\\.(CDR.*)\\.prot\\.count")
  m = as.data.frame(m[1,-1,drop=F])
  names(m) = c("SeqRunID","Time","Reaction","ReplicateID","Region")
  m$Time = as.numeric(as.character(m$Time))
  m
}

read.agseq.study.jincheng <- function(count.files,
                                      min.count=0,
                                to.ranks=T) {
  count.files = Sys.glob(count.files)
  attr = list()
  x = list()
  f_num = 1
  for(count.file in count.files) {
    d = read.agseq.jincheng(count.file,min.count=min.count,to.ranks=to.ranks)
    SampleID = sprintf("S%04i",f_num)
    d$SampleID = SampleID
    x = c(x,list(d))
    #make.global(x)
    #if(f_num>1) stop("DEBUG")
    a = agseq.count.file.name.parser.jincheng(basename(count.file))
    a$SampleID = SampleID
    attr = c(attr,list(a))
    f_num = f_num + 1
  }
  x = rbindlist(x)
  attr = rbindlist(attr)
  rownames(attr) = attr$SampleID
  return(list(data=x,attr=attr))
}

read.agseq.study.jincheng.m_a <- function(count.files,
                                      min.count=0,
                                      to.ranks=T,
                                      taxa.level="seq") {
  r = read.agseq.study.jincheng(count.files=count.files,
                                min.count=min.count,
                                to.ranks=to.ranks) 
  d = dcast.data.table(r$data,SampleID~seq,value.var = "count", fill=0)
  count = as.matrix(d[,colnames(d) != "SampleID",with=F])
  rownames(count) = d$SampleID
  attr = as.data.frame(r$attr)
  rownames(attr) = attr$SampleID
  m_a = list(count=count,attr=attr)
  make.global(m_a)
  stopifnot(all(rownames(m_a$count)==rownames(m_a$attr)))
  m_a
}

gen.tasks.agseq <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    label.base = "agseq"
    
    taxa.levels = c("seq")

    descr = "All samples"
    
    main.meta.var = "Reaction"    
    
    read.data.method=read.agseq.study.jincheng.m_a
    
    read.data.task = within(list(), {
      count.files = "input/Bioinformatics results_2015-08-20/count/Hs717*.prot.count"
      min.count = 2
    })
    
    summary.meta.method=summary.meta.diet
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.genesel=F
      do.stabsel=F
      do.glmer=F
      do.adonis=F
      do.divrich=taxa.levels
      do.divrich.pre.filter=T
      do.divrich.post.filter=F    
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.ordination=F
      do.network.features.combined=F
      do.select.samples=c()
      do.extra.method=c()
      do.aggr.after.norm=c()

      feature.ranking = "stabsel"
      
      alpha = 0.05
      
      divrich.task = within(list(),{
        n.rar.rep=4
        is.raw.count.data=T
        filtered.singletons=T
        group.attr = main.meta.var
        counts.glm.task = within(list(),{
          formula.rhs = main.meta.var
        })
        beta.task = within(list(),{
          method="-1"
          betadisper.task=list()
          adonis.task=NULL #will be replaced by global adonis.task
        })
        ## Same structure as the task-wide genesel.task; if NULL,
        ## this will be taken from task-wide structure
        counts.genesel.task = NULL
        do.plot.profiles = T
        ## Computing beta-diversity matrix on multiple rarefications can take a while
        do.beta = F
        do.accum = F
        do.incidence = F
      })
      
      count.filter.feature.options = within(list(), {
        min_quant_mean_frac=0.25
        min_quant_incidence_frac=0.25
        #min_max=30
        min_mean=10
      })
      
      norm.count.task = within(norm.count.task, {
        method="norm.ihs.prop"
        method.args = list(theta=2000)
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        #method="norm.clr"
        #method.args=list()
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        hmap.width=1000
        hmap.height=hmap.width*0.8
        attr.annot.names=c(main.meta.var)
        clustering_distance_rows="pearson"
        km.abund=0
        km.diversity=0
      })

      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"Region"))
        do.profile=T
        do.feature.meta=F
      })
      
    })
  })
  return (list(task0))
}


## number of cores to use on multicore machines
options(mc.cores=4)
options(boot.ncpus=4)
## parallel backend
options(boot.parallel="snow")


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

library("BiocParallel")
register(SnowParam(4))

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

report <- PandocAT$new(author="noone@mail.com",
                       title="Analysis of AgSeq annotation data",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.agseq
)

report$save()

