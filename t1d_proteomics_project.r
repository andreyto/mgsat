read.pieper.t1d <- function(count.file="Original Collapesed APEX (All information).AT.tsv",taxa.level=NULL) {
  
  file_name = count.file
  #data =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  data =read.delim(file_name, header=T,sep="\t",stringsAsFactors=FALSE)
  make.global(data)
  #row.names(meta) = meta$ID
  id.ind = 13
  feat.attr.ind.last = 15
  feat.attr = data[,c(1:feat.attr.ind.last)]
  row.names(feat.attr) = data[,id.ind]
  row.names(data) = data[,id.ind]
  data = data[,-c(1:feat.attr.ind.last)]
  data = t(data)
  make.global(data)
  meta = as.data.frame(t(as.data.frame(strsplit(row.names(data),"[._]"))))
  row.names(meta) = row.names(data)
  names(meta) = c("Group","FamilyID","TechReplRelID")
  meta = within(meta,{
    SubjectID = factor(paste(Group,FamilyID,sep="."))
  })
  meta$SampleID = row.names(meta)
  attr.names = names(meta)
  m_a = list(count=data,attr=meta,feat.attr=feat.attr)
  make.global(m_a)
  return(m_a)
}


## This function must generate a list with analysis tasks

gen.tasks.t1d.prot <- function() {
  
  norm.count.task.drop.other = within(list(), {
    method = "norm.ident"
  })
  
  
  mgsat.proteomics.task.template = within(mgsat.16s.task.template, {
    
    label.base = "prot"
    
    read.data.method=read.pieper.t1d
    
    read.data.task = list(
      count.file=NULL
    )
    
    taxa.levels = c("prot")
    
    count.filter.sample.options=list()
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.genesel=T
      do.stabsel=F
      do.glmer=F
      do.adonis=F
      do.divrich=c()
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.select.samples=c()
      
      count.filter.feature.options=list(min_mean=10)
      #count.filter.feature.options=list()
      
      norm.count.task = norm.count.task.drop.other
      
      plot.profiles.abund.task = within(plot.profiles.abund.task, {
        norm.count.task = norm.count.task.drop.other  
      })    
      
    })
    
  })
  
  
  get.taxa.meta.aggr.base<-function(m_a) { 
    m_a = aggregate.by.meta.data.m_a(m_a,
                                           group_col="SubjectID",
                                           count_aggr=mean)    
    report$add.p(paste("After aggregating by averaging technical replicates per subject:",nrow(m_a$count)))
    
    #gmr = group.mean.ratio(m_a$count,m_a$attr$Group)
    
    #m_a$count = m_a$count[,gmr <= 1/1.5 | gmr >= 1.5]
    
    #report$add.p(paste("After filtering for group mean fold change between cases and controls 
    #                   to be at least 1.5 in either direction,",
    #                   ncol(m_a$count),"features left"))
    
    return (m_a)
  }

  task0 = within( mgsat.proteomics.task.template, {
      
  descr = "All samples, aggregated by SubjectID"
  
  read.data.task = within(read.data.task, {
    count.file="Original Collapesed APEX (All information).AT.tsv"
  })
  
  get.taxa.meta.aggr<-function(m_a) { 
    return (get.taxa.meta.aggr.base(m_a))
  }
  
  })
  
  task1 = within( task0, {
    
  do.summary.meta = F
  
  do.tests = T

  summary.meta.task = within(summary.meta.task, {
    meta.x.vars = NULL
    group.var = main.meta.var
  })

  test.counts.task = within(test.counts.task, {
    
    genesel.task = within(genesel.task, {
      
      genesel.param = within(genesel.param, {
        block.attr = "FamilyID"
        type="paired"
        replicates=10
        maxrank=100
        samp.fold.ratio=0.5
        comp.log.fold.change=T
      })
      
    })    
    
    adonis.task = within(adonis.task, {
      
      dist.metr="euclidean"
      col.trans="standardize"
      data.descr="normalized counts"
      
      tasks = list(
        list(formula.rhs=main.meta.var,
             strata=NULL,
             descr="Association with the patient/control status unpaired"),
        list(formula.rhs=main.meta.var,
             strata="FamilyID",
             descr="Association with the patient/control status paired by family")
      )
      
    })
    
  })
    
})

return (list(task1))
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
                       title="Analysis of T1D proteomics data",
                       out.file.md="report.md",
                       incremental.save=T)


res = proc.project(
  task.generator.method=gen.tasks.t1d.prot
)

report$save()

