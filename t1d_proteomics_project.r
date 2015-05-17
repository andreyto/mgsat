load.meta.t1d.pieper <- function(file.name) {
  
  meta =read.delim(file.name, header=TRUE,stringsAsFactors=T, sep="\t")
  make.global(meta)
  meta$ProteomicsID = factor(sub("-",".",as.character(meta$ProteomicsID),fixed=T))
  rownames(meta) = meta$ProteomicsID
  ##Setting Control to base level plays nice with DESeq2 defaults
  meta$T1D = relevel(meta$T1D,"Control")
  return (meta)
}

read.pieper.t1d <- function(count.file="Original Collapesed APEX (All information).AT.tsv",
                            meta.file="metadata.Family_groups_T1D-C.txt",
                            load.meta.method=NULL,
                            taxa.level=NULL) {
  require(data.table)
  
  file_name = count.file
  #data = read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  #data = read.delim(file_name, header=T,sep="\t",stringsAsFactors=FALSE)
  data = fread(file_name, header=T, sep="\t", stringsAsFactors=FALSE, data.table=F)
  
  #row.names(meta) = meta$ID
  id.ind = 13
  feat.attr.ind.last = 15
  feat.attr = data[,c(1:feat.attr.ind.last)]
  row.names(feat.attr) = data[,id.ind]
  row.names(data) = data[,id.ind]
  data = data[,-c(1:feat.attr.ind.last)]
  data = t(data)
  make.global(data)
  meta = as.data.frame(t(as.data.frame(strsplit(row.names(data),"[-_]"))))
  make.global(meta)
  row.names(meta) = row.names(data)
  names(meta) = c("Group","DummyID","TechReplRelID")
  meta = within(meta,{
    ProteomicsID = factor(paste(Group,DummyID,sep="."))
  })
  meta$SampleID = row.names(meta)
  meta.bio = do.call(load.meta.method,
                     list(file.name=meta.file)
                     )
  meta = join(meta,meta.bio,by="ProteomicsID",match="first")
  meta = within(meta,{
    UnitID = factor(paste(FamilyID,Group,sep="."))
  })  
  m_a = list(count=data,attr=meta,feat.attr=feat.attr)
  make.global(m_a)
  return(m_a)
}

## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.t1d.proteomics <- function(m_a) {
  require(doBy)
  
  report$add.header("Summary of metadata variables")
  
  meta = m_a$attr
  
  report$add.printed(summary(meta),caption="Summary of metadata variables")
  
  #xtabs.formulas = list("Group~","Group~FamilyID","FamilyID~.")
  xtabs.formulas = list("TechReplRelID~FamilyID+Group+SubjectID",
                        "SubjectID~FamilyID+Group",
                        "Group~FamilyID"
  )
  fact.xtabs = meta
  #xtabs.formulas = list("~SubjectID","~FamilyID+SubjectID","~SubjectID+TechReplRelID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = summaryBy(as.formula(xtabs.formula),data=fact.xtabs,FUN=length,keep.names=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    #report$add.printed(summary(fact.xtabs))
  }
  
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
      count.file=NULL,
      meta.file=NULL,
      load.meta.method=load.meta.t1d.pieper
    )
    
    summary.meta.method=summary.meta.t1d.proteomics
    
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
      
      feature.ranking = "genesel"
      
      count.filter.feature.options=list(min_incidence_frac=0.25)
      #count.filter.feature.options=list()
      
      norm.count.task = norm.count.task.drop.other
      
      plot.profiles.abund.task = within(plot.profiles.abund.task, {
        norm.count.task = norm.count.task.drop.other  
      })    
      
    })
    
  })
  
  
  get.taxa.meta.aggr.base<-function(m_a) { 
    m_a = aggregate.by.meta.data.m_a(m_a,
                                     group_col="UnitID",
                                     count_aggr=mean)    
    report$add.p(paste("After aggregating by averaging samples per family/condition:",nrow(m_a$count)))
    
    return (m_a)
  }
  
  task0 = within( mgsat.proteomics.task.template, {
    
    descr = "All samples, aggregated by UnitID (family/condition)"
    
    read.data.task = within(read.data.task, {
      count.file="Original Collapesed APEX (All information).AT.tsv"
      meta.file="metadata.Family_groups_T1D-C.txt"
    })
    
    get.taxa.meta.aggr<-function(m_a) { 
      return (get.taxa.meta.aggr.base(m_a))
    }
    
  })
  
  task1 = within( task0, {
    
    descr = "All samples, aggregated by UnitID (family/condition), paired Wilcoxon test"
    
    do.summary.meta = F
    
    do.tests = T
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = NULL
      group.vars = main.meta.var
    })    
    
    test.counts.task = within(test.counts.task, {
      
      do.stabsel=F
      do.adonis=F
      
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      
      genesel.task = within(genesel.task, {
        
        genesel.param = within(genesel.param, {
          block.attr = "FamilyID"
          type="paired"
          replicates=400
          maxrank=20
          samp.fold.ratio=0.5
          comp.log.fold.change=T
        })
        
      })    
      
      adonis.task = within(adonis.task, {
        
        #dist.metr="euclidean"
        #col.trans="standardize"
        norm.count.task=NULL
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
  
  task2 = within( task1, {
    
    descr = "All samples, aggregated by UnitID (family/condition), unpaired Wilcoxon test"
    
    do.summary.meta = T
    
    do.tests = T
    
    test.counts.task = within(test.counts.task, {  
      
      do.genesel=T
      do.stabsel=T
      do.adonis=T
      
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      do.extra.method = taxa.levels
      
      genesel.task = within(genesel.task, {
        
        genesel.param = within(genesel.param, {
          type="unpaired"
        })
        
      })
      
      extra.method.task = within(extra.method.task, {
      
        func = function(m_a,m_a.norm,res.tests) {
          test.dist.matr.within.between(m_a=m_a.norm,
                                        group.attr="Group",
                                        block.attr="FamilyID",
                                        n.perm=4000)
        }
      })
    })
    
  })

  task.verify.biomarkers.power = within( task1, {
    
    descr = "All samples, aggregated by UnitID (family/condition), power analysis for biomarker verification study"
    
    do.summary.meta = F
    
    do.tests = T
    
    test.counts.task = within(test.counts.task, {  
      
      do.genesel=F
      do.stabsel=F
      do.adonis=F
      
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      
      do.extra.method = taxa.levels
            
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,id.markers) {
          verification.power(m_a=m_a.norm,
                                        group.attr="Group",
                                        id.markers=id.markers)
        }
        
        id.markers = c("P17050","O00754","P53634",
                       "P20774","Q9BTY2","P02750",
                       "Q14393","P40197","Q96HD9",
                       "P14151","Q13421","Q6UXB8",
                       "P19320","Q9BYF1","P08195",
                       "Q13740","P33151","P04066",
                       "P13473","Q9H0E2")
        
      })
    })
    
  })

  task.assoc.power = within( task1, {
    
    descr = "All samples, aggregated by UnitID (family/condition), power analysis for association study"
    
    do.summary.meta = F
    
    do.tests = T
    
    test.counts.task = within(test.counts.task, {  
      
      do.genesel=F
      do.stabsel=F
      do.adonis=F
      
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      
      do.extra.method = taxa.levels
      
      count.filter.feature.options=list(min_incidence_frac=0.25)
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,id.markers) {
          wilcox.power(m_a=m_a.norm,
                             group.attr="Group",
                             id.markers=id.markers,R=2000,n=200)
        }
        
        id.markers = c("P17050","O00754","P53634",
                       "P20774","Q9BTY2","P02750",
                       "Q14393","P40197","Q96HD9",
                       "P14151","Q13421","Q6UXB8",
                       "P19320","Q9BYF1","P08195",
                       "Q13740","P33151","P04066",
                       "P13473","Q9H0E2")
        
      })
    })
    
  })  

  return (list(task.assoc.power))
  return (list(task.verify.biomarkers.power))
  return (list(task1,task2))
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
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.t1d.prot
)

report$save()

