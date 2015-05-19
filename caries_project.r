
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)

load.meta.caries <- function(file.name,batch=NULL,aggr.var=NULL) {
  
  meta =read.delim(file.name, header=TRUE,stringsAsFactors=T, sep="\t")
  
  allnames = replace.col.names(names(meta),
                               c("zygosity","jcvi_sex","sample_site","collection_date","jcvi_twin_id"),
                               c("Zygosity","Gender","OralSite","CollectionDate","TwinIndex"))
  
  names(meta) = allnames  
  
  meta = within(meta,{
    OralSite = factor(substr(OralSite,1,2))
  })
  
  row.names(meta) = meta$SampleID
  
  ##pass filter argument and filter here before aggregation 
  
  if(!is.null(aggr.var)) {
    
    meta = aggregate.by.meta.data(meta,
                                  group_col=aggr.var,
                                  col_ignore=names(meta))  
  }
  
  return (meta)
}


## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.caries <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr
  
  report$add.printed(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~OralSite","~Gender","~Zygosity",
                        "~Gender+Zygosity")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(meta,{
    
    report$add(aov(age~OralSite),
               caption="ANOVA for age and collection site")
    report$add(qplot(OralSite,age,geom="violin"),
               caption="Violin plot for age and collection site")
    
  })
  
}

## This function must generate a lits with analysis tasks

gen.tasks.caries <- function() {
  
  get.taxa.meta.aggr.base<-function(m_a) { 
    return (m_a)
  }
  
  task0 = within( mgsat.16s.task.template, {
    #DEBUG: 
    taxa.levels = c(2,3,4,5,6,"otu")
    #taxa.levels = c(2)
    
    descr = "All samples"
    
    main.meta.var = "OralSite"
    
    read.data.task.yap = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="input/EDITED*.shared"
      cons.taxonomy.file="input/EDITED*AT.taxonomy"
      taxa.summary.file.otu = NA
    })
    
    read.data.task = within(read.data.task.yap, {
      meta.file="input/EDITED_metadata_tab.txt"
      load.meta.method=load.meta.caries
      load.meta.options=list()
      
      otu.count.filter.options=list()
      count.filter.options = list()
      
    })
    
    get.taxa.meta.aggr<-function(taxa.meta) { 
      return (get.taxa.meta.aggr.base(taxa.meta))
    }
    
    summary.meta.method=summary.meta.caries
    
    test.counts.task = within(test.counts.task, {
      
      count.filter.feature.options = within(count.filter.feature.options, {
        #min_mean_frac=0.00005
        min_quant_mean_frac=0.25
        min_quant_incidence_frac=0.25
        min_max=20
        #min_mean=10
      })    
      
      
      norm.count.task = within(norm.count.task, {
        #method="norm.ihs.prop"
        #method.args = list(theta=2000)
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        method="norm.clr"
        method.args=list()
      })
      
      adonis.task = within(adonis.task, {
        
        dist.metr="euclidean"
        col.trans=NULL
        norm.count.task=NULL
        data.descr="normalized counts"
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        trans.clust=NULL
        stand.clust=NULL
        dist.metr="euclidian"
        trans.show=NULL
        stand.show="range"
        cluster.row.cuth=10
      })
      
      ordination.task = within(ordination.task, {
        distance="euclidean"
        ord.tasks = list(
          list(
            ordinate.task=list(
              method="RDA"
              ##other arguments to phyloseq:::ordinate
            ),
            plot.task=list(
              type="samples",
              color=main.meta.var
              ##other arguments to phyloseq:::plot_ordination
            )
          ),
          list(
            ordinate.task=list(
              method="RDA",
              formula=main.meta.var
              ##other arguments to phyloseq:::ordinate
            ),
            plot.task=list(
              type="samples",
              color=main.meta.var
              ##other arguments to phyloseq:::plot_ordination
            )
          )          
        )
      })
            
    })
    
  })

  task1 = within( task0, {
    
    descr = "All samples, general review"
    
    do.summary.meta = T
    
    do.tests = T
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = "age"
      group.vars = c(main.meta.var)
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.divrich = c(6,"otu")
      
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      do.ordination=T
      
      divrich.task = within(divrich.task,{
        #n.rar.rep=4
        is.raw.count.data=T
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = NULL
        })      
        do.plot.profiles=T
      })      
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(),c(main.meta.var),c("Zygosity"),c(main.meta.var,"Zygosity"))
        feature.meta.x.vars=c("CollectionDate","age")
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"age")
      })
      
    })
    
  })
  
  
  task1.1 = within( task1, {
    
    descr = "All samples, aggregated by SubjectID"
    
    do.summary.meta = F
    
    do.tests = T

    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = aggregate.by.meta.data.m_a(m_a,group_col="SubjectID")
      m_a$OralSite = factor(m_a$OralSite)
      return(m_a)
    }    
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      
      do.divrich = c()
      
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      
      do.extra.method = taxa.levels
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests) {
          test.dist.matr.within.between(m_a=m_a.norm,
                                        group.attr="TwinIndex",
                                        block.attr="FamilyID",
                                        col.trans=NULL,
                                        dist.metr="euclidian",
                                        n.perm=4000)
        }
      })
      
    })
    
  })

  task1.2 = within( task1, {
    
    descr = "Only samples with defined mouth quarters"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(substr(m_a$attr$OralSite,1,1) == "Q"))
      return(m_a)
    }    
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = T
      do.glmer = F
      
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      
      do.extra.method = taxa.levels
      
      divrich.task = within(divrich.task,{
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })      
        do.plot.profiles = F
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("SubjectID+%s",main.meta.var)
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
        genesel.param = within(genesel.param, {
          block.attr = "SubjectID"
          type="paired"
        })        
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with the site")
        )
        
      })
      
      glmer.task = within(glmer.task, {
        
        tasks = list(list(
          descr.extra = "",
          formula.rhs = paste(main.meta.var,"(1|SubjectID/SampleID)",sep="+"),
          linfct=c("OralSiteQ4 = 0")
        ))
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests) {
          test.dist.matr.within.between(m_a=m_a.norm,
                                        group.attr=main.meta.var,
                                        block.attr="SubjectID",
                                        col.trans=NULL,
                                        dist.metr="euclidian",
                                        n.perm=4000)
        }
      })
      
    })
    
  })
  
  
  return (list(task1,task1.1,task1.2))
  #return (list(task1))
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
set_trace_options(try.debug=F)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="atovtchi@jcvi.org",
                       title="Analysis of 16S data for oral microbiome",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.caries
)

report$save()

