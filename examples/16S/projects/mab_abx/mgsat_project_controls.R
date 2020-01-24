
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load_meta_controls <- function(file.name) {
    return (NULL)
}
  
## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.controls <- function(m_a) {
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")

}


new_ordination.task <- function(main.meta.var,norm.method,label=NULL,size=NULL) {
  if(norm.method=="prop") {
    ord.method = "CCA"
    distance.0="bray"
  }
  else {
    ord.method = "RDA"
    distance.0="euclidean"
  }
  within(mgsat.16s.task.template$test.counts.task$ordination.task, {
    distance=distance.0
    ord.tasks = list(
      list(
        ordinate.task=list(
          method=ord.method
          ##other arguments to phyloseq:::ordinate
        ),
        plot.task=list(
          type="samples",
          color=main.meta.var,
          label = label,
          size = size
          ##other arguments to phyloseq:::plot_ordination
        )
      ),
      list(
        ordinate.task=list(
          method=ord.method,
          formula=main.meta.var
          ##other arguments to phyloseq:::ordinate
        ),
        plot.task=list(
          type="samples",
          color=main.meta.var,
          label = label,
          size = size
          ##other arguments to phyloseq:::plot_ordination
        )
      )          
    )
  })            
}


## This function must generate a list with analysis tasks

gen.tasks.controls <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    taxa.levels = c(2,5,6,"otu")
    #taxa.levels = c(6)
    
    descr = "All samples, no aggregation, no tests here, only plots"
    
    main.meta.var = "Group"    
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="otu.shared"
      cons.taxonomy.file="cons.taxonomy"
      meta.file="StudyDesign.txt"
      load.meta.method=load_meta_controls
      load.meta.options=list()
      count.filter.options = list()    
      otu.count.filter.options=list()      
    })
    
    count.filter.sample.options=list() 
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      m_a$attr$Group = substr(as.character(m_a$attr$Group),1,2)
      m_a$attr$SampleType = substr(as.character(m_a$attr$Group),1,1)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Group %in% c("M1","M2","NT")))
      return(m_a)
    }
    
    summary.meta.method=summary.meta.controls
    
    test.counts.task = within(test.counts.task, {
      
      do.ordination=T
      do.network.features.combined=T      
      
      count.filter.feature.options = within(list(), {
        min_quant_mean_frac=0
        min_quant_incidence_frac=0
        #min_max=30
        min_mean=0
      })
      
      norm.count.task = within(norm.count.task, {
        method="norm.prop"
        method.args=list()
        #method="norm.ihs.prop"
        #method.args = list(theta=1)
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        #method="norm.clr"
        #method.args=list()
      })
      
      adonis.task = within(adonis.task, {
        
        dist.metr="bray"
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
      
      heatmap.combined.task = within(heatmap.combined.task, {
        hmap.width=1000
        hmap.height=hmap.width*0.8
        attr.annot.names=c(main.meta.var)
        clustering_distance_rows="pearson"
        km.abund=0
        km.diversity=0
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="ihs.prop",label="SampleID")      
      
    })
    
  })
  
  
  task1 = within( task0, {

    descr = "All samples, no aggregation, no tests, only plots"    
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      return(m_a)
    }    
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = NULL
      group.vars = c("Group")
    })
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c("otu")
      do.deseq2 = F
      do.adonis = T
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })        
        do.plot.profiles = T
        n.rar.rep=40
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(main.meta.var)
        feature.meta.x.vars=c()
        do.profile=T
        do.feature.meta=F
        #show.profile.task=within(show.profile.task,{
        #  facet_wrap_ncol=3
        #})
        
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               #strata="MouseID",
               ## strata will not work - mice are nested within treatments
               ## also, there are missing samples between time points
               descr="Association with group")
        )
        
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=main.meta.var
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="prop",label="SampleID") 
      
    })
    
  })

  
  return (list(task1))
}

threads=16

## number of cores to use on multicore machines
options(mc.cores=threads)
options(boot.ncpus=threads)
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
register(SnowParam(threads))

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=F)

report <- PandocAT$new(author="John Doe",
                       title="Analysis of control samples",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.controls
)

report$save()

