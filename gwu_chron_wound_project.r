
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)

load.meta.gwu_cw <- function(file.name,batch=NULL,aggr.var=NULL) {
  
  meta =read.delim(file.name, header=TRUE,stringsAsFactors=T, sep="\t")
  make.global(meta)
  names(meta) =
    replace.col.names(
      names(meta),
      c("new_id","patient_id","sample_type","Hidradenitis.suppurativa"),
      c("SampleID","PatientID","SampleType","Hidradenitis")
    )
  
  meta$SampleID = paste("S_",meta$SampleID,sep="")
  meta = meta[!duplicated(meta$SampleID),]
  
  to.num.fields = c("WSA.collection","WSA.one.month","Healing.rate.collection","Healing.rate.one.month")
  for (field in to.num.fields) {
    meta[,field] = as.numeric(as.character(meta[,field]))
    meta[,paste(field,"quant",sep=".")] = quantcut.ordered(meta[,field])
  }
  
  meta$PatientID = paste("P",as.character(meta$PatientID),sep="")
  
  meta = arrange(meta,PatientID,age,SampleType)
  
  meta$Diabetes = toupper(as.character(meta$Diabetes))
  meta = meta[meta$Diabetes %in% c("YES","NO"),]
  meta$Diabetes = factor(meta$Diabetes)
  
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

summary.meta.gwu_cw <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr
  
  report$add.printed(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~WP.collection","~WP.one.month","~WP.collection+WP.one.month",
                        "~Hidradenitis+WP.one.month","~PatientID+SampleID",
                        "~Hidradenitis+Region","~SampleType","~Diabetes")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  #   with(meta,{
  #     
  #     report$add(aov(Healing.rate.one.month~Hidradenitis),
  #                caption="ANOVA for rate of healing and Hidradenitis")
  #     report$add(qplot(Hidradenitis,Healing.rate.one.month,geom="violin"),
  #                caption="Violin plot for healing and Hidradenitis")
  #     
  #   })
  
}

## This function must generate a list with analysis tasks

gen.tasks.gwu_cw <- function() {
  
  get.taxa.meta.aggr.base<-function(m_a) { 
    return (m_a)
  }

  gen.ordination.task <- function(ordination.task=list(),main.meta.var,distance="euclidean",method="RDA") {
  ordination.task = within(ordination.task, {
    distance=distance
    ord.tasks = list(
      list(
        ordinate.task=list(
          method=method
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
          method=method,
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
  ordination.task
  }
  
  task0 = within( mgsat.16s.task.template, {
    #DEBUG: 
    taxa.levels = c(2,3,4,6,"otu")
    #taxa.levels = c(2)
    
    descr = "One sample per patient"
    
    main.meta.var = "WP.collection"
    
    read.data.task.yap = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="yap/*.shared"
      cons.taxonomy.file="yap/*.taxonomy"
      taxa.summary.file.otu = NA
    })
    
    read.data.task = within(read.data.task.yap, {
      meta.file="new_id.data.merged.v7.txt" #"simplified-wound-meta-Apr9.txt"
      load.meta.method=load.meta.gwu_cw
      load.meta.options=list()
      
      otu.count.filter.options=list()
      count.filter.options = list()
      
    })
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      meta = m_a$attr
      #take unique records
      meta = meta[meta$Hidradenitis!="Yes",,drop=F]
      meta = meta[meta$SampleType=="q",,drop=F]
      meta = meta[!duplicated(meta$PatientID),,drop=F]  
      m_a = subset.m_a(m_a,subset = rownames(meta))
      return (m_a)
    }
    
    summary.meta.method=summary.meta.gwu_cw
    
    test.counts.task = within(test.counts.task, {
      
      do.extra.method = c()
      do.ordination = T
      
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
        dist.metr="euclidean"
        trans.show=NULL
        stand.show="range"
        cluster.row.cuth=10
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        hmap.width=1000
        hmap.height=hmap.width*0.8
        attr.annot.names=c(main.meta.var,"Year","Region","Diabetes",
                           "WP.one.month","Hidradenitis",
                           "SampleType","Healing.rate.collection.quant")
        clustering_distance_rows="pearson"
        km.abund=0
        km.diversity=3
        show_row_names=T
      })
      
      ordination.task = gen.ordination.task(ordination.task,main.meta.var=main.meta.var,
                                            distance="euclidean")
            
    })
    
  })
  
  task1 = within( task0, {
        
    do.summary.meta = T
    
    do.tests = T
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = "Healing.rate.collection"
      group.vars = c(main.meta.var)
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = T
      do.glmer = F
      do.network.features.combined=T
      #do.divrich = c(6,"otu")
      
      do.plot.profiles.abund=F
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        #n.rar.rep=4
        is.raw.count.data=T
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = paste("Diabetes",main.meta.var,sep="*")
        })      
        do.plot.profiles=T
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = main.meta.var
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata=NULL,
               descr="Association with the healing status")
        )
        
      })
      
      glmer.task = within(glmer.task, {
        
        tasks = list(list(
          descr.extra = "",
          formula.rhs = paste(main.meta.var,"(1|PatientID/SampleID)",sep="+"),
          linfct=c("WP.one.monthH = 0")
        ))
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"SampleType"),c(main.meta.var,"Diabetes"))
        feature.meta.x.vars=c("Healing.rate.collection","age")
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"SampleType","Diabetes")
      })
      
    })
    
  })
  
  task1.with_repeats = within( task1, {

    descr = "All samples per patient"

    taxa.levels = c(6)
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      return (m_a)
    }    
    
    do.summary.meta = T
    
    do.tests = T
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = "Healing.rate.collection"
      group.vars = c(main.meta.var)
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = T
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.network.features.combined=F
      do.divrich = c() #c(6,"otu")
      
      do.plot.profiles.abund=F
      do.heatmap.abund=T

      norm.count.task = within(norm.count.task, {
        method="norm.prop"
        method.args = list()
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        #method="norm.clr"
        #method.args=list()
      })      
      
      heatmap.combined.task = within(heatmap.combined.task, {
        km.abund=0
        km.diversity=0
        clustering_distance_rows="manhattan"
        })

      ordination.task = gen.ordination.task(ordination.task,main.meta.var=main.meta.var,
                                            distance="manhattan",method="NMDS")
      
    })
    
  })
  
  task2 = within( task0, {
    
    do.summary.meta = F
    
    do.tests = T
    
    descr = "One sample per patient, association with healing rate"
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset = (!is.na(m_a$attr$Healing.rate.collection)))
      return (m_a)
    }
    
    main.meta.var = "Healing.rate.collection.quant"
    main.meta.var.cont = "Healing.rate.collection" 
    
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = c("age")
      group.var = c(main.meta.var)
    })
    
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = T
      do.adonis = T
      do.genesel = F
      do.stabsel = T
      do.glmer = F
      do.extra.method = c()
      #do.divrich = c()
      
      do.plot.profiles.abund=T
      do.heatmap.abund=T

      norm.count.task = within(norm.count.task, {
        method="norm.ihs.prop"
        method.args = list(theta=2000)
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        #method="norm.clr"
        #method.args=list()
      })
      
      divrich.task = within(divrich.task,{
        #n.rar.rep=4
        is.raw.count.data=T
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var.cont
        })      
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = main.meta.var.cont
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var.cont
        args.fitfun = within(args.fitfun, {
          family="gaussian"
        })
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var.cont,
               strata=NULL,
               descr=sprintf("Association with %s",main.meta.var.cont))
          #           ,
          #           list(formula.rhs=paste(main.meta.var.cont,"Hidradenitis",sep="*"),
          #                strata=NULL,
          #                descr=sprintf("Association with the %s and Hidradenitis status",main.meta.var.cont)) 
        )
      })
      
      ordination.task = gen.ordination.task(ordination.task,main.meta.var=main.meta.var.cont,
                                            distance="euclidean")
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var))
        feature.meta.x.vars=c("age")
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var.cont)
      })
      
    })
    
  })
  return (list(task1.with_repeats))
  return (list(task1,task1.with_repeats,task2))
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
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=F)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="atovtchi@jcvi.org",
                       title="Analysis of 16S data for chronic wound microbiome",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.gwu_cw
)

report$save()

