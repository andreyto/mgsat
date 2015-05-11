
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load.meta.diet <- function(file.name) {
  ## We are calling default simple table loading here,
  ## but you can do something else. Result must have
  ## rownames() set to SampleID field
    meta = load.meta.default(file.name)
    ## we need to create this as ordered quantile to show
    ## in the right order on plots (and possible used
    ## in models that care about ordered factors)
    meta$age.quant = quantcut.ordered(meta$age)
    return (meta)
}
  

## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.diet <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~Sample.type+DietStatus","~Drug.Before.Diet + Sample.type",
                        "~Complaints + Sample.type",
                        "~Sample.type+visit","~MatchedGroupID","~Sample.type.1","~SubjectID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(meta,{
    report$add.printed(summary(aov(age~Sample.type)),
                       caption="ANOVA for age and sample type")
    report$add(qplot(Sample.type,age,geom="violin"),
               caption="Violin plot for age and sample type")
  })
  
  with(meta,{
    report$add(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit")
    
  })
  with(meta[meta$Sample.type=="patient",],{
    report$add(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit, patients only")
    
  })
  
  #summary(glht(lm(age~visit,data=meta[meta$Sample.type=="patient",]),linfct="visit=0"))
  #summary(glht(lmer(age~visit+(visit|Sample.type),data=meta),linfct="visit=0"))
  report$add(ggplot(meta,aes(x=visit,y=age,color=Sample.type))+
               geom_point()+
               stat_smooth(method="loess", se = T,degree=1,size=1),
             caption="Plot for age and visit with Loess trend line")
  
}

## This function must generate a list with analysis tasks

gen.tasks.diet <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    taxa.levels = c(2,4,6,"otu")
    #taxa.levels = c(2)
    
    descr = "All samples, no aggregation, no tests here, only plots"
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = "example_taxa.summary.file"
      otu.shared.file="example_otu.shared.file"
      cons.taxonomy.file="example_cons.taxonomy.file"
      meta.file="example_meta.txt"
      load.meta.method=load.meta.diet
      load.meta.options=list()
      count.filter.options = list()    
      otu.count.filter.options=list()      
    })
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      ## We already have all fields, so not aggregating anything here.
      
      ##Throw away all visits with very few samples
      m_a = subset.m_a(m_a,subset=(m_a$attr$visit<=4))
      return(m_a)
    }
    
    summary.meta.method=summary.meta.diet
    
    test.counts.task = within(test.counts.task, {
      
      count.filter.feature.options = within(list(), {
        min_quant_mean_frac=0.25
        min_quant_incidence_frac=0.25
        #min_max=30
        min_mean=10
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
      meta.x.vars = c("visit")
      group.vars = c("Sample.type","visit")
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = NULL
        counts.glm.task = NULL
        do.plot.profiles = T
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Sample.type","visit"),
                            c("Sample.type.Drug.Before","visit"),
                            c("Drug.Before.Diet","Sample.type.1"))
        feature.meta.x.vars=c("visit")
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c("Sample.type","visit","Drug.Before.Diet")
      })
      
    })
    
  })

  
  task2 = within( task0, {
    
    main.meta.var = "Sample.type"
    
    descr = "Patient/control samples before diet aggregated by SubjectID"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$DietStatus=="before.diet"))
      m_a = aggregate.by.meta.data.m_a(m_a,group_col="SubjectID")
      return(m_a)
    }

    summary.meta.task = within(summary.meta.task, {
      group.vars = c(main.meta.var)
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = T
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })
        do.plot.profiles = T
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
               descr="Association with the patient/control status unpaired"),
          list(formula.rhs=paste("age.quant *", main.meta.var),
               strata=NULL,
               descr="Association with the age quartiles and patient/control status")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"age.quant"))
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"age.quant")
      })
      
    })
    
  })

  
  task2.1 = within( task2, {
    
    descr = paste(descr,"Additional tests")
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = task2$get.taxa.meta.aggr(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$has.matched.subject)) 
      return(m_a)
    }    
    
    test.counts.task = within(test.counts.task, {
      
      do.divrich = c()
      do.deseq2 = F
      do.adonis = T
      do.genesel = T
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.extra.method = taxa.levels

      genesel.task = within(genesel.task, {
        genesel.param = within(genesel.param, {
          block.attr = "MatchedGroupID"
          type="paired"
          #replicates=0
        })
      })

      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata="MatchedGroupID",
               descr="Association with the patient/control status paired by family")
        )
        
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="Sample.type",
                                        block.attr="MatchedGroupID",
                                        n.perm=4000,
                                        dist.metr="euclidian",
                                        norm.count.task=norm.count.task.extra
          )
        }
        
        norm.count.task.extra = within(norm.count.task, {
        })
        
      })
      
      
    })
    
  })
  
  
  task3 = within( task0, {
    
    main.meta.var = "DietStatus"
    
    descr = "Patients' samples at visits 1 (before diet) and 2 (after diet), only paired samples"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient" 
                                   & m_a$attr$visit <= 2 
                                   & m_a$attr$visit.1
                                   & m_a$attr$visit.2))
      return(m_a)
    }
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c()
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = T
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      do.extra.method = taxa.levels
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })
        do.plot.profiles = T
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = main.meta.var
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
        do.plot.profiles = T
        genesel.param = within(genesel.param, {
          block.attr = "SubjectID"
          type="paired"
          #replicates=0
        })
      })
      
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               descr="Association with diet status unpaired by subject"),
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with diet status paired by subject"),
          list(formula.rhs=paste("Drug.Before.Diet * ", main.meta.var),
               strata=NULL,
               descr="Association with Drug use before diet and diet status")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"Drug.Before.Diet"))
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Drug.Before.Diet")
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="DietStatus",
                                        block.attr="SubjectID",
                                        n.perm=8000,
                                        dist.metr="euclidian",
                                        col.trans="ident",
                                        norm.count.task=norm.count.task.extra
                                        )
        }
        norm.count.task.extra = within(norm.count.task, {
        })
        
      })
      
    })
    
  })

  task3.1 = within( task3, {
    
    descr = paste(descr,"Additional tests")
    
    do.summary.meta = F
    
    do.tests = T
        
    test.counts.task = within(test.counts.task, {
      
      do.divrich = c()
      do.deseq2 = T
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.extra.method = c()
            
      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("Drug.Before.Diet+%s",main.meta.var)
      })
    })
    
  })
  
  
  task4 = within( task0, {
    
    main.meta.var = "visit"
    
    descr = "Patients samples"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient"))
      return(m_a)
    }
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c()
      do.deseq2 = T
      do.adonis = T
      do.genesel = F
      do.stabsel = T
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })
        do.plot.profiles = T
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = main.meta.var
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
        args.fitfun = within(args.fitfun, {
          family="gaussian"
          standardize=T                                     
        })
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with visit paired by subject")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var,"Drug.Before.Diet"))
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Drug.Before.Diet")
      })
      
    })
    
  })
  
  #return (list(task2.1))
  return (list(task1,task2,task2.1,task3,task3.1,task4))
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

report <- PandocAT$new(author="noone@mail.com",
                       title="Analysis of Dieting study 16S data",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.diet
)

report$save()

