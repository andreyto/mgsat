
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load.meta.choc <- function(file.name) {
  meta = read.delim(file.name,header=T,stringsAsFactors=T)
  make.global(meta)
  
  allnames = replace.col.names(names(meta),
                               c("Subject.ID..blinded.","JCVI.Sample.ID","Sample_Tube_ID..revised.","Subject.s.Gender","Subject.s.YEAR.of.birth","Had.subject.fever.at.the.time.of.sample.collection."),
                               c("SubjectID","SampleID","SampleID.compound","Gender","YearOfBirth","Fever.descr"))
  
  names(meta) = allnames
  
  meta$SampleID = factor(paste("CHOC",meta$SampleID,sep=""))
  
  ## Format of SubjectID from Raja:
  ## The third letter in the letter set describes whether the sample is 
  ## from the patient (P) or sibling (S). The first set of numbers is 
  ## the date of consent without year. The second set of two numbers 
  ## links patient and sibling together. For example, ARP-0328-01 and 
  ## CRS-0329-01 are family set 01
  
  meta$FamilyID = as.factor(substring(meta$SubjectID,10))
  
  meta$Sample.type = gsub(" ",".",meta$Sample.type)
  
  #Therapy.Status
  meta$Sample.type.1 = meta$Sample.type
  
  meta$Sample.type = as.factor(unlist(apply(meta,
                                            1,
                                            function(row) {switch(row["Sample.type.1"],
                                                                  sibling ="sibling",
                                                                  patient.after.chemo="patient",
                                                                  patient.before.chemo="patient")})))
  meta$TherapyStatus = as.factor(unlist(apply(meta,
                                              1,
                                              function(row) {switch(row["Sample.type.1"],
                                                                    sibling ="before.chemo",
                                                                    patient.after.chemo="after.chemo",
                                                                    patient.before.chemo="before.chemo")})))
  
  meta$Antibiotic = tolower(meta$Antibiotic.treatment.within.the.last.1.months.) != "no"
  meta$Antibiotic[str_blank(meta$Antibiotic.treatment.within.the.last.1.months.)] = NA
  meta$Antibiotic[tolower(meta$Antibiotic.treatment.within.the.last.1.months.)=="unknown"] = NA
  ## R gotcha: if logical converted to factor, it loses is.logical() status, and gives
  ## unexpected results in various statements that assume logical inputs
  #meta$Antibiotic = as.factor(meta$Antibiotic)
  
  meta$Fever = tolower(meta$Fever.descr) != "no"
  meta$Fever[str_blank(meta$Fever.descr)] = NA
  meta$Fever[tolower(meta$Fever.descr)=="unknown"] = NA
  ## R gotcha: if logical converted to factor, it loses is.logical() status, and gives
  ## unexpected results in various statements that assume logical inputs
  #meta$Fever = as.factor(meta$Fever)
  
  ## ggplot needs Date object
  Specimen.Collection.Date = 
    as.Date(as.character(meta$Specimen.Collection.Date),format = "%m/%d/%Y")
  
  ## arithmetics work better with date object
  Specimen.Collection.Date.d = 
    as.date(as.character(meta$Specimen.Collection.Date),order="mdy")
  meta$age = as.numeric(
    Specimen.Collection.Date.d
    - 
      as.date(paste(meta$YearOfBirth,"-06-01",sep=""),order="ymd")
  )/365
  
  meta$age.quant = quantcut(meta$age)
  
  meta$Specimen.Collection.Date = Specimen.Collection.Date
  
  meta$visit = as.numeric(laply(as.character(meta$SampleID.compound),
                                function(x) unlist(strsplit(x,".v",fixed=T))[2]))
  meta$SampleID.1 = as.factor(paste(meta$SubjectID,meta$Sample.type.1,sep="."))
  meta$Sample.type = as.factor(meta$Sample.type)
  meta$Sample.type.1 = as.factor(meta$Sample.type.1)
  meta$TherapyStatus = as.factor(meta$TherapyStatus)
  meta$SampleID = as.factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  meta$Sample.type.1 = relevel(meta$Sample.type.1,"sibling")
  meta$Sample.type = relevel(meta$Sample.type,"sibling")
  meta$TherapyStatus = relevel(meta$TherapyStatus,"before.chemo")
  
  meta.aggr = join(meta,
                   ddply(meta,"SubjectID",summarise,
                         Antibiotic.Before.Therapy=any(ifelse(visit==1,as.logical(Antibiotic),F))),
                   by="SubjectID",
                   match="first")
  stopifnot(nrow(meta.aggr)==nrow(meta))

  meta = meta.aggr

  ## ignore antibiotic status in siblings - this field is for faceted plots
  meta$Sample.type.Antibio.Before = factor(with(meta,
                                     ifelse(Sample.type=="patient",
                                            paste(Sample.type,Antibiotic.Before.Therapy,"."),
                                            as.character(Sample.type))
  ))
  
  make.global(meta)
  return (meta)
}


## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.choc <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~Sample.type+TherapyStatus","~Antibiotic.Before.Therapy + Sample.type",
                        "~Fever + Sample.type",
                        "~Sample.type+visit","~FamilyID","~Sample.type.1","~SubjectID")
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

gen.tasks.choc <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    #DEBUG:
    taxa.levels = c(2,3,4,5,6,"otu")
    #taxa.levels = c(2)
    
    descr = "All samples, no aggregation"
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = "2b256d57d3fd3114dc1c0391cc87c2f8.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
      otu.shared.file="b98469d29ce9694f9f0254f415aab6f6.files_x1.sorted.0.03.shared"
      cons.taxonomy.file="71f294ce8d2a4b33bf6656ef72b0d5f6.files_x1.sorted.0.03.cons.taxonomy"
      meta.file="CHOC_ALL_Samples_Metadata_Nov-19-2014.txt"
      load.meta.method=load.meta.choc
      load.meta.options=list()    
    })
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      meta = m_a$attr
      meta.aggr = join(meta,
                            ddply(meta,"SubjectID",summarise,
                                  visit.max=max(visit),
                                  visit.min=min(visit),
                                  visit.1=any(visit==1),
                                  visit.2=any(visit==2)),
                            by="SubjectID",
                            match="first")
      stopifnot(!any(is.na(meta.aggr$visit.max)) && 
                  nrow(meta.aggr)==nrow(meta))
      
      meta = meta.aggr
      meta.aggr = join(meta,
                       ddply(meta,"FamilyID",summarise,
                             has.sibling=("patient" %in% Sample.type) && ("sibling" %in% Sample.type)),
                       by="FamilyID",
                       match="first")
      stopifnot(!any(is.na(meta.aggr$has.sibling)) && 
                  nrow(meta.aggr)==nrow(meta))
            
      m_a$attr = meta.aggr
      
      ##As of 2014-11-05, there are only 6 samples at visit 5, and less in higher visits
      ##and their profiles look similar to visit 4
      m_a = subset.m_a(m_a,subset=(m_a$attr$visit<=4))
      return(m_a)
    }
    
    summary.meta.method=summary.meta.choc
    
    test.counts.task = within(test.counts.task, {
      
      norm.count.task = within(norm.count.task, {
        method="norm.ihs.prop"
      })
      
    })
    
  })
  
  
  task1 = within( task0, {

    descr = "All samples, no aggregation"    
    
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
                            c("Sample.type.Antibio.Before","visit"),
                            c("Antibiotic.Before.Therapy","Sample.type.1"))
        clade.meta.x.vars=c("visit")
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c("Sample.type","visit","Antibiotic.Before.Therapy")
      })
      
    })
    
  })

  task1.1 = within( task1, {
    
    descr = paste(descr,"Additional tests")
    
    taxa.levels = c(2)
    
    do.summary.meta = F
    
    do.tests = T
    
    
    test.counts.task = within(test.counts.task, {
      
      do.divrich = c()
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.extra.method = taxa.levels
      
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task) {
          require(vegan)
          
          report$add.header('Testing that patients move closer to sibling profiles over time')
          
          if(!is.null(norm.count.task)) {
            m_a <- norm.count.report(m_a,
                                     descr="Profile time trend",
                                     norm.count.task=norm.count.task)
          }
          
          make.global(m_a)
          
        }
        norm.count.task = within(norm.count.task, {
          method="norm.prop"
          #drop.features = list()
        })
        
      })
      
      
    })
    
  })
  
  
  task2 = within( task0, {
    
    main.meta.var = "Sample.type"
    
    descr = "Patient/Sibling samples before therapy aggregated by SubjectID"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TherapyStatus=="before.chemo"))
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
        do.clade.meta=F
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
      m_a = subset.m_a(m_a,subset=(m_a$attr$has.sibling)) 
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
          block.attr = "FamilyID"
          type="paired"
          #replicates=0
        })
      })

      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata="FamilyID",
               descr="Association with the patient/control status paired by family")
        )
        
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="Sample.type",
                                        block.attr="FamilyID",
                                        n.perm=4000,
                                        norm.count.task=norm.count.task.extra
          )
        }
        norm.count.task.extra = within(norm.count.task, {
          method="norm.prop"
          #drop.features = list()
        })
        
      })
      
      
    })
    
  })
  
  
  task3 = within( task0, {
    
    main.meta.var = "TherapyStatus"
    
    descr = "Patients' samples at visits 1 (before therapy) and 2 (after therapy), only paired samples"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient" 
                                   & m_a$attr$visit <= 2 
                                   & m_a$attr$visit.1
                                   & m_a$attr$visit.2))
      #DEBUG: scrambling SubjectID of before.chemo to see how paired tests behave on random pairings
      #m_a$attr$SubjectID[m_a$attr$TherapyStatus=="before.chemo"] = 
      #  sample(m_a$attr$SubjectID[m_a$attr$TherapyStatus=="before.chemo"])
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
               descr="Association with therapy status unpaired by subject"),
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with therapy status paired by subject"),
          list(formula.rhs=paste("Antibiotic.Before.Therapy * ", main.meta.var),
               strata=NULL,
               descr="Association with Antibiotic use before therapy and therapy status")
        )
        
        #dist.metr="euclidian"
        #col.trans="standardize"

        #norm.count.task = within(norm.count.task, {
        #  method="norm.clr"
        #  drop.features = list()
        #})
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"Antibiotic.Before.Therapy"))
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic.Before.Therapy")
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="TherapyStatus",
                                        block.attr="SubjectID",
                                        n.perm=8000,
                                        #dist.metr="euclidian",
                                        col.trans="ident",
                                        norm.count.task=norm.count.task.extra
                                        )
        }
        norm.count.task.extra = within(norm.count.task, {
          method="norm.prop"
          #drop.features = list()
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
        formula.rhs = sprintf("Antibiotic.Before.Therapy+%s",main.meta.var)
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
        id.vars.list = list(c(main.meta.var,"Antibiotic.Before.Therapy"))
        do.profile=T
        do.clade.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic.Before.Therapy")
      })
      
    })
    
  })
  
  #return (list(task1.1))
  return (list(task1,task2,task2.1,task3,task3.1,task4))
}


'
> library(coin)
> m_a$attr$SubjectID = factor(m_a$attr$SubjectID)
> wilcoxsign_test(Roseburia.Otu2637~TherapyStatus|SubjectID,cbind(m_a$count,m_a$attr))

Asymptotic Wilcoxon-Signed-Rank Test (zeros handled a la Pratt)

data:  y by x (neg, pos) 
stratified by block
Z = 2.3863, p-value = 0.01702
alternative hypothesis: true mu is not equal to 0

Warning message:
In wilcoxsign_test.IndependenceProblem(object = <S4 object of class "IndependenceProblem">) :
Handling of zeros defaults to \'Pratt\' in newer versions of coin
#GeneSelector::WilcoxonRanking paired gave Roseburia.Otu2637   stat=167 	p.value=0.017181

> wilcox_test(Roseburia.Otu2637~TherapyStatus,cbind(m_a$count,m_a$attr))

Asymptotic Wilcoxon Mann-Whitney Rank Sum Test

data:  Roseburia.Otu2637 by
TherapyStatus (before.chemo, after.chemo)
Z = -2.7655, p-value = 0.005683
alternative hypothesis: true mu is not equal to 0
#GeneSelector::WilcoxonRanking unpaired gave Roseburia.Otu2637   stat=295.0 	p.value=0.008712
'

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
                       title="Analysis of CHOC ALL 16S data",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.choc
)

report$save()

