
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
  meta$Antibiotic = as.factor(meta$Antibiotic)
  
  meta$Fever = tolower(meta$Fever.descr) != "no"
  meta$Fever[str_blank(meta$Fever.descr)] = NA
  meta$Fever[tolower(meta$Fever.descr)=="unknown"] = NA
  meta$Fever = as.factor(meta$Fever)
  
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
  
  meta$Specimen.Collection.Date = Specimen.Collection.Date
  
  meta$visit = as.numeric(laply(as.character(meta$SampleID.compound),
                                function(x) unlist(strsplit(x,".v",fixed=T))[2]))
  meta.visit.max = join(meta,
                        ddply(meta,"SubjectID",summarise,visit.max=max(visit)),
                        by="SubjectID",
                        match="first")
  stopifnot(!any(is.na(meta.visit.max$visit.max)) && 
              nrow(meta.visit.max)==nrow(meta))
  meta = meta.visit.max
  meta$SampleID.1 = as.factor(paste(meta$SubjectID,meta$Sample.type.1,sep="."))
  meta$Sample.type = as.factor(meta$Sample.type)
  meta$Sample.type.1 = as.factor(meta$Sample.type.1)
  meta$TherapyStatus = as.factor(meta$TherapyStatus)
  meta$SampleID = as.factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  meta$Sample.type.1 = relevel(meta$Sample.type.1,"sibling")
  meta$Sample.type = relevel(meta$Sample.type,"sibling")
  meta$TherapyStatus = relevel(meta$TherapyStatus,"before.chemo")
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
  
  xtabs.formulas = list("~Sample.type+TherapyStatus","~Antibiotic + Sample.type",
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
    
    taxa.levels = c(2,6,"otu")
    
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
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      return(m_a)
    }    
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = c("visit")
      group.var = c("Sample.type","visit")
    })
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c(2)
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
        id.vars.list = list(c("Sample.type","visit"))
        clade.meta.x.vars=c("visit")
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c("Sample.type","visit","Antibiotic")
      })
      
    })
    
  })
  
  
  task2 = within( task0, {
    
    main.meta.var = "Sample.type"
    
    descr = "Samples before therapy aggregated by SubjectID"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TherapyStatus=="before.chemo"))
      m_a = aggregate.by.meta.data.m_a(m_a,group_col="SubjectID")
      return(m_a)
    }
    
    test.counts.task = within(test.counts.task, {
      
      #do.divrich = c()
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = T
      do.glmer = F
      do.plot.profiles.abund=F
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
        resp.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata=NULL,
               descr="Association with the patient/control status unpaired"),
          list(formula.rhs=main.meta.var,
               strata="FamilyID",
               descr="Association with the patient/control status paired by family"),
          list(formula.rhs=paste(main.meta.var,"* age"),
               strata="FamilyID",
               descr="Association with the patient/control status and age  paired by family")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"Antibiotic"))
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic")
      })
      
    })
    
  })

  task3 = within( task0, {
    
    main.meta.var = "TherapyStatus"
    
    descr = "Patients samples at visits 1 (before therapy) and 2 (after therapy)"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient" & m_a$attr$visit <= 2))
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
        resp.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs=main.meta.var,
               strata="SubjectID",
               descr="Association with therapy status paired by subject"),
          list(formula.rhs=c("Antibiotic * ", main.meta.var),
               strata="SubjectID",
               descr="Association with Antibiotic use and therapy status paired by subject")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var),c(main.meta.var,"Antibiotic"))
        do.profile=T
        do.clade.meta=F
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic")
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
        resp.attr = main.meta.var
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
        id.vars.list = list(c(main.meta.var))
        do.profile=T
        do.clade.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic")
      })
      
    })
    
  })
  
  
  return (list(task3))
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
                       title="Analysis of CHOC ALL 16S data",
                       out.file.md="report.md",
                       incremental.save=T)


res = proc.project(
  task.generator.method=gen.tasks.choc
)

report$save()

