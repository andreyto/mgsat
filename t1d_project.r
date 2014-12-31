
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)

load.meta.t1d <- function(file.name,batch=NULL,aggr.var=NULL) {
  
  meta =read.delim(file.name, header=TRUE,stringsAsFactors=T, sep="\t")
  make.global(meta)
  names(meta) =
    replace.col.names(
      names(meta),
      c("YAP_Aliq_ID","Sample.Name","family.ID", "Subject.ID", 
        "gender", "autoantibody.status", "T1D.status",
        "year.of.birth","date.collected","date.of.diagnosis",
        "timestamp","sample.type"),
      c("AliquotID","SampleID",      "FamilyID",            "SubjectID",            
        "gender",           "AA",                  "T1D",
        "YearOfBirth","Specimen.Collection.Date","Date.of.Diagnosis",
        "Timestamp","Specimen.Type")
    )
  
  meta = meta[!duplicated(meta$SampleID),]
  
  meta$gender = unlist( lapply( meta$gender, substr, 1,1  ))
  
  meta$Batch = 4
  ## We set Batch=0 downstream as a special value, it should not be
  ## present already
  stopifnot(sum(meta$Batch==0) == 0)
  meta$BatchRepeat=as.factor(1)
  
  meta$Specimen.Collection.Date = 
    as.date(as.character(meta$Specimen.Collection.Date),order="mdy")
  meta$age = as.numeric(
    meta$Specimen.Collection.Date
    - 
      as.date(paste(meta$YearOfBirth,"-06-01",sep=""),order="ymd")
  )/365
  
  ##day is fake
  meta$Date.of.Diagnosis = as.date(as.character(meta$Date.of.Diagnosis),order="mdy")
  
  meta$YearsSinceDiagnosis = as.numeric(
    meta$Specimen.Collection.Date
    -
      meta$Date.of.Diagnosis
  )/365
  
  meta$YearsSinceDiagnosis.quant = quantcut(meta$YearsSinceDiagnosis)
  
  meta$T1D[meta$T1D=="Unknown"] = "Control"
  
  meta$Timestamp = 
    strptime(as.character(meta$Timestamp),format = "%m/%d/%Y %H:%M")
  
  meta$BioSampleID = paste(meta$SubjectID,meta$Specimen.Type,as.numeric(meta$Timestamp),sep="_")
  
  make.global(meta)
  
  meta$age.quant = quantcut(meta$age)
  meta$A1C = as.double(as.character(meta$A1C))
  #meta$A1C[is.na(meta$A1C)] = 6.2
  meta$A1C.quant = quantcut(as.double(as.character(meta$A1C)))
  
  meta = arrange(meta,SubjectID,Timestamp)
  
  ## Filtering and aggregation of rows has to be done before we start marking repeats for throwing
  ## out repeats
  if(!is.null(batch)) {
    report$add.p(paste(c("Filtering metadata by batch:",batch),collapse=" "))
    batch.mask = (meta$Batch %in% batch)
    meta = meta[batch.mask,]
  }  
  
  ## That would leave only AA defined samples or T1D patients (under
  ## an assumption that T1D patients are always AA positive even if AA is not measured)
  #meta = meta[meta$AA!="Unknown" | meta$T1D == "T1D",]
  
  ## Subject is at least first repeat, but bio sample is first occurence, in time order 
  ## (e.g. second visit)
  isBioRepeatFirst = duplicated(meta$SubjectID) & ! duplicated(meta$BioSampleID)
  
  BioSampleIDRep = meta$BioSampleID[isBioRepeatFirst]
  
  meta$isBioRepeat = meta$BioSampleID %in% BioSampleIDRep
  
  ## Now do the same as above, but ignore failed samples when looking
  ## for bio repeats
  meta_qa = meta[meta$Sample.QA=="PASS",]
  
  isBioRepeatFirst = duplicated(meta_qa$SubjectID) & ! duplicated(meta_qa$BioSampleID)
  
  BioSampleIDRep = meta_qa$BioSampleID[isBioRepeatFirst]
  
  ## mark all failed samples or bio repeats following at least one good sample
  meta$isBioRepeatOrFailed = (meta$BioSampleID %in% BioSampleIDRep) | 
    ! (meta$Sample.QA=="PASS")
  
  SampleIDRep = meta_qa$SampleID[duplicated(meta_qa$SubjectID)]
  
  ## mark all failed samples or repeats following at least one good sample
  meta$isRepeatOrFailed = (meta$SampleID %in% SampleIDRep) | 
    ! (meta$Sample.QA=="PASS")
  
  ##
  ## This section builds a difference in months between next samples and first sample for every subject
  ##
  meta_keys = meta[,c("SubjectID","SampleID","BioSampleID","Timestamp","isBioRepeat","isBioRepeatOrFailed")]  
  meta_keys_non_rep = meta_keys[!meta$isBioRepeatOrFail,]
  #join() will keep duplicate names, so we rename the names, except SubjectID merge key
  names(meta_keys_non_rep) = paste(names(meta_keys_non_rep),".first",sep="")
  names(meta_keys_non_rep)[1] = "SubjectID"
  meta_keys_1 = join(meta_keys,meta_keys_non_rep,by=c("SubjectID"),match="first",type="left")
  stopifnot(meta$SampleID==meta_keys_1$SampleID)
  stopifnot(is.na(meta_keys_1$Timestamp.first) |
              (meta_keys_1$Timestamp >= meta_keys_1$Timestamp.first) |
              ! (meta$Sample.QA=="PASS")) #should not see otherwise
  meta_keys_1$MonthsAfterFirstBioSample = as.numeric(meta_keys_1$Timestamp - meta_keys_1$Timestamp.first,units="days")/30
  ## NA will be at records that did not pass QA. Set to 0.
  meta_keys_1$MonthsAfterFirstBioSample[is.na(meta_keys_1$MonthsAfterFirstBioSample)] = 0
  ##
  ## uncheck those bio repeats that passed QA and were made less than 9 mon after the first sampling - they
  ## are testing intentional repeats and should be aggregated just like batch repeats when analyzing 
  ## unpaired group differences
  ##
  meta_keys_1$isAnnualRepeatOrFailed = meta_keys_1$isBioRepeatOrFailed & 
    ! (
      (meta$Sample.QA=="PASS") & 
        meta_keys_1$MonthsAfterFirstBioSample < 9
    )
  make.global(meta_keys_1)
  meta$isAnnualRepeatOrFailed = meta_keys_1$isAnnualRepeatOrFailed
  meta$MonthsAfterFirstBioSample = meta_keys_1$MonthsAfterFirstBioSample
  meta$TimestampMonth = as.numeric(meta$Timestamp - min(meta$Timestamp),units="days")/30  
  meta$Timestamp = as.numeric(meta$Timestamp)
  meta$TimestampDate = as.Date(timeDate(meta$Timestamp))
  
  row.names(meta) = meta$SampleID
  
  ##pass filter argument and filter here before aggregation but after repeats etc are marked
  
  if(!is.null(aggr.var)) {
    
    meta = aggregate.by.meta.data(meta,
                                  group_col=aggr.var,
                                  col_ignore=names(meta))  
  }
  
  ##Setting Control to base level plays nice with DESeq2 defaults
  meta$T1D = relevel(meta$T1D,"Control")
  return (meta)
}


## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.t1d <- function(m_a) {
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr
  
  report$add.printed(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~T1D","~T1D+FamilyID","~FamilyID","~SubjectID","~AA+T1D")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(meta,{
    
    report$add(aov(age~T1D),
                       caption="ANOVA for age and cohort")
    report$add(qplot(T1D,age,geom="violin"),
               caption="Violin plot for age and cohort")
    
    #summary(glht(lm(age~visit,data=meta[meta$Sample.type=="patient",]),linfct="visit=0"))
    #summary(glht(lmer(age~visit+(visit|Sample.type),data=meta),linfct="visit=0"))
    report$add(ggplot(meta,aes(x=age,y=YearsSinceDiagnosis,color=T1D))+
                 geom_point()+
                 stat_smooth(method="loess", se = T,degree=1,size=1),
               caption="Plot for age and time since diagnosis with Loess trend line")
    
  })
  
  with(meta,{
    report$add(cor.test(YearsSinceDiagnosis,
                                Timestamp,
                                method="spearman"),
                       caption="Spearman RHO for metadata variables")
    report$add(cor.test(YearsSinceDiagnosis,
                                age,
                                method="spearman"),
                       caption="Spearman RHO for metadata variables")
    report$add(cor.test(YearsSinceDiagnosis,
                                A1C,
                                method="spearman"),
                       caption="Spearman RHO for metadata variables")
    
  })
  
  
}

## This function must generate a lits with analysis tasks

gen.tasks.t1d <- function() {
  
  get.taxa.meta.aggr.base<-function(m_a) { 
    m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.QA=="PASS"))

    report$add.p(paste("After filtering for QAed samples:",nrow(m_a$count)))
    
    m_a = subset.m_a(m_a,subset=(!m_a$attr$isRepeatOrFailed))

    report$add.p(paste("After filtering for repeated samples per subject:",nrow(m_a$count)))
    
    return (m_a)
  }

  task0 = within( mgsat.16s.task.template, {
    
  taxa.levels = c(2,6,"otu")
  
  descr = "All samples, aggregated by AliquotID"
  
  main.meta.var = "T1D"
  
  read.data.task.yap = within(read.data.task, {
    taxa.summary.file = "a93eeedeef2878be17d30f27b1b0de1c.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
    otu.shared.file="d9f44114ac6369cc66978002228bc5d3.files_x1.sorted.0.03.shared"
    cons.taxonomy.file="4272870500ad07a7270f9772581fe7fe.files_x1.sorted.0.03.cons.taxonomy"
  })

  read.data.task.mothur.cdhit = within(read.data.task, {
    taxa.summary.file = "stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
  })
  
  read.data.task = within(read.data.task.yap, {
    meta.file="aliq_id_to_metadata_for_T1D_YAP_run_20140922_batches.tsv"
    load.meta.method=load.meta.t1d
    load.meta.options=list(aggr.var="AliquotID")    
  })
  
  get.taxa.meta.aggr<-function(taxa.meta) { 
    return (get.taxa.meta.aggr.base(taxa.meta))
  }
  
  summary.meta.method=summary.meta.t1d
  
  })
  
  task1 = within( task0, {
    
  do.summary.meta = T
  
  do.tests = T

  summary.meta.task = within(summary.meta.task, {
    meta.x.vars = c("Timestamp")
    group.var = c(main.meta.var)
  })

  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    do.deseq2 = T
    do.adonis = T
    do.genesel = T
    do.stabsel = T
    do.glmer = F
    #do.divrich = c()

    do.plot.profiles.abund=T
    do.heatmap.abund=T
    
    divrich.task = within(divrich.task,{
      #n.rar.rep=4
      is.raw.count.data=T
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })      
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
      
      #dist.metr="euclidean"
      #col.trans="standardize"
      data.descr="normalized counts"
      
      tasks = list(
        list(formula.rhs="T1D",
             strata=NULL,
             descr="Association with the patient/control status unpaired"),
        list(formula.rhs="T1D",
             strata="FamilyID",
             descr="Association with the patient/control status paired by family"),
        list(formula.rhs="T1D*age",
             strata=NULL,
             descr="Association with the patient/control status and age unpaired")
        
      )
      
    })
    
    glmer.task = within(glmer.task, {
      
      tasks = list(list(
        descr.extra = "",
        formula.rhs = paste(main.meta.var,"(1|SampleID)",sep="+"),
        linfct=c("T1DT1D = 0")
      ))
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
      clade.meta.x.vars=c("YearsSinceDiagnosis","TimestampDate","age")
      do.profile=T
      do.clade.meta=F
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var)
    })
    
  })
  
})


task2 = within( task0, {
  
  do.summary.meta = F
  
  do.tests = T

  descr = "Patient samples only, association with years since diagnosis"
  
  get.taxa.meta.aggr<-function(m_a) { 
    m_a = get.taxa.meta.aggr.base(m_a)
    m_a = subset.m_a(m_a,subset=(m_a$attr$T1D=="T1D"))
    return (m_a)
  }
  
  main.meta.var = "YearsSinceDiagnosis.quant"
  main.meta.var.cont = "YearsSinceDiagnosis" 
  
  
  summary.meta.task = within(summary.meta.task, {
    meta.x.vars = c("Timestamp")
    group.var = c(main.meta.var)
  })
  
  
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    do.deseq2 = T
    do.adonis = T
    do.genesel = F
    do.stabsel = T
    do.glmer = F
    #do.divrich = c()
    
    do.plot.profiles.abund=T
    do.heatmap.abund=T
    
    divrich.task = within(divrich.task,{
      #n.rar.rep=4
      is.raw.count.data=T
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
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
      )
    })
        
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
      clade.meta.x.vars=c("age")
      do.profile=T
      do.clade.meta=F
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var.cont)
    })
    
  })
  
})

task3 = within( task2, {
  
  do.summary.meta = F
  
  do.tests = T
  
  descr = "Patient samples only, association with A1C"
  
  main.meta.var = "A1C.quant"
  main.meta.var.cont = "A1C" 
  
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    divrich.task = within(divrich.task,{
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
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
    })
    
    adonis.task = within(adonis.task, {
      
      tasks = list(
        list(formula.rhs=main.meta.var.cont,
             strata=NULL,
             descr=sprintf("Association with %s",main.meta.var.cont))        
      )
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var.cont)
    })
    
  })
  
})

task3.1 = within( task3, {
  
  do.summary.meta = F
  
  do.tests = T
  
  descr = "Patient samples only, association with age and A1C"
  
  main.meta.var = "A1C.quant"
  main.meta.var.cont = "A1C"
  control.meta.var = "age.quant"
  control.meta.var.cont = "age"
  
  test.counts.task = within(test.counts.task, {

    do.deseq2 = T
    do.adonis = F
    do.genesel = F
    do.stabsel = F
    do.divrich = c()
    
    do.plot.profiles.abund=F
    do.heatmap.abund=F
    
    divrich.task = within(divrich.task,{
      group.attr = NULL
      counts.glm.task = within(list(),{
        formula.rhs = paste(control.meta.var.cont,main.meta.var.cont,sep="+")
      })      
    })
    
    deseq2.task = within(deseq2.task, {
      formula.rhs = paste(control.meta.var.cont,main.meta.var.cont,sep="+")
    })
    
    adonis.task = within(adonis.task, {
      
      tasks = list(
        list(formula.rhs=paste(control.meta.var.cont,main.meta.var.cont,sep="+"),
             strata=NULL,
             descr=sprintf("Association with %s and %s",control.meta.var.cont,main.meta.var.cont))
      )
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(control.meta.var,main.meta.var))
    })
    
  })
  
})

task4 = within( task1, {
  
  do.summary.meta = F
  
  do.tests = T
  
  descr = "Control samples only, association with Auto Antibody status (AA)"
  
  get.taxa.meta.aggr<-function(m_a) { 
    m_a = get.taxa.meta.aggr.base(m_a)
    m_a = subset.m_a(m_a,subset=(m_a$attr$T1D!="T1D" & m_a$attr$AA %in% c("AA+","AA-")))
    return (m_a)
  }
    
  main.meta.var = "AA"
  main.meta.var.cont = "age" 
  
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    divrich.task = within(divrich.task,{
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })      
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
             descr=sprintf("Association with %s",main.meta.var.cont))        
      )
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var)
    })
    
  })
  
})

task4.1 = within( task4, {
  
  do.summary.meta = F
  
  do.tests = T
  
  descr = "All subjects where AA is defined, association with Auto Antibody status (AA)"
  
  get.taxa.meta.aggr<-function(m_a) { 
    m_a = get.taxa.meta.aggr.base(m_a)
    m_a = subset.m_a(m_a,subset=(m_a$attr$AA %in% c("AA+","AA-")))
    m_a$attr$AA = factor(m_a$attr$AA)
    return (m_a)
  }
  
  test.counts.task = within(test.counts.task, {
    
    do.deseq2 = T
    do.adonis = F
    do.genesel = F
    do.stabsel = F
    do.divrich = c()
    do.plot.profiles.abund=F
    do.heatmap.abund=F    
  })
})

task1.1 = within( task1, {
  
  do.summary.meta = F
  
  do.tests = T
  
  descr = "All samples, Adonis test for association with T1D"
  
  main.meta.var.cont = "age" 
  
  test.counts.task = within(test.counts.task, {
    
    do.deseq2 = F
    do.adonis = T
    do.genesel = F
    do.stabsel = F
    do.divrich = c()
    do.plot.profiles.abund=F
    do.heatmap.abund=F    
    
    deseq2.task = within(deseq2.task, {
      formula.rhs = sprintf("age.quant+%s",main.meta.var)
    })
        
  })
  
})


task5.1 = within( task1, {
  descr = "Control samples with Auto Antibody status (AA) AA- vs T1D samples with A1C <= 7.3"
  
  get.taxa.meta.aggr<-function(m_a) { 
    m_a = get.taxa.meta.aggr.base(m_a)
    m_a = subset.m_a(m_a,subset=((m_a$attr$T1D!="T1D" & m_a$attr$AA == "AA-") | 
                                   (m_a$attr$T1D=="T1D" & m_a$attr$A1C <= 7.3)))
    return (m_a)
  }  
})

task5.2 = within( task1, {
  descr = "Control samples with Auto Antibody status (AA) AA- vs T1D samples with AA+"
  
  get.taxa.meta.aggr<-function(m_a) { 
    m_a = get.taxa.meta.aggr.base(m_a)
    m_a = subset.m_a(m_a,subset=((m_a$attr$T1D!="T1D" & m_a$attr$AA == "AA-") | 
                                   (m_a$attr$T1D=="T1D" & m_a$attr$AA == "AA+")))
    return (m_a)
  }  
})

extra.tasks = foreach(task=list(task5.1,task5.2)) %do% {
  within( task, {
  
  do.summary.meta = F
  
  do.tests = T
  
  taxa.levels = c(2,6)
  
  main.meta.var = "T1D"
  main.meta.var.cont = "A1C" 
  
  test.counts.task = within(test.counts.task, {
    
    do.deseq2 = T
    do.adonis = F
    do.genesel = F
    do.stabsel = T
    do.divrich = c()
    do.plot.profiles.abund=F
    do.heatmap.abund=F    
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    divrich.task = within(divrich.task,{
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })      
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
             descr=sprintf("Association with %s",main.meta.var))        
      )
    })
    
    plot.profiles.task = within(plot.profiles.task, {
      id.vars.list = list(c(main.meta.var))
    })
    
    heatmap.abund.task = within(heatmap.abund.task,{
      attr.annot.names=c(main.meta.var)
    })
    
  })
  
})
}

return (list(task4.1))
return (list(task1,task2,task3,task3.1,task4))
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
                       title="Analysis of T1D 16S data",
                       out.file.md="report.md",
                       incremental.save=T)


res = proc.project(
  task.generator.method=gen.tasks.t1d
)

report$save()

