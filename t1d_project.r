
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
    
    meta = aggregate_by_meta_data(meta,
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
  
  xtabs.formulas = list("~T1D","~T1D+FamilyID","~FamilyID","~SubjectID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add.table(fact.xtabs,show.row.names=T,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(meta,{
    
    report$add.printed(summary(aov(age~T1D)),
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
    report$add.printed(cor.test(YearsSinceDiagnosis,
                                Timestamp,
                                method="spearman"),
                       caption="Spearman RHO for metadata variables")
    report$add.printed(cor.test(YearsSinceDiagnosis,
                                age,
                                method="spearman"),
                       caption="Spearman RHO for metadata variables")
    report$add.printed(cor.test(YearsSinceDiagnosis,
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
  
  descr = "All samples, aggregated by AliquotID"
  
  main.meta.var = "T1D"
  
  read.data.task = within(read.data.task, {
    taxa.summary.file = "a93eeedeef2878be17d30f27b1b0de1c.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
    otu.shared.file="d9f44114ac6369cc66978002228bc5d3.files_x1.sorted.0.03.shared"
    cons.taxonomy.file="4272870500ad07a7270f9772581fe7fe.files_x1.sorted.0.03.cons.taxonomy"
    meta.file="aliq_id_to_metadata_for_T1D_YAP_run_20140922.tsv"
    load.meta.method=load.meta.t1d
    load.meta.options=list(aggr.var="AliquotID")    
  })
  
  get.taxa.meta.aggr<-function(taxa.meta) { 
    return (get.taxa.meta.aggr.base(taxa.meta))
  }
  
  summary.meta.method=summary.meta.t1d
  
  })
  
  task1 = within( task0, {
  
  taxa.levels = c(6)
    
  do.summary.meta = T
  
  do.tests = T

  summary.meta.task = within(summary.meta.task, {
    meta.x.vars = c("Timestamp")
    group.var = c(main.meta.var)
  })
  
  #DEBUG:
  #get.taxa.meta.aggr<-function(taxa.meta) { 
  #  d = get.taxa.meta.aggr.base(taxa.meta)
  #  d$data = d$data[sample.int(nrow(d$data),40),]
  #  d
  #}

  #get.taxa.meta.aggr<-function(taxa.meta) { 
  #  taxa.meta.aggr = get.taxa.meta.aggr.base(taxa.meta)
  #  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,T1D=="T1D")  
  #  return (taxa.meta.aggr)
  #}
    
  test.counts.task = within(test.counts.task, {
    
    norm.count.task = within(norm.count.task, {
      method="norm.ihs.prop"
    })
    
    do.deseq2 = F
    do.adonis = F
    do.genesel = T
    do.stabsel = F
    do.glmer = F
    #do.divrich = c()

    do.plot.profiles.abund=T
    do.heatmap.abund=T
    
    divrich.task = within(divrich.task,{
      #n.rar.rep=4
      is.raw.count.data=T
      pool.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })      
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

return (list(task1))
}

hide_debug <- function() {


task2 = within(
  task1,
{
  
  descr = "Patient samples only"
  
  do.summary.meta = F
  do.genesel=F
  
  get.taxa.meta.aggr<-function(taxa.meta) { 
    taxa.meta.aggr = get.taxa.meta.aggr.base(taxa.meta)
    taxa.meta.aggr$data = subset(taxa.meta.aggr$data,T1D=="T1D")  
    return (taxa.meta.aggr)
  }
  
  stability.resp.attr = "YearsSinceDiagnosis" 
  stability.model.family = "gaussian"
  
  genesel.resp.attr = stability.resp.attr
  
  adonis.tasks = list(
    list(formula_rhs="YearsSinceDiagnosis",
         strata=NULL,
         descr="Association with the years since diagnosis"),
    list(formula_rhs="age*YearsSinceDiagnosis",
         strata=NULL,
         descr="Association with the age and years since diagnosis"),
    list(formula_rhs="A1C",
         strata=NULL,
         descr="Association with A1C")    
  )
  
  glmer.task = list(
    formula.rhs = "YearsSinceDiagnosis + (1|FamilyID/SubjectID)",
    linfct=c("YearsSinceDiagnosis = 0"),
    descr=sprintf(glmer.descr.tpl,"YearsSinceDiagnosis")
  )  
  
  plot.group = list(
    c("YearsSinceDiagnosis.quant")
  )            
  
  clade.meta.x.vars=c("YearsSinceDiagnosis")
  
  heatmap.task = list(
    attr.annot.names=c("YearsSinceDiagnosis"),
    attr.row.labels=NULL
  )
  
}
)



return (list(task1,task2))
#return (list(task5))
#return (list(task2))
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
                       title="16S Analysis of T1D cohort",
                       out.file.md="report.md",
                       incremental.save=T)


res = proc.project(
  task.generator.method=gen.tasks.t1d
)

report$save()

