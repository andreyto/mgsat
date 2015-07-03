
compute.has.sibling <- function(meta) {
  fam.with.sibl = ddply(meta,"FamilyID",summarise,
                        has.sibling=("patient" %in% Sample.type) & ("sibling" %in% Sample.type))
  
  
  fam.with.sibl = fam.with.sibl[fam.with.sibl$has.sibling,,drop=F]
  
  return (meta$FamilyID %in% fam.with.sibl$FamilyID)
}

as.Date.date.time <- function(date,time,date.format="%Y-%m-%d",default.time="12:00") {
  stopifnot(length(date)==length(time))
  time = toupper(time)
  time[grepl(".*UNK.*",time)] = NA
  time[grepl(".*NR.*",time)] = NA
  time[grepl(".*NA.*",time)] = NA
  isel = nchar(time)==4 & !grepl(".*.M$",time) & !grepl(":",time)
  time[isel] = paste(substr(time[isel],1,2),substring(time[isel],3),sep=":")  
  isel = nchar(time)==3 & !grepl(".*.M$",time) & !grepl(":",time)
  time[isel] = paste(substr(time[isel],1,1),substring(time[isel],2),sep=":")
  isel = nchar(time)==3 & grepl(".*.M$",time) & !grepl(":",time)
  time[isel] = paste(substr(time[isel],1,1),":00",substring(time[isel],2),sep="")
  isel = !is.na(time) & nchar(time)==2 & !grepl(".*.M$",time) & !grepl(":",time)
  time[isel] = paste(substr(time[isel],1,1),":",substring(time[isel],2),sep="")
  time[is.na(time)] = default.time
  #dt = time
  dt = paste(date,time,sep=" ")
  dt.24 = as.POSIXct(dt,format = paste(date.format,"%H:%M"))
  dt.pm = as.POSIXct(dt,format = paste(date.format,"%I:%M%p"))
  isel = !grepl(".*.M$",time)
  dt.pm[isel] = dt.24[isel]
  return (dt.pm)
}

## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load.meta.choc <- function(file.name) {
  meta = read.delim(file.name,header=T,stringsAsFactors=F)
  make.global(meta)
  
  allnames = replace.col.names(names(meta),
                               c("SubjectID.blinded.","JCVI.Sample.ID","Sample_Tube_ID..revised.","Subject.s.Gender","Subject.s.YEAR.of.birth","Had.subject.fever.at.the.time.of.sample.collection."),
                               c("SubjectID","SampleID","SampleID.compound","Gender","YearOfBirth","Fever.descr"))
  
  names(meta) = allnames
  
  meta$SampleID = paste("CHOC",meta$SampleID,sep="")
  
  ## Format of SubjectID from Raja:
  ## The third letter in the letter set describes whether the sample is 
  ## from the patient (P) or sibling (S). The first set of numbers is 
  ## the date of consent without year. The second set of two numbers 
  ## links patient and sibling together. For example, ARP-0328-01 and 
  ## CRS-0329-01 are family set 01
  
  meta$FamilyID = substring(meta$SubjectID,10)
  
  meta$Sample.type = gsub(" ",".",meta$Sample.type)
  
  #Therapy.Status
  meta$Sample.type.1 = meta$Sample.type
  
  meta$Sample.type = unlist(apply(meta,
                                            1,
                                            function(row) {switch(row["Sample.type.1"],
                                                                  sibling ="sibling",
                                                                  patient.after.chemo="patient",
                                                                  patient.before.chemo="patient")}))
  meta$TherapyStatus = unlist(apply(meta,
                                              1,
                                              function(row) {switch(row["Sample.type.1"],
                                                                    sibling ="before.chemo",
                                                                    patient.after.chemo="after.chemo",
                                                                    patient.before.chemo="before.chemo")}))
  
  
  meta$Diagnosis[meta$Diagnosis=="#N/A"] = "Healthy"
  meta$Diagnosis[is.na(meta$Diagnosis) & meta$Sample.type=="sibling"] = "Healthy"
  #meta$Diagnosis[is.na(meta$Diagnosis) & meta$Sample.type=="patient"] = "leukemia"
  if(any(meta$Diagnosis=="Healthy" & meta$Sample.type=="patient")) {
    stop("Detected pateint with 'healthy' diagnosis value")
  }
  
  diag.sel = (meta$Diagnosis == "Healthy") | grepl("leukemia",meta$Diagnosis,ignore.case = T)
  
  meta = meta[diag.sel,,drop=F]
  
  subj.dual.sample.types = ddply(meta,"SubjectID",summarise,
                         dual.type=(length(unique(Sample.type))>1))
  
  if(any(subj.dual.sample.types$dual.type)) {
    stop(sprintf("Some subjects have multiple sample types - they should be either 
         only patient or sibling but not both: %s",
                 paste(subj.dual.sample.types$SubjectID[subj.dual.sample.types$dual.type],collapse=", ")))
  }
  
  meta$has.sibling = compute.has.sibling(meta)
  
  #meta = meta[meta$has.sibling,,drop=F]
  
  #meta$SampleID = factor(meta$SampleID)
  
  meta$Study.participation.status[str_blank(meta$Study.participation.status) | meta$Study.participation.status == "#N/A"] = "Active"
  
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
  meta$Specimen.Collection.Date = 
    as.Date.date.time(as.character(meta$Specimen.Collection.Date),
                      as.character(meta$Specimen.Collection.Time),
                      date.format = "%m/%d/%y")
  
  meta$age = as.numeric(
    difftime(
    meta$Specimen.Collection.Date
    ,
    as.POSIXct(paste(meta$YearOfBirth,"-06-01",sep=""),format="%Y-%m-%d"),
    units="days"
    )
  )/365
  
  ## we need to create this as ordered quantile to show
  ## in the right order on plots (and possible used
  ## in models that care about ordered factors)
  meta$age.quant = quantcut.ordered(meta$age)
  
  meta$Sampling.Visit = as.numeric(laply(as.character(meta$SampleID.compound),
                                function(x) unlist(strsplit(x,".v",fixed=T))[2]))
  meta$SampleID.1 = factor(paste(meta$SubjectID,meta$Sample.type.1,sep="."))
  meta$Sample.type = factor(meta$Sample.type)
  meta$Sample.type.1 = factor(meta$Sample.type.1)
  meta$TherapyStatus = factor(meta$TherapyStatus)
  meta$SampleID = factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  meta$Sample.type.1 = relevel(meta$Sample.type.1,"sibling")
  meta$Sample.type = relevel(meta$Sample.type,"sibling")
  meta$TherapyStatus = relevel(meta$TherapyStatus,"before.chemo")
  
  meta.aggr = join(meta,
                   ddply(meta,"SubjectID",summarise,
                         Antibiotic.Before.Therapy=any(ifelse(Sampling.Visit==1,as.logical(Antibiotic),F)),
                         Antibiotic.Min.Sampling.Visit=min(Sampling.Visit[Antibiotic])),
                   by="SubjectID",
                   type="inner")
  stopifnot(nrow(meta.aggr)==nrow(meta))

  meta = meta.aggr  
  rownames(meta) = meta$SampleID
  
  meta$Antibiotic.Before.Sample = (meta$Sampling.Visit >= meta$Antibiotic.Min.Sampling.Visit)
  meta$Antibiotic.Before.Sample[is.na(meta$Antibiotic.Before.Sample)] = T
  
  meta.aggr = join(meta,
                   ddply(meta,"SubjectID",summarise,
                         First.Specimen.Date=Specimen.Collection.Date[Sampling.Visit==1]),
                   by="SubjectID",
                   type="inner")
  stopifnot(nrow(meta.aggr)==nrow(meta))  
  
  meta = meta.aggr
  rownames(meta) = meta$SampleID
  
  meta$Days.In.Study = as.numeric(difftime(meta$Specimen.Collection.Date,meta$First.Specimen.Date,units="days"))
  meta$Days.In.Study.Bin = quantcut.ordered(meta$Days.In.Study,q=seq(0,1,by=0.20))
  
  report$add(qplot(x=Sampling.Visit,y=Days.In.Study,color=Antibiotic.Before.Therapy,shape=Sample.type,data=meta),
             caption="Relation between Sampling.Visits and dates")

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
  
  make.global(m_a)
  
  report$add.header("Summary of metadata variables")
  
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~Sample.type+TherapyStatus","~Antibiotic.Before.Therapy + Sample.type",
                        "~Antibiotic.Before.Therapy + Sample.type + Sampling.Visit",
                        "~Fever + Sample.type",
                        "~Sample.type+Sampling.Visit","~FamilyID","~Sample.type.1","~SubjectID",
                        "~SubjectID+Study.participation.status")
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
                                Sampling.Visit,
                                method="spearman"),
                       caption="Spearman RHO for age and Sampling.Visit")
    
  })
  with(meta[meta$Sample.type=="patient",],{
    report$add(cor.test(age,
                                Sampling.Visit,
                                method="spearman"),
                       caption="Spearman RHO for age and Sampling.Visit, patients only")
    report$add(glm(Antibiotic.Before.Therapy~Days.In.Study,family="binomial"),
               caption="How patients with initial anitbioitc treatment distributed between therapy periods")
  })
  report$add(ggplot(data=meta[meta$Sample.type=="patient",],
                    aes(x=Days.In.Study,
                        y=as.numeric(Antibiotic.Before.Therapy))) + 
               geom_point() + 
               stat_smooth(method="glm",family="binomial"),
             caption="How patients with initial anitbioitc treatment distributed between therapy periods")
  
  report$add(qplot(x=Sampling.Visit,y=Days.In.Study,color=Antibiotic.Before.Therapy,shape=Sample.type,data=meta),
             caption="Are Sampling.Visits evenly spaced for all patients (colored by initial antibiotic treatment)?")
  report$add(qplot(x=Sampling.Visit,y=Days.In.Study,color=Antibiotic.Before.Sample,shape=Sample.type,data=meta),
             caption="Are Sampling.Visits evenly spaced for all patients (colored by antibiotic treatment prior to current sample)?")

  #summary(glht(lm(age~Sampling.Visit,data=meta[meta$Sample.type=="patient",]),linfct="Sampling.Visit=0"))
  #summary(glht(lmer(age~Sampling.Visit+(Sampling.Visit|Sample.type),data=meta),linfct="Sampling.Visit=0"))
  report$add(ggplot(meta[meta$Sample.type=="patient",],aes(x=Sampling.Visit,y=age,color=Antibiotic.Before.Therapy))+
               geom_point()+
               stat_smooth(method="loess", se = T,degree=1,size=1),
             caption="Plot for age and Sampling.Visit with Loess trend line")

  aggr = ddply(meta,"SubjectID",summarise,
        Num.Sampling.Visits = length(Sampling.Visit),
        Sample.type=Sample.type[1],
        Antibiotic.Before.Therapy=Antibiotic.Before.Therapy[1],
        Antibiotic.Ever=any(Antibiotic),
        age=min(age),
        Max.Days.In.Study=max(Days.In.Study)
  )
  aggr$Frequency.of.Sampling.Visits = aggr$Num.Sampling.Visits/aggr$Max.Days.In.Study
  
  report$add(ggplot(aggr[aggr$Sample.type=="patient",],aes(x=Antibiotic.Before.Therapy,y=Num.Sampling.Visits,color=Max.Days.In.Study))+
               geom_jitter(position = position_jitter(width = .1)) + stat_summary(fun.data = "mean_cl_boot", 
                                            geom = "crossbar",
                              colour = "red", width = 0.2),
             caption="Plot for total number of Sampling.Visits")

  report$add(ggplot(aggr[aggr$Sample.type=="patient",],aes(x=Antibiotic.Before.Therapy,y=Frequency.of.Sampling.Visits,color=Max.Days.In.Study))+
               geom_jitter(position = position_jitter(width = .1)) + stat_summary(fun.data = "mean_cl_boot", 
                                            geom = "crossbar",
                                            colour = "red", width = 0.2),
             caption="Plot for total number of Sampling.Visits")
  report$add(ggplot(meta[meta$Sample.type=="patient",],aes(x=Days.In.Study,y=Sampling.Visit,color=Antibiotic.Before.Therapy,group=SubjectID))+
               geom_line() + scale_y_continuous(breaks=1:20),
             caption="Cumulative number of Sampling.Visits over therapy period")
  # + geom_text(aes(label=SubjectID,x=Specimen.Collection.Date,y=Sampling.Visit),size=4,angle=45)
  report$add(ggplot(meta[meta$Sample.type=="patient",],aes(x=Specimen.Collection.Date,y=Sampling.Visit,color=Antibiotic.Before.Therapy,group=SubjectID))+
               geom_line() + scale_y_continuous(breaks=1:20),
             caption="Cumulative number of Sampling.Visits vs date of collection colored by pre-study antibioitc use")
  report$add(ggplot(meta[meta$Sample.type=="patient",],aes(x=Specimen.Collection.Date,y=Sampling.Visit,color=Study.participation.status,group=SubjectID))+
               geom_line() + scale_y_continuous(breaks=1:20),
             caption="Cumulative number of Sampling.Visits vs date of collection colored by Study.participation.status")
  
}

new_emrevents <- function(...) {
  x = data.frame(...)
  class(x) <- append(class(x),"emrevents",0)
  return(x)
}

extract.by.name.emrevents <- function(x,names,conv=as.character,is.regex=F,ignore.case=T) {
  if(!is.regex) {
    res = x[x$Name %in% names,,drop=F]
  }
  else {
    res = x[grepl(names,x$Name,ignore.case=ignore.case),,drop=F]
  }
  res$Value = do.call(conv,list(res$Value))
  return (res)
}

extract.by.type.emrevents <- function(x,name) {
  if(name == "Temperature") {
    y = extract.by.name.emrevents(x,names="^TEMPERATURE.*",is.regex=T,conv=as.numeric)
    y = y[y$Value>=35,]
    y = y[y$Value<=45,]
  }
  else if(name == "Fever") {
    y = extract.by.type.emrevents(x,"Temperature")
    ## http://www.seattlechildrens.org/medical-conditions/symptom-index/fever/
    conv = c(TEMPERATUREORALC=37.8,TEMPERATUREAXILLARYC=37.2)
    y$Value = factor(ifelse(y$Value >= conv[y$Name], "Fever", "No.Fever"),levels=c("No.Fever","Fever"))
    y$Name = "Fever"
  }
  y = y[!is.na(y$Value),]
  return (y)
}

join.with.biom.attr.emrevents <- function(x,biom.attr) {
  biom.attr = biom.attr[!is.na(biom.attr$Sampling.Visit) & biom.attr$Sampling.Visit==1,]
  stopifnot(!any(duplicated(biom.attr$SubjectID)))
  x = merge(x,biom.attr[,c("SubjectID","First.Specimen.Date")],by="SubjectID")
  x$Days.In.Study = as.numeric(difftime(x$Date,x$First.Specimen.Date,units="days"))
  x$First.Specimen.Date = NULL
  x
}

load.choc.emrevents <- function(emr.measure.file) {
  events = read.delim(emr.measure.file,sep=",",stringsAsFactors=F,header=T)
  allnames = replace.col.names(names(events),
                               c("STUDYID",    "GROUP",      "EVENTDTTM",  "DISPLAYKEY", "DISPLAY",    "RESULT",     "RUNIT"),
                               c("SubjectID",    "Group",  "Date",  "Name", "Name.Descr",    "Value",     "Units"))
  
  names(events) = allnames  
  events$Date = as.POSIXct(events$Date,format="%m/%d/%y %H:%M")
  ret = new_emrevents(events)
  return (ret)
}

load.and.align.choc.emr.data <- function(emr.measure.file,biom.attr) {
  events = load.choc.emrevents(emr.measure.file)
  events = join.with.biom.attr.emrevents(events,biom.attr)
  return (events)
}

plot.emrevents <- function(x,m_a) {
  attr = m_a$attr[m_a$attr$SubjectID %in% x$SubjectID,]
  p = ggplot(x,aes(x=Days.In.Study,y=Value)) + 
    geom_point() +
    facet_wrap(~SubjectID,scale="free_x") +
    geom_vline(aes(xintercept=Days.In.Study), data=attr, color="black")
  return(p)
}

describe.choc.emr <- function(m_a) {
  events = load.and.align.choc.emr.data("emr.2015-05-15/microbiome_ce.csv",m_a$attr)
  y = extract.by.type.emrevents(events,"Fever")
  #print(plot.emrevents(y,m_a))
  #y = extract.by.type.emrevents(events,"Temperature")
  
  return (y)
}

## This function must generate a list with analysis tasks

gen.tasks.choc <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    #DEBUG:
    #taxa.levels = c(2,3,4,5,6,"otu")
    taxa.levels = c(2)
    
    descr = "All samples, no aggregation"
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="yap_2015-05-11/*.0.03.shared"
      cons.taxonomy.file="yap_2015-05-11/*.0.03.cons.taxonomy.seq.taxonomy"
      meta.file="CHOC_ALL_Samples_Metadata_April_2015.txt"
      load.meta.method=load.meta.choc
      load.meta.options=list()    
    })
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##some extremely young siblings will be outliers in microbiome compositon, remove
      m_a = subset.m_a(m_a,subset=(m_a$attr$age>=1.5))
#       m_a.pat = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient"))
#       m_a.sib = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="sibling"))
#       m_a.sib = subset.m_a(m_a.sib,subset=(!duplicated(m_a$attr$SubjectID)))
#       m_a = cbind.m_a(list(m_a.sib,m_a.pat),batch.attr="Sample.type",col.match=T)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient") | (!duplicated(m_a$attr$SubjectID)))
      
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      meta = m_a$attr
      meta.aggr = join(meta,
                            ddply(meta,"SubjectID",summarise,
                                  Sampling.Visit.max=max(Sampling.Visit),
                                  Sampling.Visit.min=min(Sampling.Visit),
                                  Sampling.Visit.1=any(Sampling.Visit==1),
                                  Sampling.Visit.2=any(Sampling.Visit==2)),
                            by="SubjectID",
                            type="inner")
      stopifnot(!any(is.na(meta.aggr$Sampling.Visit.max)) && 
                  nrow(meta.aggr)==nrow(meta))

      meta = meta.aggr
      rownames(meta) = meta$SampleID
      
      meta.aggr = join(meta,
                       ddply(meta,"FamilyID",summarise,
                             has.sibling.sample=("patient" %in% Sample.type) & ("sibling" %in% Sample.type)),
                       by="FamilyID",
                       type="inner")
      stopifnot(!any(is.na(meta.aggr$has.sibling.sample)) && 
                  nrow(meta.aggr)==nrow(meta))
      meta = meta.aggr
      rownames(meta) = meta$SampleID
      
      m_a$attr = meta
      
      ##As of 2014-11-05, there are only 6 samples at Sampling.Visit 5, and less in higher Sampling.Visits
      ##and their profiles look similar to Sampling.Visit 4
      #m_a = subset.m_a(m_a,subset=(m_a$attr$Sampling.Visit<=4))
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

    descr = "All samples (up to Sampling.Visit 4), no aggregation"    
    
    taxa.levels = c(2)
    
    do.summary.meta = T
    
    do.tests = F
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      return(m_a)
    }    
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = c("Days.In.Study")
      group.vars = c("Sample.type","Days.In.Study.Bin")
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=F
      do.network.features.combined = F
      do.divrich = c()
      
      do.plot.profiles.abund=T
      do.heatmap.abund=F
      
      
      divrich.task = within(divrich.task,{
        group.attr = NULL
        counts.glm.task = NULL
        do.plot.profiles = T
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Sample.type","Days.In.Study.Bin"),
                            c("Sample.type.Antibio.Before","Days.In.Study.Bin"),
                            c("Antibiotic.Before.Therapy","Sample.type.1"))
        feature.meta.x.vars=c("Days.In.Study")
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c("Sample.type","Sampling.Visit","Antibiotic.Before.Therapy")
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
      m_a = subset.m_a(m_a,subset=(m_a$attr$has.sibling.sample)) 
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
    
    descr = "Patients' samples at Sampling.Visits 1 (before therapy) and 2 (after therapy), only paired samples"
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Sample.type=="patient" 
                                   & m_a$attr$Sampling.Visit <= 2 
                                   & m_a$attr$Sampling.Visit.1
                                   & m_a$attr$Sampling.Visit.2))
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
        do.feature.meta=F
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
    
    main.meta.var = "Sampling.Visit"
    
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
               descr="Association with Sampling.Visit paired by subject")
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(main.meta.var,"Antibiotic.Before.Therapy"))
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic.Before.Therapy")
      })
      
    })
    
  })
  
  return (list(task1))
  #return (list(task1,task2,task2.1,task3,task3.1,task4))
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
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="atovtchi@jcvi.org",
                       title="Analysis of CHOC ALL 16S data",
                       incremental.save=F)

if(T) {
res = proc.project(
  task.generator.method=gen.tasks.choc
)

report$save()
}
