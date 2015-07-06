
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
  meta = meta[!is.na(meta$Specimen.Collection.Date),]
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
  
  pat.vis.2 = meta[meta$Sample.type=="patient" & meta$Sampling.Visit==2,]
  stopifnot(!any(duplicated(pat.vis.2$FamilyID)))
  
  pat.vis.2 = cbind(pat.vis.2[,c("FamilyID"),drop=F],PatientID=pat.vis.2$SubjectID,First.Therapy.Date=pat.vis.2$Specimen.Collection.Date)
  
  meta = merge(meta,pat.vis.2,by="FamilyID",all.x=T)
  
  imp.first.ther = meta$First.Specimen.Date + mean(meta$First.Therapy.Date-meta$First.Specimen.Date,na.rm=T)
  
  meta$First.Therapy.Date[is.na(meta$First.Therapy.Date)] = imp.first.ther[is.na(meta$First.Therapy.Date)]
  
  rownames(meta) = meta$SampleID
  
  meta$Days.In.Therapy = as.numeric(difftime(meta$Specimen.Collection.Date,meta$First.Therapy.Date,units="days"))
  
  ## some siblings had their virst visit around the second visit of the patient, so ther day in therapy are positive at this point; reset
  ## days in therapy of all siblings to the mean of all patients at their visit 1
  meta$Days.In.Therapy[meta$Sample.type == "sibling"] = mean(meta$Days.In.Therapy[meta$Sample.type == "patient" & meta$Sampling.Visit==1],na.rm=T)
  stopifnot(all(meta$Days.In.Therapy[meta$Sample.type == "sibling"]<0))
  
  qt = quantile(meta$Days.In.Therapy[meta$Days.In.Therapy>=0],na.rm=T)
  qt[length(qt)] = qt[length(qt)] + 1
  meta$Days.In.Therapy.Bin = cut(meta$Days.In.Therapy,
                                 c(min(meta$Days.In.Therapy,na.rm=T)-1,
                                   0,
                                   qt[-1]),
                                 right=F,
                                 ordered_result=T)
  lev = levels(meta$Days.In.Therapy.Bin)
  levels(meta$Days.In.Therapy.Bin) = paste("St.",seq_along(lev),"_",lev,sep="")
  meta$Study.Stage = as.character(meta$Days.In.Therapy.Bin)
  meta$Study.Stage[meta$Sample.type=="sibling"] = "St.0_sibling"
  meta$Study.Stage = ordered(meta$Study.Stage,levels=c("St.0_sibling",levels(meta$Days.In.Therapy.Bin)))
  
  ## ignore antibiotic status in siblings - this field is for faceted plots
  meta$Sample.type.Antibio.Before = factor(with(meta,
                                                ifelse(Sample.type=="patient",
                                                       paste(as.character(Sample.type),
                                                             ifelse(Antibiotic.Before.Therapy,
                                                                    "Antibiotic.Before",
                                                                    "No.Antibiotic.Before"),
                                                             sep="."),
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
    report$add(glm(Antibiotic.Before.Therapy~Days.In.Therapy,family="binomial"),
               caption="How patients with initial anitbioitc treatment distributed between therapy periods")
  })
  report$add(ggplot(data=meta[meta$Sample.type=="patient",],
                    aes(x=Days.In.Therapy,
                        y=as.numeric(Antibiotic.Before.Therapy))) + 
               geom_point() + 
               stat_smooth(method="glm",family="binomial"),
             caption="How patients with initial anitbioitc treatment distributed between therapy periods")
  
  report$add(qplot(x=Sampling.Visit,y=Days.In.Therapy,color=Antibiotic.Before.Therapy,shape=Sample.type,data=meta),
             caption="Are Sampling.Visits evenly spaced for all patients (colored by initial antibiotic treatment)?")
  report$add(qplot(x=Sampling.Visit,y=Days.In.Therapy,color=Antibiotic.Before.Sample,shape=Sample.type,data=meta),
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
               Max.Days.In.Therapy=max(Days.In.Therapy)
  )
  aggr$Frequency.of.Sampling.Visits = aggr$Num.Sampling.Visits/aggr$Max.Days.In.Therapy
  
  report$add(ggplot(aggr[aggr$Sample.type=="patient",],aes(x=Antibiotic.Before.Therapy,y=Num.Sampling.Visits,color=Max.Days.In.Therapy))+
               geom_jitter(position = position_jitter(width = .1)) + stat_summary(fun.data = "mean_cl_boot", 
                                                                                  geom = "crossbar",
                                                                                  colour = "red", width = 0.2),
             caption="Plot for total number of Sampling.Visits")
  
  report$add(ggplot(aggr[aggr$Sample.type=="patient",],aes(x=Antibiotic.Before.Therapy,y=Frequency.of.Sampling.Visits,color=Max.Days.In.Therapy))+
               geom_jitter(position = position_jitter(width = .1)) + stat_summary(fun.data = "mean_cl_boot", 
                                                                                  geom = "crossbar",
                                                                                  colour = "red", width = 0.2),
             caption="Plot for total number of Sampling.Visits")
  report$add(ggplot(meta[meta$Sample.type=="patient",],aes(x=Days.In.Therapy,y=Sampling.Visit,color=Antibiotic.Before.Therapy,group=SubjectID))+
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

extract.by.group.emrevents <- function(x,groups,is.regex=F,ignore.case=T) {
  if(!is.regex) {
    res = x[x$Group %in% groups,,drop=F]
  }
  else {
    res = x[grepl(groups,x$Group,ignore.case=ignore.case),,drop=F]
  }
  return (res)
}

extract.by.name.emrevents <- function(x,names,conv=as.character,is.regex=F,ignore.case=T,invert=F) {
  if(invert) op = function(x) (!x)
  else op = function(x) (x)
  if(!is.regex) {
    res = x[op(x$Name %in% names),,drop=F]
  }
  else {
    res = x[op(grepl(names,x$Name,ignore.case=ignore.case)),,drop=F]
  }
  res$Value = do.call(conv,list(res$Value))
  return (res)
}

extract.other.infections.raw.emrevents <- function(x) {
  y = extract.by.group.emrevents(x,groups="MICROBIOLOGY",is.regex=F)
  y = extract.by.name.emrevents(y,names="^CDIFF.*",is.regex=T,invert=T)
  y
}


mean.na.rm <- function(x) mean(x,na.rm=T)

extract.by.type.emrevents <- function(x,name,round.into.days=F,round.aggr.fun=NULL) {
  if(name == "BMI") {
    ##both weight and height are not always recorded on the same day
    w = extract.by.name.emrevents(x,names="WEIGHT",is.regex=F,conv=as.numeric)
    w = round.numeric.into.days.emrevents(w,aggr.fun=mean.na.rm)
    h = extract.by.name.emrevents(x,names="HEIGHT",is.regex=F,conv=as.numeric)
    h = round.numeric.into.days.emrevents(h,aggr.fun=mean.na.rm)
    h = h[,c("SubjectID","Days.In.Therapy","Value")]
    h$Value.H = h$Value
    h$Value = NULL
    y = merge(w,h,by=c("SubjectID","Days.In.Therapy"))
    y$Value = y$Value/(y$Value.H/100)**2
    y$Name = "BMI"
    y$Name.Descr = y$Name
  }
  else if(name == "Temperature") {
    y = extract.by.name.emrevents(x,names="^TEMPERATURE.*",is.regex=T,conv=as.numeric)
    y = y[y$Value>=35,]
    y = y[y$Value<=45,]
    if(round.into.days) {
      y = round.numeric.into.days.emrevents(y,aggr.fun=round.aggr.fun)
    }
  }
  else if(name == "Fever") {
    y = extract.by.type.emrevents(x,"Temperature")
    ## http://www.seattlechildrens.org/medical-conditions/symptom-index/fever/
    #conv = c(TEMPERATUREORALC=37.8,TEMPERATUREAXILLARYC=37.2)
    #Keri's response on Oral and my guess on axillary
    conv = c(TEMPERATUREORALC=38.0,TEMPERATUREAXILLARYC=37.8)
    y$Value = factor(ifelse(y$Value >= conv[y$Name], "Fever.Pos", "Fever.Neg"),levels=c("Fever.Neg","Fever.Pos"))
    y$Name = name
    if(round.into.days) {
      y = round.binary.into.days.emrevents(y)
    }
  }
  else if(name == "Diarrhea") {
    y = extract.by.name.emrevents(x,names="STOOLCONSISTENCY")
    y$Value = factor(ifelse(grepl("Liquid",y$Value,ignore.case=T), "Diarrhea.Pos", "Diarrhea.Neg"))
    y$Value = relevel(y$Value,"Diarrhea.Neg")
    y$Name = name
    if(round.into.days) {
      y = round.binary.into.days.emrevents(y)
    }
  }
  else if(name == "Neutropenia") {
    y = extract.by.name.emrevents(x,names="^ABSOLUTENEUTROPHILCOUNT.*",is.regex=T,conv=as.numeric)
    y$Value = factor(ifelse(y$Value < 500, "Neutropenia.Pos", "Neutropenia.Neg"))
    y$Value = relevel(y$Value,"Neutropenia.Neg")
    y$Name = name
    if(round.into.days) {
      y = round.binary.into.days.emrevents(y)
    }
  }
  else if(name == "C.Diff") {
    y = extract.by.name.emrevents(x,names="^CDIFF.*",is.regex=T)
    y$Value = factor(ifelse(grepl("^POS.*",y$Value,ignore.case=T), "C.Diff.Pos", "C.Diff.Neg"))
    y$Value = relevel(y$Value,"C.Diff.Neg")
    y$Name = name
    if(round.into.days) {
      y = round.binary.into.days.emrevents(y)
    }
  }
  else if(name == "Other.Infections.Test") {
    y = extract.other.infections.raw.emrevents(x)
    y$Value = factor(ifelse(grepl("^POS.*",y$Value,ignore.case=T), "Other.Infections.Test.Pos", "Other.Infections.Test.Neg"))
    y$Value = relevel(y$Value,"Other.Infections.Test.Neg")
    y$Name = name
    if(round.into.days) {
      y = round.binary.into.days.emrevents(y)
    }
  }
  
  y = y[!is.na(y$Value),]
  return (y)
}

join.with.biom.attr.emrevents <- function(x,biom.attr) {
  biom.attr = biom.attr[!is.na(biom.attr$Sampling.Visit) & biom.attr$Sampling.Visit==1,]
  stopifnot(!any(duplicated(biom.attr$SubjectID)))
  x = merge(x,biom.attr[,c("SubjectID","First.Therapy.Date")],by="SubjectID")
  x$Days.In.Therapy = as.numeric(difftime(x$Date,x$First.Therapy.Date,units="days"))
  x$First.Therapy.Date = NULL
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

is.two.level.factor <- function(x) {
  is.factor(x) && (length(levels(x)) == 2)
}

get.presence.level <- function(x) {
  lev = levels(x)
  stopifnot(length(lev)==2)
  return (lev[2])
}

binary.aggregate.by.presence.emrevents <- function(x,keys,min.presence=1) {
  x$Value = factor(x$Value)
  lev = get.presence.level(x$Value)
  all.lev = levels(x$Value)  
  x = ddply(x,keys,here(summarise),Value=(sum(Value==lev)>=min.presence))
  x$Value = factor(as.numeric(x$Value),labels=all.lev)
  return (x)  
}

round.binary.into.days.emrevents <- function(x,min.presence=1) {
  x$Days.In.Therapy = round(x$Days.In.Therapy)
  x = binary.aggregate.by.presence.emrevents(x,c("SubjectID","Days.In.Therapy"),min.presence=min.presence)
  return (x)
}

round.numeric.into.days.emrevents <- function(x,aggr.fun=NULL) {
  if(is.null(aggr.fun)) {
    aggr.fun = mean.na.rm
  }
  x$Days.In.Therapy = round(x$Days.In.Therapy)
  x = ddply(x,c("SubjectID","Days.In.Therapy"),here(summarise),Value=aggr.fun(Value))
  return (x)
}

two.level.to.numeric <- function(x) {
  lev = get.presence.level(x)
  as.numeric(x == lev)
}

two.level.to.numeric.emrevents <- function(x) {
  x$Value = two.level.to.numeric(x$Value)
  x
}

rolling.window.smooth.emrevents <- function(x,lower=-5,upper=5) {
  y = x
  y$Value = rolling.window.irreg(data=x,bylist="SubjectID",dates="Days.In.Therapy",
                                 target="Value",lower=-5,upper=5,incbounds=T,
                                 stat=mean,na.rm=T,cores=1)
  y
}

plot.emrevents <- function(x,m_a,str.ylab=NULL,facet.by.subject=T,to.binomial=F,
                           smooth=F,smooth.lower=-5,smooth.upper=5,jitter=T,stat.smooth=T,
                           mixed.model=F) {
  ret = list()
  val.two.level = is.two.level.factor(x$Value)
  val.numeric = is.numeric(x$Value)
  if(smooth) {
    to.binomial = T
  }
  if(to.binomial) {
    if(val.two.level) {
      x = two.level.to.numeric.emrevents(x)
    }
  }
  p = ggplot(x,aes(x=Days.In.Therapy,y=Value))
  if(smooth) {
    y = rolling.window.smooth.emrevents(x,lower=smooth.lower,upper=smooth.upper)
    p = p + geom_point(data=y,color="blue")
    color.data = "cyan"
    alpha.data = 0.7
  }
  else {
    color.data = "black"
    alpha.data = 1 
  }
  if(jitter) {
    p = p + geom_jitter(position = position_jitter(height = .03),color=color.data,alpha=alpha.data)
  }
  else {
    p = p + geom_point(color=color.data,alpha=alpha.data)
  }
  if(!is.null(ylab)) {
    p = p + ylab(str.ylab)
  }
  if(facet.by.subject) {
    attr = m_a$attr[m_a$attr$SubjectID %in% x$SubjectID,]    
    p = p + facet_wrap(~SubjectID,scale="free_x") + geom_vline(aes(xintercept=Days.In.Therapy), data=attr, color="green")
  }
  if(stat.smooth) {
    method = "glm"
    if(to.binomial && val.two.level) {
      family = "binomial"
      lmer.method = "glmer"
    }
    else if(val.numeric) {
      family = "gaussian"
      lmer.method = "lmer"
    }
    else {
      family = NA
    }
    if(!is.na(family)) {
      p = p + stat_smooth(method=method,family=family)
      if(mixed.model) {
        require(lme4)
        if(facet.by.subject) re.form = NULL
        else re.form = NA
        if(lmer.method=="glmer") {
          m = glmer(Value ~ 1 + Days.In.Therapy + (1|SubjectID) + (0+Days.In.Therapy | SubjectID), data=x, family=family)
        }
        else {
          m = lmer(Value ~ 1 + Days.In.Therapy + (1|SubjectID) + (0+Days.In.Therapy | SubjectID), data=x)
          require(afex)
          ret$mixed.model.test = mixed(m,data=x)
        }
        ret$mixed.model=m
        ret$mixed.model.family=family
        pr.x = cbind(x,ModelValue=predict(m,type="response",re.form=re.form))
        p = p + geom_point(aes(x=Days.In.Therapy,y=ModelValue),data=pr.x,color="red")
      }
    }
  }
  ret$main.plot = p
  return (ret)
}

xtabs.report <- function(form,data) {
  fact.xtabs = xtabs(as.formula(form),data=data,drop.unused.levels=T)
  report$add(fact.xtabs,caption=paste("Sample cross tabulation",form))
  report$add.table(as.data.frame(fact.xtabs),caption=paste("Sample cross tabulation as table",form),show.row.names=T)
  report$add.printed(summary(fact.xtabs))
}

describe.emrevents <- function(x,m_a,str.ylab,mixed.model=F,smooth=F) {
  two.level = is.two.level.factor(x$Value)  
  if(two.level) {
    to.binomial = T
  }
  else {
    to.binomial = F
  }
  
  res = plot.emrevents(x,m_a,
                       str.ylab=str.ylab,
                       facet.by.subject=T,
                       to.binomial=to.binomial,
                       smooth=smooth,
                       stat.smooth=T,
                       mixed.model=mixed.model)
  
  width = 1200
  height = width*0.8
  report$add(res$main.plot,
             caption = str.ylab,
             width=width,
             height=height,
             hi.res.width=width,
             hi.res.height=height)
  
  if(mixed.model && !is.null(res$mixed.model)) {
    if(res$mixed.model.family=="binomial") {
      require(arm)
      report$add(binnedplot(predict(res$mixed.model,type="response",re.form=NA),
                            resid(res$mixed.model,type="response"),
                            xlab="Expected Values", ylab="Average residual"),
                 caption="Binned plot without random effects")
    }
  }
  
  res = plot.emrevents(x,m_a,
                       str.ylab=str.ylab,
                       facet.by.subject=F,
                       to.binomial=to.binomial,
                       smooth=smooth,
                       stat.smooth=T,
                       mixed.model=mixed.model)
  
  report$add(res$main.plot,
             caption = str.ylab)
  
  if(mixed.model && !is.null(res$mixed.model)) {
    if(res$mixed.model.family=="binomial") {
      require(arm)
      report$add(binnedplot(predict(res$mixed.model,type="response",re.form=NULL),
                            resid(res$mixed.model,type="response"),
                            xlab="Expected Values", ylab="Average residual"),
                 caption="Binned plot with random effects")
    }
    report$add.printed(res$mixed.model)
    report$add.printed(summary(res$mixed.model))
    anv = anova(res$mixed.model)
    report$add.printed(anv)
    if(!is.null(res$mixed.model.test)) {
      report$add.printed(res$mixed.model.test)
    }
    report$add(plot(res$mixed.model))
  }
  
}

describe.choc.emr <- function(m_a,mixed.model=F) {
  events = load.and.align.choc.emr.data("emr.2015-05-15/microbiome_ce.encod.csv",m_a$attr)
  make.global(events)
  
  y = extract.by.type.emrevents(events,"Fever",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="Fever", mixed.model=mixed.model) 
  
  y = extract.by.type.emrevents(events,"Diarrhea",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="Diarrhea", mixed.model=mixed.model) 
  
  y = extract.by.type.emrevents(events,"Neutropenia",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="Neutropenia", mixed.model=mixed.model) 
  
  y = extract.by.type.emrevents(events,"C.Diff",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="C.Diff", mixed.model=mixed.model)   
  
  y = extract.by.type.emrevents(events,"BMI",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="BMI", mixed.model=mixed.model)   
  
  y = extract.by.type.emrevents(events,"Other.Infections.Test",round.into.days=T)
  describe.emrevents(y,m_a=m_a,str.ylab="Other.Infections.Test", mixed.model=mixed.model)   
  y.raw = extract.by.type.emrevents(events,"Other.Infections.Test",round.into.days=F)
  y.raw.aggr = binary.aggregate.by.presence.emrevents(y.raw,c("SubjectID","Name.Descr"),1)
  xtabs.report("~Name.Descr+Value",y.raw.aggr)
  #return (y)
}

new_ordination.task <- function(main.meta.var,norm.method) {
  if(norm.method=="prop") {
    ord.method = "CCA"
  }
  else {
    ord.method = "RDA"
  }
  within(mgsat.16s.task.template$test.counts.task$ordination.task, {
    distance="euclidean"
    ord.tasks = list(
      list(
        ordinate.task=list(
          method=ord.method
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
          method=ord.method,
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
}


## This function must generate a list with analysis tasks

gen.tasks.choc <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    #DEBUG:
    taxa.levels = c(2,3,4,5,6,"otu")
    #taxa.levels = c(2,4,6,"otu")
    
    norm.method = "ihs.prop"
    
    descr = "All samples, no aggregation"
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="yap_2015-05-11/*.0.03.shared"
      cons.taxonomy.file="yap_2015-05-11/*.0.03.cons.taxonomy.seq.taxonomy"
      meta.file="CHOC_ALL_Samples_Metadata_April_2015.txt"
      load.meta.method=load.meta.choc
      load.meta.options=list()
      count.filter.options = list()    
      #       count.filter.options = within(list(), {
      #         drop.unclassified=T
      #         min_quant_mean_frac=0.25
      #         min_quant_incidence_frac=0.25
      #         #min_max=30
      #         min_mean=10
      #       })
      
      otu.count.filter.options=list()
      
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
      
      do.ordination=T
      do.network.features.combined=T
      #do.divrich = c()
      
      count.filter.feature.options = within(list(), {
        drop.unclassified=T
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
        if(norm.method=="prop") {
          method="norm.prop"
          method.args=list()
        }
        else if(norm.method=="ihs.prop") {
          method="norm.ihs.prop"
          method.args = list(theta=1000)
        }
        else if(norm.method=="clr") {
          method="norm.clr"
          method.args=list(offset=1)
        }
        else if(norm.method=="alr") {
          method="norm.alr"
          method.args=list(offset=1,divcomp="other")
        }        
        else if(norm.method=="fr") {
          method="norm.fr"
          method.args=list(offset=0.1,divcomp=1,offset.is.prop=T)
        }                
        else {
          stop(sprintf("Unplanned norm.method %s",norm.method))
        }
      })
      
      divrich.task = within(divrich.task,{
        do.beta = F
        do.accum = F
        do.incidence = F
      })
            
      genesel.task = within(genesel.task, {
        ## use common m_a.norm
        norm.count.task = NULL
      })
      
      adonis.task = within(adonis.task, {
        if(norm.method=="prop") {
          dist.metr="bray"
        }
        else {
          dist.metr="euclidean"          
        }
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
        if(norm.method=="prop") {
          clustering_distance_rows="manhattan"
        }
        else {
          clustering_distance_rows="pearson"
        }
        km.abund=0
        km.diversity=0
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method)
      
    })
    
  })
  
  
  task1 = within( task0, {
    
    descr = "All samples, no aggregation, only visualization"    
    
    #taxa.levels = c(2)
    
    main.meta.var = "Study.Stage"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      return(m_a)
    }    
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = c("Days.In.Therapy")
      group.vars = c("Sample.type","Days.In.Therapy.Bin")
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      #do.network.features.combined = F
      #do.divrich = c()      
      
      divrich.task = within(divrich.task,{
        group.attr = NULL
        counts.glm.task = NULL
        do.plot.profiles = T
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Sample.type"),c("Study.Stage"))
        #                            c("Antibiotic.Before.Therapy","Study.Stage"))
        feature.meta.x.vars=c("Days.In.Therapy")
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Days.In.Therapy.Bin","Antibiotic.Before.Therapy")
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        #attr.annot.names=c(main.meta.var,"age","Antibiotic.Before.Therapy")
        attr.annot.names=c(main.meta.var,"Sample.type.1")
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method)
      
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
      
    })
    
  })
  
  
  task2 = within( task0, {
    
    main.meta.var = "Sample.type"
    
    descr = "Patient/Sibling samples before therapy aggregated by SubjectID"
    
    do.summary.meta = F
    
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
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c(main.meta.var,"age")
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method)
      
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
      do.deseq2 = T
      do.adonis = T
      do.genesel = T
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.heatmap.abund=F
      do.extra.method = taxa.levels
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("FamilyID+%s",main.meta.var)
      })
      
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
          if(norm.method=="prop") {
            dist.metr = "bray"
          }
          else {
            dist.metr = "euclidean"
          }
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="Sample.type",
                                        block.attr="FamilyID",
                                        n.perm=4000,
                                        dist.metr=dist.metr,
                                        col.trans="ident",
                                        norm.count.task=norm.count.task.extra
          )
        }
        norm.count.task.extra = within(norm.count.task, {
          #method="norm.prop"
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
        formula.rhs = sprintf("SubjectID+%s",main.meta.var)
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
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c(main.meta.var,"SubjectID")
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method)
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          if(norm.method=="prop") {
            dist.metr = "bray"
          }
          else {
            dist.metr = "euclidean"
          }
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr="TherapyStatus",
                                        block.attr="SubjectID",
                                        n.perm=8000,
                                        dist.metr=dist.metr,
                                        col.trans="ident",                                        
                                        norm.count.task=norm.count.task.extra
          )
        }
        norm.count.task.extra = within(norm.count.task, {
          #method="norm.prop"
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
    
    main.meta.var = "Days.In.Therapy"
    main.meta.var.bin = "Days.In.Therapy.Bin"
    
    descr = "Patients samples over time in study"
    
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
        group.attr = main.meta.var.bin
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })
        do.plot.profiles = T
      })
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("SubjectID+%s",main.meta.var)
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
               descr=sprintf("Association with %s paired by subject",main.meta.var))
        )
        
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c(),c(main.meta.var.bin),c(main.meta.var.bin,"Antibiotic.Before.Therapy"))
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=T
      })
      
      heatmap.abund.task = within(heatmap.abund.task,{
        attr.annot.names=c(main.meta.var,"Antibiotic.Before.Therapy")
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c(main.meta.var.bin,"SubjectID")
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method)
      
    })
    
  })

  
  task4.infl.up = task4
  task4.infl.up$get.taxa.meta.aggr = function(m_a) { 
    m_a = task4$get.taxa.meta.aggr(m_a)
    m_a = subset.m_a(m_a,subset=(m_a$attr$Study.Stage <= "St.3_[25.6,83.2)"))
    return(m_a)
  }

  task4.infl.down = task4
  task4.infl.down$get.taxa.meta.aggr = function(m_a) { 
    m_a = task4$get.taxa.meta.aggr(m_a)
    m_a = subset.m_a(m_a,subset=(m_a$attr$Study.Stage >= "St.3_[25.6,83.2)"))
    return(m_a)
  }
    
  #return (list(task1))
  return (list(task1,task2,task2.1,task3,task3.1,task4,task4.infl.up,task4.infl.down))
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
source(paste(MGSAT_SRC,"rolling_window_irreg.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=F)

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
