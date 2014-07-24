
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)

load.meta.example <- function(file_name,counts.row.names) {
  meta = load.meta.default(file_name,counts.row.names)
  meta.visit.max = join(meta,
                        ddply(meta,"SubjectID",summarise,visit.max=max(visit)),
                        by="SubjectID",
                        match="first")
  stopifnot(!any(is.na(meta.visit.max$visit.max)) && 
              nrow(meta.visit.max)==nrow(meta))
  meta = meta.visit.max
  return (meta)
}


## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.example <- function(taxa.meta) {
  
  report$add.header("Summary of metadata variables")
  
  
  xtabs.formulas = list("~Sample.type+DietStatus","~Sample.type+visit","~MatchedGroup","~Sample.type.1","~SubjectID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=taxa.meta$data,drop.unused.levels=T)
    report$add.table(fact.xtabs,show.row.names=T,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }
  
  with(taxa.meta$data,{
    report$add.printed(summary(aov(age~Sample.type)),
                       caption="ANOVA for age and sample type")
    report$add(qplot(Sample.type,age,geom="violin"),
               caption="Violin plot for age and sample type")
  })
  
  with(taxa.meta$data,{
    report$add.printed(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit")
    
  })
  with(taxa.meta$data[taxa.meta$data$Sample.type=="patient",],{
    report$add.printed(cor.test(age,
                                visit,
                                method="spearman"),
                       caption="Spearman RHO for age and visit, patients only")
    
  })
  
  report$add(ggplot(taxa.meta$data,aes(x=visit,y=age,color=Sample.type))+
               geom_point()+
               stat_smooth(method="loess", se = T,degree=1,size=1),
             caption="Plot for age and visit with Loess trend line")
  
}

## This function must generate a lits with analysis tasks

gen.tasks.example <- function(taxa.meta) {
  
  glmer.descr.tpl = "The random effect terms in the model formula describe that: 
  - intercept varies among matched groups
  and among individuals within matched groups (nested random effects);
  - when there are multiple observations (samples) per
  individual, we add a random effect for each observation to account
  for the overdispersion;
  The fixed effect is %s. The binomial family
  is used to build a set of univariate models, with each
  model describing the observed counts of one clade.
  P-values are estimated from the model under a null hypothesis
  of zero coefficients and a two-sided alternative. 
  Benjamini & Hochberg (1995) method is used 
  for multiple testing correction, and the significant clades
  are reported."
  
  task = within(
    list(),
{
  do.std.plots = T
  do.clade.meta=T
  do.profile=T
  
  do.heatmap = T
  do.std.tests = T
  
  do.stability=T
  do.tests=T
  do.genesel=T
  do.glmnet=T
  do.glmer=T
  do.adonis=T
  
  stability.steps = 600
  n.adonis.perm = 4000  
  
  stability.resp.attr = "NULL" 
  stability.model.family = "binomial"
  
  genesel.resp.attr = stability.resp.attr
  
  stability.transform.counts="ihs"
  
  adonis.tasks = list()
  
  glmer.task = list()
  
  plot.group = list()            
  
  clade.meta.x.vars=c()
  
  heatmap.task = list()
  
  
}
  )

task1 = within(
  task,
{
  
  descr = "All samples, no aggregation"
  
  do.stability=T
  do.tests=F
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  stability.resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = stability.resp.attr
  
  plot.group = list(
    c("DietStatus","Sample.type"),
    c("Sample.type","visit")
  )            
  
  clade.meta.x.vars=c("visit","age")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type.1","visit","age","Drug"),
    attr.row.labels="SampleID"
  )
  
}
)

task2 = within(
  task,
{
  
  descr = "Samples before Diet aggregated by SubjectID"
  
  do.clade.meta=F
  
  taxa.meta.aggr = aggregate_by_meta_data(subset(taxa.meta$data,DietStatus=="before.diet"),
                                          "SubjectID",
                                          taxa.meta$attr.names)
  stability.resp.attr = "Sample.type" 
  stability.model.family = "binomial"
  
  genesel.resp.attr = stability.resp.attr
  
  adonis.tasks = list(
    list(formula_rhs="Sample.type",
         strata=NULL,
         descr="Association with the patient/control status unpaired"),
    list(formula_rhs="Sample.type",
         strata="MatchedGroup",
         descr="Association with the patient/control status paired")
  )
  
  glmer.task = list(
    formula.rhs = "Sample.type + (1|MatchedGroup/SubjectID)",
    linfct=c("Sample.typecontrol = 0"),
    descr=sprintf(glmer.descr.tpl,"Sample.type")
  )
  
  plot.group = list(
    c("Sample.type")
  )            
  
  clade.meta.x.vars=c("visit")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type","Drug"),
    attr.row.labels="SampleID"
  )
  
}
)


task3 = within(
  task,
{
  
  descr = "Patients only, no aggregation"
  
  do.std.plots = F
  do.clade.meta=T
  do.profile=T
  
  do.heatmap = F
  do.std.tests = T
  
  do.stability=T
  do.tests=T
  do.genesel=T
  do.glmnet=T
  do.glmer=T
  do.adonis=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="patient")
  
  stability.resp.attr = "visit" 
  stability.model.family = "gaussian"
  
  genesel.resp.attr = "DietStatus" 
  
  adonis.tasks = list(
    list(formula_rhs="visit",
         strata=NULL,
         descr="Association with the visit number unpaired"),
    list(formula_rhs="visit",
         strata="SubjectID",
         descr="Association with the visit number paired by subject")
  )
  
  glmer.task = list(
    formula.rhs = "visit + (1|SubjectID) + (1|SampleID)",
    linfct=c("visit = 0"),
    descr=sprintf(glmer.descr.tpl,"visit")
  )
  
}
)

task4 = within(
  task,
{
  
  descr = "All controls and patients with two or more visits, no aggregation"
  
  do.stability=F
  do.tests=F
  do.clade.meta=F
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="control" | visit.max>=2)
  
  stability.resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = stability.resp.attr
  
  plot.group = list(
    c("Sample.type","visit")
  )            
  
  clade.meta.x.vars=c("visit","age")
  
  heatmap.task = list(
    attr.annot.names=c("Sample.type.1","visit","age","Drug"),
    attr.row.labels="SampleID"
  )
  
}
)

task5 = within(
  task,
{
  
  descr = "Patients only, no aggregation. Checking influence of the drug treatment confounder."
  
  do.stability=F
  do.tests=F
  do.clade.meta=T
  
  do.genesel=F
  do.glmnet=T
  
  taxa.meta.aggr = taxa.meta
  
  taxa.meta.aggr$data = subset(taxa.meta.aggr$data,Sample.type=="patient")
  
  stability.resp.attr = "Sample.type.1" 
  stability.model.family = "multinomial"
  
  genesel.resp.attr = stability.resp.attr
  
  plot.group = list(
    c("Drug","visit")
  )            
  
  clade.meta.x.vars=c("visit")
  
  heatmap.task = list(
    attr.annot.names=c("visit","Drug"),
    attr.row.labels="SampleID"
  )
  
}
)


#return (list(task1,task2,task3,task4,task5))
#return (list(task5))
return (list())
}




## number of cores to use on multicore machines
options(mc.cores=2)
options(boot.ncpus=2)
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

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=F)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="my_email@somewhere.com",
                       title="Example of 16S Analysis using MGSAT",
                       out.file.md="report.md",
                       incremental.save=T)

proc.project(
  taxa.summary.file=paste(MGSAT_SRC,"example.seq.taxsummary",sep="/"),
  meta.file=paste(MGSAT_SRC,"example_metadata.tsv",sep="/"),
  summary.meta.method=summary.meta.example,
  task.generator.method=gen.tasks.example,
  load.meta.method=load.meta.example
)

