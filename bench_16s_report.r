
cbind.m_a <- function(m_a.list,batch.attr,col.match=T) {
  cols.map = unlist(lapply(m_a.list,function(x) colnames(x$count)))
  batch.id.rows = unlist(lapply(m_a.list,function(x) (x$attr[,batch.attr])))
  batch.id.cols = unlist(lapply(m_a.list,function(x) {
    b.attr = unique(x$attr[,batch.attr])
    stopifnot(length(b.attr)<=1)
    rep(b.attr,ncol(x$count))
  }))
  names(cols.map) = paste(cols.map,ifelse(batch.id.cols!="",".",""),batch.id.cols,sep="")
  if(col.match) {
    cols = unique(cols.map)
  }
  else {
    cols.map[] = names(cols.map)
    cols = cols.map
  }
  m_a = foreach(m_a=m_a.list,
                .final=function(m_a.list) {
                  m_a = list()
                  m_a$count = do.call(rbind,lapply(m_a.list,function(x) {x$count}))
                  m_a$attr = do.call(rbind,lapply(m_a.list,function(x) {x$attr}))
                  m_a
                }) %do% {
                  batch.attr.val = unique(m_a$attr[,batch.attr])
                  stopifnot(length(batch.attr.val) <= 1)
                  if(batch.attr.val != "") {
                    sep="."
                  }
                  else {
                    sep=""
                  }
                  cols.keys = paste(colnames(m_a$count),batch.attr.val,sep=sep)
                  rows = paste(rownames(m_a$count),batch.attr.val,sep=sep)
                  count = matrix(0.,
                                 nrow=nrow(m_a$count),
                                 ncol=length(cols),
                                 dimnames=list(rows,cols))
                  count[,cols.map[cols.keys]] = m_a$count
                  m_a$count = count
                  rownames(m_a$attr) = rows
                  m_a
                }
  m_a$attr$SampleID = rownames(m_a$attr)
  m_a
}

load.ground.thruth.m_a <- function(abund.file,
                                   aggr.type=c("genus.abund","genus.otus","otu"),
                                   fake.sum=50000) {
  x = read.delim(abund.file,header=T,sep="\t")
  x = x[,c("Organism","Genus","Profile","WGS.Proportion.16S.Gene","ProfileID")]
  names(x)[4] = "Prop"
  x$Prop = round(x$Prop * fake.sum)
  x$Organism = sanitize.taxa.names(x$Organism)
  x$Genus = sanitize.taxa.names(x$Genus)
  if(aggr.type %in% c("genus.abund","genus.otus")) {
    aggr = ddply(x,c("Profile","ProfileID","Genus"),summarise,Abund=sum(Prop),Count=sum(Prop>0))
    if(aggr.type=="genus.abund") {
      value.var = "Abund"
    }
    else if(aggr.type=="genus.otus") {
      value.var = "Count"
    }
    count = acast(aggr,ProfileID~Genus,value.var=value.var)
  }
  else if(aggr.type == "otu") {
    count = acast(x,ProfileID~Organism,value.var="Prop")
  }
  attr = unique(x[,c("ProfileID","Profile")])
  rownames(attr) = attr$ProfileID
  attr = attr[rownames(count),]
  m_a = list(count=count,attr=attr)
  m_a
}


load.wl.cdhit.m_a <- function(data.file,meta.file=NULL,aggr.type=c("genus.abund","genus.otus","otu")) {
  x = read.delim(data.file,header=T,sep="\t",stringsAsFactors=F)
  ## extract all fields up to total field plus genus
  x = x[,c(names(x)[seq(match("total",names(x))-1)],"family","genus")]
  
  ## make taxa names compatible with Mothur output
  mask.rec = (x$genus == "Incertae Sedis")
  x$genus[mask.rec] = paste(x$family[mask.rec],"incertae sedis",sep="_")
  mask.rec = grepl("uncultured.*",x$genus,ignore.case=T) | (!grepl("\\S+",x$genus))
  x$genus[mask.rec] = paste("Unclassified",x$family[mask.rec],sep="_")
  
  x$genus = gsub("\\s","_",x$genus)
  x$genus = gsub("-","_",x$genus) # Escherichia-Shigella
  x$genus[x$genus == "Clostridium"] = "Clostridium_sensu_stricto"
  
  x$OTU = paste(x$OTU,x$genus,sep=".")
  rownames(x) = x$OTU
  x$OTU = NULL
  x$family = NULL
  
  if(aggr.type=="otu") {
    x$genus = NULL
    x = t(x)
  }
  else {
    aggr.func = switch(aggr.type,
                       genus.abund = sum,
                       genus.otus = function(y) sum(y>0)
    )
    x = aggregate(x[,!("genus" == names(x))],list(genus = x$genus),aggr.func)
    
    rownames(x) = x$genus
    x$genus = NULL
    x = t(x)
  }
  
  meta = load.meta.default(meta.file)
  return (merge.counts.with.meta(x,meta))
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
                       title="Report on benchmarking of 16S annotation pipelines",
                       incremental.save=F)

cl = start.cluster.project()

report$add.descr("16S annotation count matrices from different pipelines and parameter combinations are
                 compared against the expected ground truth.")

meta.file.samples = "refdata/bench_bei_meta.txt"

with(mgsat.16s.task.template,{
  runs.files = list(
    list(
      RunID="YAP.comm",
      run.descr="YAP currently deployed in common area",
      read.data.task=within(read.data.task, {
        taxa.summary.file = "yap/bei/v13/common/7913675baaa4dd38c5ad2c19b10bdda5.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
        otu.shared.file="yap/bei/v13/common/5f44040ee0ae9f533d83e48deded4bac.files_x1.sorted.0.03.shared"
        cons.taxonomy.file="yap/bei/v13/common/dce008a0990dc67e8148b3a31fb0bb02.files_x1.sorted.0.03.cons.taxonomy"
        taxa.summary.file.otu = "yap/bei/v13/common/0211bc387719c7e0a4dc9d4ea0eccee9.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
        meta.file=meta.file.samples
      })
    ),  
    #   list(
    #     RunID="YAP.base",
    #     read.data.task=within(read.data.task, {
    #       taxa.summary.file = "yap/bei/v13/baseline/bd568ff2897e717e1ab6c948f9509344.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
    #       otu.shared.file="yap/bei/v13/baseline/21c401833a6ba4d011974fa46f05e476.files_x1.sorted.0.03.shared"
    #       cons.taxonomy.file="yap/bei/v13/baseline/ff4e455254bb6759e890e170151a4e37.files_x1.sorted.0.03.cons.taxonomy"
    #       taxa.summary.file.otu = "yap/bei/v13/baseline/1baf5634efbb3aa711408089eda1faf5.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"      
    #       meta.file=meta.file.samples
    #     })
    #   ),
    #   list(
    #     RunID="YAP.mc.ch1",
    #     read.data.task=within(read.data.task, {
    #       taxa.summary.file = "yap/bei/v13/make.contigs.chim_full_ref/5526056cfd0ce087ba07cb3377c90b4f.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
    #       otu.shared.file="yap/bei/v13/make.contigs.chim_full_ref/60d6c77f675da1ab6bada2c013e8a8aa.files_x1.sorted.0.03.shared"
    #       cons.taxonomy.file="yap/bei/v13/make.contigs.chim_full_ref/3ca693e09ea5fb1d1e769bb5e7cef101.files_x1.sorted.0.03.cons.taxonomy"
    #       taxa.summary.file.otu = "yap/bei/v13/make.contigs.chim_full_ref/83324aef0f7627d8195b55efa56e5fb7.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
    #       meta.file=meta.file.samples
    #     })
    #   ),  
    list(
      RunID="YAP.upd",
      run.descr="YAP updated to use Mothur make.contigs and Mothur reference DB v.10",    
      read.data.task=within(read.data.task, {
        taxa.summary.file = "yap/bei/v13/make.contigs.chim_gold/d4db63db7e566943c9ba8277f0d6b253.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary"
        otu.shared.file="yap/bei/v13/make.contigs.chim_gold/57782e4f219e8f2212aa5e3ba3f91452.files_x1.sorted.0.03.shared"
        cons.taxonomy.file="yap/bei/v13/make.contigs.chim_gold/4021347ead1b4c49a7723f682c9362a1.files_x1.sorted.0.03.cons.taxonomy"
        taxa.summary.file.otu = "yap/bei/v13/make.contigs.chim_gold/e7cdd0d2622b10ea7d2a58666eb713f4.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary"
        meta.file=meta.file.samples
      })
    ),    
    list(
      RunID="Mothur.01",
      run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP by Andrey (reference-based chimera removal)",
      read.data.task=within(read.data.task, {
        taxa.summary.file = "bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
        otu.shared.file="bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
        cons.taxonomy.file="bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
        taxa.summary.file.otu = "bei/v13/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
        meta.file=meta.file.samples
      })
    ),
    #   list(
    #   RunID="Mothur.Sarah.01",
    #   read.data.task=within(read.data.task, {
    #     taxa.summary.file = "sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
    #     otu.shared.file="sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
    #     cons.taxonomy.file="sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
    #     taxa.summary.file.otu = "sarah/mothur.even.2015-03-17/HM782D.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
    #     meta.file=meta.file.samples
    #   })
    #  )
      list(
      RunID="Mothur.Sarah.02",
      run.descr="Mothur 1.34.4 with Mothur reference DB v.10 ran through SOP by Sarah (dataset-based chimera removal)",
      read.data.task=within(read.data.task, {
        taxa.summary.file = "sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.pds.wang.tax.summary"
        otu.shared.file="sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.shared"
        cons.taxonomy.file="sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.taxonomy"
        taxa.summary.file.otu = "sarah/mothur.2015-03-19/stability.trim.contigs.good.unique.good.filter.unique.precluster.pick.an.unique_list.0.03.cons.tax.summary"
        meta.file=meta.file.samples
      })
     )
    
  )
  
  ground.truth.WGS.file = "refdata/sarah.2015-03-08/derived/bei_abund_v.5.txt"
  ground.truth.run.id = "Ground.Truth.WGS"
  ground.truth.run.descr = "Mock samples sequenced with WGS, reads mapped to reference to compute coverage
which is multiplied by 16S copy number and normalized to proportions to get expected ground truth organism
abundance in the biological sample. Proportions multiplied by a large constant to imitate 16S read counts
so this dataset can be used in diversity and abundance estimates together with actual 16S annotation runs."
  ProfileID = "HM782D"
  drop.taxa = c("Nothing") #c("Clostridium_sensu_stricto") #c("Helicobacter") #c("Streptococcus")
  
  norm.method.basic.default = "norm.prop"
  true.taxa.only = F
  compositional.transform = F
  
  vegdist.method = "euclidian"
  norm.count.task = within(test.counts.task$norm.count.task, {
    method = "norm.clr"
    method.args = list(offset=0,tol=0.005)
    drop.features=list("other")
  })
  
  report.section = report$get.section()
  
  for(task.aggr in list(
    list(aggr.type="genus.abund",
         aggr.descr="Genus relative abundance",
         norm.method.basic = norm.method.basic.default),
    list(aggr.type="genus.otus",
         aggr.descr="OTU counts per genus",
         norm.method.basic = "ident"),
    list(aggr.type="genus.otus",
         aggr.descr="OTU relative counts per genus",
         norm.method.basic = norm.method.basic.default),    
    list(aggr.type="otu",
         aggr.descr="OTU abundance",
         norm.method.basic = norm.method.basic.default)
  )) {
    
    with(task.aggr,{
      
      report$add.header(paste("Annotation level:",aggr.descr),
                        report.section=report.section,sub=T)
      
      
      count.basis = "seq"
      taxa.level = 6
      col.match = T
      do.divrich = F
      do.abund = T
      do.summary.meta = F
      
      if(aggr.type == "genus.otus") {
        count.basis = "otu"
        compositional.transform = F
        taxa.level = 6
        do.divrich = F
        do.summary.meta = T
      }
      else if(aggr.type == "otu") {
        taxa.level = "otu"
        col.match = F
        do.summary.meta = T
        do.divrich = T
        do.abund = F
      }
      
      if(!compositional.transform) {
        vegdist.method = "manhattan"
        norm.count.task$method = "ident"
      }
      
      
      test.counts.task = within(test.counts.task, {
        divrich.task = within(list(),{
          n.rar.rep=400
          is.raw.count.data=T
          group.attr = "SampleID"
          counts.glm.task = NULL
          beta.task = NULL
          counts.genesel.task = NULL
          do.plot.profiles = T
          do.incidence=F
          do.abundance=T
          do.rarefy=F
          do.accum=F
        })
        
        plot.profiles.task = within(plot.profiles.task, {
          id.vars.list = list(c("SampleID"))
          feature.meta.x.vars=NULL
          do.profile=T
          do.feature.meta=F
          show.profile.task=within(show.profile.task, {
            geoms=c("bar")
            sqrt.scale=T
          })
        })
        
      })

      report$add.header("Loading annotation files")
      report$push.section(report.section)
      
      runs.data = lapply(runs.files,
                         function(run.files) {
                           with(run.files,{
                             read.data.task$count.basis = count.basis
                             m_a = do.call(read.data.method,
                                           c(
                                             list(taxa.level=taxa.level),
                                             read.data.task
                                           )
                             )
                             m_a$attr$RunID = RunID
                             m_a$attr$run.descr = run.descr
                             m_a
                           })
                         }
      )
      
      m_a = cbind.m_a(runs.data,batch.attr="RunID",col.match=col.match)
      make.global(m_a)
      m_a.gt = load.ground.thruth.m_a(ground.truth.WGS.file,aggr.type=aggr.type)
      make.global(m_a.gt)
      m_a.gt$attr$SampleID = m_a.gt$attr$ProfileID
      m_a.gt$attr$RunID = ground.truth.run.id
      m_a.gt$attr$run.descr = ground.truth.run.descr
      m_a.gt$attr$IdSfx = m_a.gt$attr$RunID
      make.global(m_a.gt)
      
      m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.00005.txt",
                                 meta.file=meta.file.samples,
                                 aggr.type=aggr.type)
      m_a.wl$attr$RunID = "WL.00005"
      m_a.wl$attr$run.descr = "Weizhong's CD-HIT-OTU pipeline that discarded OTUs below relative abundance of 0.00005"
      m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
      m_a.wl.1 = m_a.wl
      m_a.wl = load.wl.cdhit.m_a(data.file="wl-cdhit/bei/v13/wl-cdhit.bei.v13.0.0001.txt",
                                 meta.file=meta.file.samples,
                                 aggr.type=aggr.type)
      m_a.wl$attr$RunID = "WL.0001"
      m_a.wl$attr$run.descr = "Weizhong's CD-HIT-OTU pipeline that discarded OTUs below relative abundance of 0.0001"
      m_a.wl$attr$IdSfx = m_a.wl$attr$RunID
      m_a.wl.2 = m_a.wl
      
      m_a$attr$IdSfx = ""
      #m_a.prop = norm.count.m_a(m_a,method=norm.method.basic)
      #make.global(m_a.prop)
      m_a.abs = cbind.m_a(list(m_a.gt,m_a,m_a.wl.1,m_a.wl.2),batch.attr="IdSfx",col.match=col.match)
      m_a.abs$attr$IdSfx = NULL
            
      attr.rep = unique(m_a.abs$attr[,c("RunID","run.descr")])
      attr.rep$run.descr = gsub("\n"," ",attr.rep$run.descr)
      report$add.table(attr.rep,caption="Description of annotation runs")
      
      report$pop.section()
      
      report$add.header("Data analysis")
      report$push.section(report.section)
        
      m_a.abs = subset.m_a(m_a.abs,subset=(m_a.abs$attr$ProfileID==ProfileID))
      
      #if(aggr.type == "otu") {
      #  m_a.abs = subset.m_a(m_a.abs,subset=!(m_a.abs$attr$RunID==ground.truth.run.id))
      #}
      
      m_a.abs = count.filter.m_a(m_a.abs,drop.zero=T)
      
      make.global(m_a.abs)
      
      if(do.summary.meta) {
        
        summary.meta.task = within(summary.meta.task, {
          group.vars = NULL
          show.sample.totals=T
          show.sample.means=F
          sub.report=F
        })
        
        do.call(report.sample.count.summary,c(
          list(m_a.abs),
          summary.meta.task
        )
        )
      }
      
      
      if(do.divrich) {
        for(do.rarefy in c(F,T)) {
          divrich.task = test.counts.task$divrich.task
          divrich.task$do.rarefy = do.rarefy
          divrich.task$extra.header = ifelse(do.rarefy,
                                             "with rarefication to common depth",
                                             "without rarefication to common depth")
          tryCatchAndWarn({ 
            do.call(mgsat.divrich.report,
                    c(list(m_a.abs,
                           plot.profiles.task=test.counts.task$plot.profiles.task),
                      divrich.task)
            )
          })
        }
      }
      
      if(do.abund) {
        m_a.abs = subset.m_a(m_a.abs,select.count=!(colnames(m_a.abs$count) %in% drop.taxa))
        if(true.taxa.only) {
          m_a.abs = subset.m_a(m_a.abs,
                               select.count=colnames(m_a.abs) %in% colnames(m_a.gt$count))
        }
        m_a.abs = count.filter.m_a(m_a.abs,drop.zero=T)
        #m_a.prop = subset.m_a(m_a.prop,select.count=seq(9))
        m_a.prop = norm.count.m_a(m_a.abs,method=norm.method.basic)
        
        m_a.prop <- norm.count.report(m_a.prop,
                                      res.tests=NULL,
                                      descr="Extra normalization",
                                      norm.count.task)
        
        m_a.prop = subset.m_a(m_a.prop,select.count=colSums(abs(m_a.prop$count))>0)
        make.global(m_a.prop)
        
        if(compositional.transform) {
          geoms.profiles=c("bar")
        }
        else {
          geoms.profiles=c("bar_stacked","bar")
        }
        
        plot.profiles.task = within(test.counts.task$plot.profiles.task, {
          show.profile.task=within(show.profile.task, {
            geoms=geoms.profiles
          })
        })
        
        do.call(plot.profiles,
                c(list(m_a=m_a.prop,
                       feature.order=NULL,
                       feature.descr=sprintf("Plots of %s",aggr.descr)),
                  test.counts.task$plot.profiles.task
                )
        )
        
        report$add.header("Distance measures between samples")
        report$push.section(report.section)
        
        if(norm.method.basic == "norm.prop" && norm.count.task$method == "ident") {
          report$add.table(as.matrix(dist.js(m_a.prop$count)),show.row.names=T,
                           caption="Jensen-Shannon distance")
          report$add.table(as.matrix(vegdist(sqrt(m_a.prop$count),method = "euclidian")/sqrt(2)),show.row.names=T,
                           caption="Hellinger distance")
        }
        report$add.table(as.matrix(vegdist(m_a.prop$count,method = vegdist.method)),show.row.names=T,
                         caption=sprintf("%s distance",vegdist.method))
        report$pop.section()
      }
      
      report$pop.section()
    })
    report$pop.section()
  }
})

stop.cluster.project(cl)

report$save()
