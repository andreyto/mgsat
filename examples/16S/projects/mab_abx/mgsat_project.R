
## Custom metadata loading function (define and pass to proc.project() when default
## implementation load.meta.default() is not sufficient)


load_meta_mouse_1 <- function(file.name) {
  library(data.table)
  ## We are calling default simple table loading here,
  ## but you can do something else. Result must have
  ## rownames() set to SampleID field
    meta = as.data.table(load.meta.default(file.name))
    meta[,TimePoint:=as.character(TimePoint)]
    meta[,TimePoint := factor(TimePoint,ordered = T)]
    meta[,Treatment := relevel(factor(Treatment),"Naive")]
    meta[,MouseID := factor(meta$MouseID)]
    meta[,TreatTime := sprintf("%s_%s",Treatment,TimePoint)]
    meta[,FullLabel := sprintf("%s_%s",MouseID,TimePoint)]
    meta[,TimePointLabel := sub("_"," ",TimePoint,fixed = T)]
    meta[,TimePointMouseLabel:=paste(TimePointLabel,MouseID)]
    
    meta = as.data.frame(meta)
    rownames(meta) = meta$SampleID
    return (meta)
}
  
## This function should carry out analysis specific to metadata fields by themselves, without
## relation to the abundance profiles. You can write it to do nothing (empty body).

summary.meta.mouse_1 <- function(m_a) {
  
  report$add.header("Summary of metadata variables after filtering samples")
  
  meta = m_a$attr

  report$add(summary(meta),caption="Summary of metadata variables")
  
  xtabs.formulas = list("~Treatment+TimePoint","~Treatment","~TimePoint",
                        "~MouseID")
  for(xtabs.formula in xtabs.formulas) {
    fact.xtabs = xtabs(as.formula(xtabs.formula),data=meta,drop.unused.levels=T)
    report$add(fact.xtabs,caption=paste("Sample cross tabulation",xtabs.formula))
    report$add.printed(summary(fact.xtabs))
  }

}


new_ordination.task <- function(main.meta.var,norm.method,label=NULL,size=NULL,
                                lines.args=NULL) {
  if(norm.method=="prop") {
    ord.method = "CCA"
    distance.0="bray"
  }
  else {
    ord.method = "RDA"
    distance.0="euclidean"
  }
  within(mgsat.16s.task.template$test.counts.task$ordination.task, {
    distance=distance.0
    ord.tasks = list(
      list(
        ordinate.task=list(
          method=ord.method
          ##other arguments to phyloseq:::ordinate
        ),
        plot.task=list(
          type="samples",
          color=main.meta.var,
          label = label,
          size = size,
          lines.args = lines.args
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
          color=main.meta.var,
          label = label,
          size = size,
          lines.args = lines.args
          ##other arguments to phyloseq:::plot_ordination
        )
      )          
    )
  })            
}

test_dist_to_baseline <- function(beta.dist,attr,unit.attr="MouseID",group.attr="Treatment",time.attr="TimePoint",time.baselev="Day_00") {
  beta.m = as.matrix(beta.dist)
  stopifnot(all(labels(beta.dist)==rownames(attr)))
  beta.m = beta.m[attr[time.attr]==time.baselev,]
  res = c()
  for(time.point in unique(attr[[time.attr]])) { 
    if(time.point != time.baselev) {
      beta.time = beta.m[,attr[[time.attr]]==time.point]
      attr.time.rows = attr[rownames(beta.time),,drop=F]
      unit.map.rows = attr.time.rows[[unit.attr]]
      unit.map.cols = attr[colnames(beta.time),unit.attr]
      dimnames(beta.time) = list(unit.map.rows,unit.map.cols)
      ind = rownames(beta.time)
      ind = as.matrix(data.frame(rows=ind,cols=ind))
      beta.val = beta.time[ind]
      beta.data = data.table(beta.dist=beta.val,group=attr.time.rows[[group.attr]])
      groups = unique(beta.data$group)
      groups.one = c()
      groups.two = c()
      pvals = c()
      for(group.one in groups) {
        for(group.two in groups) {
          if(group.two > group.one) {
            beta.data.pair = beta.data[beta.data$group %in% c(group.one,group.two)]
            mod = t.test(beta.dist~group,data = beta.data.pair)
            pvals = c(pvals,mod$p.value)
            groups.one = c(groups.one,group.one)
            groups.two = c(groups.two,group.two)
          }
        }
      }
      res_time = data.table(group.one=groups.one,group.two=groups.two,pvalue=pvals)
      res_time$time.point = time.point
      res = rbind(res,res_time)
      report$add(show.distr.group(beta.data$beta.dist,beta.data$group),
                 caption = sprintf("Empirical distribution at %s of %s groups dissimilarity to
                                         baseline %s within each experimental unit",time.point,group.attr,time.baselev))
    }
  }
  report$add.table(res,caption = sprintf("Unadjusted p-values from a t-test comparing %s groups pairwise by dissimilarity to
                                         baseline %s within each experimental unit",group.attr,time.baselev))
  res
}

## This function must generate a list with analysis tasks

gen.tasks.mouse_1 <- function() {
  
  task0 = within( mgsat.16s.task.template, {
    
    taxa.levels = c(2,5,6,"otu")
    #taxa.levels = c(6)
    
    descr = "All samples, no aggregation, no tests here, only plots"
    
    main.meta.var = "TreatTime"    
    
    read.data.task = within(read.data.task, {
      taxa.summary.file = NA
      otu.shared.file="otu.shared"
      cons.taxonomy.file="cons.taxonomy"
      meta.file="StudyDesign.txt"
      load.meta.method=load_meta_mouse_1
      load.meta.options=list()
      count.filter.options = list()    
      otu.count.filter.options=list()      
    })
    
    get.taxa.meta.aggr.base<-function(m_a) { 
      ##any aggregated attributes have to be computed here,
      ##after the available count samples have been joined,
      ##as opposed to in the load.meta() function.
      ## We already have all fields, so not aggregating anything here.
      
      return(m_a)
    }
    
    summary.meta.method=summary.meta.mouse_1
    
    test.counts.task = within(test.counts.task, {
      
      do.ordination=T
      do.network.features.combined=T      
      
      count.filter.feature.options = within(list(), {
        min_quant_mean_frac=0.25
        min_quant_incidence_frac=0.25
        #min_max=30
        min_mean=10
      })
      
      norm.count.task = within(norm.count.task, {
        method="norm.prop"
        method.args=list()
        #method="norm.ihs.prop"
        #method.args = list(theta=1)
        #method="norm.rlog.dds"
        #method.args=list(dds=NA) #signals to pull Deseq2 object
        #method="norm.clr"
        #method.args=list()
      })
      
      adonis.task = within(adonis.task, {
        
        dist.metr="bray"
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
        clustering_distance_rows="pearson"
        km.abund=0
        km.diversity=0
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="ihs.prop",label="SampleID")      
      
    })
    
  })
  
  
  task1 = within( task0, {

    descr = "All samples, no aggregation, no tests, only plots"    
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07","Day_14")))
      m_a$attr$Treatment = relevel(factor(m_a$attr$Treatment),"Naive")
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07") & m_a$attr$Treatment != "Naive"))
      return(m_a)
    }
    
    summary.meta.task = within(summary.meta.task, {
      meta.x.vars = c("TimePoint")
      group.vars = c("Treatment","TimePoint")
    })
    
    test.counts.task = within(test.counts.task, {
      
      do.divrich = c("otu")
      do.deseq2 = F
      do.adonis = T
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = "Treatment*TimePoint"
        })        
        do.plot.profiles = T
        n.rar.rep=40
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Treatment","TimePoint"),c("TimePoint","Treatment"),"TreatTime")
        feature.meta.x.vars=c("TimePoint")
        do.profile=T
        do.feature.meta=F
        #show.profile.task=within(show.profile.task,{
        #  facet_wrap_ncol=3
        #})
        
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs="Treatment*TimePoint",
               #strata="MouseID",
               ## strata will not work - mice are nested within treatments
               ## also, there are missing samples between time points
               descr="Association with treatment and timepoint, with interaction")
        )
        
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c("TreatTime","Treatment","TimePoint")
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="prop",label="FullLabel") 
      
    })
    
  })

  task1_baseline = within( task1, {
    
    descr = "Samples only on day 00, only plots"    
    
    main.meta.var = "Treatment"    
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00")))
      return(m_a)
    }
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = T
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      
      divrich.task = within(divrich.task,{
        group.attr = main.meta.var
        counts.glm.task = within(counts.glm.task,{
          formula.rhs = main.meta.var
        })        
        do.plot.profiles = T
        n.rar.rep=40
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(main.meta.var)
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
        #show.profile.task=within(show.profile.task,{
        #  facet_wrap_ncol=3
        #})
        
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs="Treatment",
               descr="Association with treatment")
        )
        
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c(main.meta.var)
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="prop",label="FullLabel") 
    })
  })

  task1_Bret_GordonConf = within( task1, {
    
    descr = "Samples only on days 00 and 07, only plots, aggregated for Bret's talk at Gordon Conference"    
    
    main.meta.var = "TreatTime"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07")))
      m_a$attr$Treatment = relevel(factor(m_a$attr$Treatment),"Naive")
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07") & m_a$attr$Treatment != "Naive"))
      return(m_a)
    }
    
    do.summary.meta = F
    
    do.tests = T
    
    taxa.levels = 6
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=F
      do.divrich = c()
      do.aggr.after.norm = taxa.levels
      
      aggr.after.norm.task = within(list(), {
        func = function(m_a,m_a.norm,m_a.abs,res.tests,...) 
        {
          m_a.norm = aggregate.by.meta.data.m_a(m_a.norm,group_col="TreatTime")
          m_a.norm$attr$Treatment = factor(m_a.norm$attr$Treatment,levels = c("Naive","LC10","R347","LVX","LZD","Van"))
          m_a.norm$attr$TimePointLabel = sub("_"," ",m_a.norm$attr$TimePoint,fixed = T)
          ## generate m_a as m_a.norm with equal library size to make sure that we
          ## are averaging after proprotions and not after summing raw counts downstream
          ## of this function
          m_a = m_a.norm
          m_a$count = m_a$count*10000
          return(list(m_a=m_a,m_a.norm=m_a.norm,m_a.abs=m_a))
        }
        ##possibly other arguments to func()
      })
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Treatment"))
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
        show.profile.task=within(show.profile.task,{
          facet_wrap_ncol=3
          geoms=c("bar_stacked")
          stat_summary.fun.y="identity"
          legend.title="Taxa"
          record.label="TimePointLabel"
          theme_font_size=1.1
        })
        
      })
      
    })
  })

  task1_means = within( task1, {
    
    descr = "Plots, aggregated into means over all samples in each cell of study design"    
    
    main.meta.var = "TreatTime"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07","Day_14")))
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07") & m_a$attr$Treatment != "Naive"))
      return(m_a)
    }
    
    do.summary.meta = F
    
    do.tests = T
    
    taxa.levels = 6
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=F
      do.ordination = T
      do.divrich = c()
      do.aggr.after.norm = taxa.levels
      
      aggr.after.norm.task = within(list(), {
        func = function(m_a,m_a.norm,m_a.abs,res.tests,...) 
        {
          m_a.norm = aggregate.by.meta.data.m_a(m_a.norm,group_col="TreatTime",count_aggr=function(x) Gmedian::Weiszfeld(x)$median,colwise = F)
          #m_a.norm$attr$Treatment = factor(m_a.norm$attr$Treatment,levels = c("Naive","LC10","R347","LVX","LZD","Van"))
          m_a.norm$attr$Treatment = relevel(factor(m_a.norm$attr$Treatment),"Naive")
          m_a.norm$attr$TimePointLabel = sub("_"," ",m_a.norm$attr$TimePoint,fixed = T)
          ## generate m_a as m_a.norm with equal library size to make sure that we
          ## are averaging after proprotions and not after summing raw counts downstream
          ## of this function
          m_a = m_a.norm
          m_a$count = m_a$count*10000
          return(list(m_a=m_a,m_a.norm=m_a.norm,m_a.abs=m_a))
        }
        ##possibly other arguments to func()
      })
      
      ordination.task = new_ordination.task("Treatment",norm.method="prop",label="FullLabel",
                                            lines.args = list(line.group="Treatment",
                                                              line.order="TimePoint")) 
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Treatment"))
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
        show.profile.task=within(show.profile.task,{
          facet_wrap_ncol=3
          geoms=c("bar_stacked")
          stat_summary.fun.y="identity"
          legend.title="Taxa"
          record.label="TimePointLabel"
          theme_font_size=1.1
        })
        
      })
      
    })
  })

  task1_means_pre_submission_2019_07_22 = within( task1_means, {
    
    descr = "Plots, aggregated into means over all samples in each cell of study design"    
    
    main.meta.var = "TreatTime"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07","Day_14")))
      m_a = subset.m_a(m_a,subset=(m_a$attr$Treatment != "KPN42"))
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07") & m_a$attr$Treatment != "Naive"))
      return(m_a)
    }
    
    do.summary.meta = F
    
    do.tests = T
    
    taxa.levels = c(2,6)
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=F
      do.ordination = T
      do.divrich = c()
      do.aggr.after.norm = taxa.levels

      count.filter.feature.options = within(list(), {
        min_quant_mean_frac=0.
        min_quant_incidence_frac=0.
        #min_max=30
        min_mean=10
      })
      
      aggr.after.norm.task = within(list(), {
        func = function(m_a,m_a.norm,m_a.abs,res.tests,...) 
        {
          m_a.norm = aggregate.by.meta.data.m_a(m_a.norm,group_col="TreatTime",count_aggr=function(x) Gmedian::Weiszfeld(x)$median,colwise = F)
          #m_a.norm$attr$Treatment = factor(m_a.norm$attr$Treatment,levels = c("Naive","LC10","R347","LVX","LZD","Van"))
          m_a.norm$attr$Treatment = relevel(factor(m_a.norm$attr$Treatment),"Naive")
          m_a.norm$attr$TimePointLabel = sub("_"," ",m_a.norm$attr$TimePoint,fixed = T)
          ## generate m_a as m_a.norm with equal library size to make sure that we
          ## are averaging after proprotions and not after summing raw counts downstream
          ## of this function
          m_a = m_a.norm
          m_a$count = m_a$count*10000
          return(list(m_a=m_a,m_a.norm=m_a.norm,m_a.abs=m_a))
        }
        ##possibly other arguments to func()
      })
    })
  })
  
  task1_color_like_means = within( task1_means, {
    
    descr = "Plots, all points, colored by Treatment"    
    
    main.meta.var = "TreatTime"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = task1_means$get.taxa.meta.aggr(m_a)
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07")))
      #m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07") & m_a$attr$Treatment != "Naive"))
      return(m_a)
    }
    
    do.summary.meta = F
    
    do.tests = T
    
    taxa.levels = 6
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.ordination = T
      do.heatmap.abund=F
      do.divrich = c()
      do.aggr.after.norm = c()
      
      aggr.after.norm.task = within(list(), {
        func = function(m_a,m_a.norm,m_a.abs,res.tests,...) 
        {
          m_a.norm$attr$Treatment = relevel(factor(m_a.norm$attr$Treatment),"Naive")
          m_a.norm$attr$TimePointLabel = sub("_"," ",m_a.norm$attr$TimePoint,fixed = T)
          return(list(m_a=m_a,m_a.norm=m_a.norm,m_a.abs=m_a))
        }
        ##possibly other arguments to func()
      })
      
      ordination.task = new_ordination.task("Treatment",norm.method="prop",label="FullLabel",
                                            lines.args = list(line.group="MouseID",
                                                              line.order="TimePoint")) 
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Treatment","TimePoint"))
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
        show.profile.task=within(show.profile.task,{
          facet_wrap_ncol=3
          geoms=c("bar_stacked")
          stat_summary.fun.y="identity"
          legend.title="Taxa"
          record.label="MouseID" #"TimePointMouseLabel"
          theme_font_size=1.1
          height = 600
          width = 800
        })
        
      })
      
    })
  })
      
  task1_focused = within( task1, {
    
    descr = "Samples only on days with Abx sampling, only plots"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07","Day_14"))) #,"Day_14"
      m_a$attr$Treatment = relevel(factor(m_a$attr$Treatment),"Naive")
      return(m_a)
    }
    
  })

  
  task1_focused_pre_submission_2019_07_22 = within( task1_focused, {
    
    descr = "Samples only on days with Abx sampling, only plots"
    
    taxa.levels = c("otu")
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07","Day_14"))) #,"Day_14"
      m_a = subset.m_a(m_a,subset=(m_a$attr$Treatment != "KPN42"))
      m_a$attr$Treatment = relevel(factor(m_a$attr$Treatment),"Naive")
      return(m_a)
    }

    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = F
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=F
      do.ordination = F
      do.heatmap.abund=F
      do.divrich = c("otu")
      do.aggr.after.norm = c()
      do.extra.method=c("otu")
      
      divrich.task = within(divrich.task, {
        do.plot.profiles = T
        do.plot.profiles.glm = F
        ## Computing beta-diversity matrix on multiple rarefications can take a while
        do.beta = T
        do.accum = F
        do.incidence = F
      })

      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(c("Treatment","TimePoint"))
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
        show.profile.task=within(show.profile.task,{
          facet_wrap_ncol=3
          geoms=c("line")
          show.samp.n = F
          line.legend.thickness=rel(5)
          theme_font_size = 0.6
          width = 600
          height = 400
        })
        
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          test_dist_to_baseline(beta.dist=get.beta.dist(res.tests),
                                attr=m_a$attr,
                                unit.attr="MouseID",
                                group.attr="Treatment",
                                time.attr="TimePoint",
                                time.baselev="Day_00")
        }
        
        norm.count.task.extra = within(norm.count.task, {
        })
        
      })
      
    })
    
  })
  
  task1_max_change = within( task1, {
    
    descr = "Timepoints with largest changes due to Abx treatment"
    
    do.summary.meta = T
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = get.taxa.meta.aggr.base(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$TimePoint %in% c("Day_00","Day_07"))) #,"Day_14"
      m_a$attr$Treatment = relevel(factor(m_a$attr$Treatment),"Naive")
      return(m_a)
    }

    summary.meta.task = within(summary.meta.task, {
      group.vars = c(main.meta.var)
    })
    
    test.counts.task = within(test.counts.task, {
      do.deseq2 = F
      do.adonis = T
      do.genesel = F
      do.stabsel = F
      do.glmer = F
      do.plot.profiles.abund=T
      do.heatmap.abund=T
      do.network.features.combined=F
      
      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("MouseID + %s",main.meta.var)
      })
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
      })
      
      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs="Treatment*TimePoint",
               #strata="MouseID",
               ## strata will not work - mice are nested within treatments
               ## there are also missing samples between time points
               descr="Association with treatment and timepoint, with interaction")
        )
        
      })
      
    })
    
  })

  
  task3 = within( task1_max_change, {
    
    descr = paste(descr,"Changes in 3902 between days 0 and 7")
    
    do.summary.meta = F
    
    do.tests = T
    
    get.taxa.meta.aggr<-function(m_a) { 
      m_a = task2$get.taxa.meta.aggr(m_a)
      m_a = subset.m_a(m_a,subset=(m_a$attr$Treatment=="3902"))
      return(m_a)
    }    
    
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
      
      plot.profiles.task = within(plot.profiles.task, {
        id.vars.list = list(main.meta.var)
        feature.meta.x.vars=c(main.meta.var)
        do.profile=T
        do.feature.meta=F
      })
      
      heatmap.combined.task = within(heatmap.combined.task, {
        attr.annot.names=c(main.meta.var)
      })
      
      ordination.task = new_ordination.task(main.meta.var,norm.method="prop",label="FullLabel") 
      
      
      genesel.task = within(genesel.task, {
        group.attr = main.meta.var
        genesel.param = within(genesel.param, {
          block.attr = "MouseID"
          type="paired"
          #replicates=0
        })
      })

      deseq2.task = within(deseq2.task, {
        formula.rhs = sprintf("MouseID + %s",main.meta.var)
      })

      stabsel.task = within(stabsel.task, {
        resp.attr=main.meta.var
      })
      
      adonis.task = within(adonis.task, {
        
        tasks = list(
          list(formula.rhs="TreatTime",
               strata="MouseID",
               descr="Association with treatment, within mouse strata")
        )
        
      })
      
      extra.method.task = within(extra.method.task, {
        
        func = function(m_a,m_a.norm,res.tests,norm.count.task.extra) {
          group.attr = main.meta.var
          test.dist.matr.within.between(m_a=m_a,
                                        group.attr=group.attr,
                                        block.attr="MouseID",
                                        n.perm=4000,
                                        dist.metr="manhattan",
                                        norm.count.task=norm.count.task.extra
          )
          require(d3heatmap)
          require(graphics)
          ord = order(m_a.norm$attr[,group.attr])
          labRow = as.character(m_a.norm$attr$SampleID[ord])
          p = d3heatmap(m_a.norm$count[ord,],
                        Rowv=F,
                        Colv=F,
                        labRow = labRow,
                        yaxis_width=label.size.points(labRow),
                        color="Reds")
          report$add.widget(p,caption="Dynamic Heatmap of normalized abundance")
          
        }
        
        norm.count.task.extra = within(norm.count.task, {
        })
        
      })
      
      
    })
    
  })

  #return (list(task1_focused_pre_submission_2019_07_22))
  return (list(task1_means_pre_submission_2019_07_22,task1_focused_pre_submission_2019_07_22))
  return (list(task1_focused_pre_submission_2019_07_22))
  return (list(task1,task1_max_change,task1_means,task1_color_like_means))
  return (list(task1_means,task1,task1_baseline,task2))
  return (list(task1,task1_focused,task2,task3))
}

threads=16

## number of cores to use on multicore machines
options(mc.cores=threads)
options(boot.ncpus=threads)
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

library("BiocParallel")
register(SnowParam(threads))

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"g_test.r",sep="/"),local=T)

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

evalsOptions("graph.output","svg")
#evalsOptions("width",600) #800 #1000
#evalsOptions("height",400) #640 #840

report <- PandocAT$new(author="John Doe",
                       title="Analysis of mouse mAb vs Abx treatment 16S data",
                       incremental.save=F)


res = proc.project(
  task.generator.method=gen.tasks.mouse_1
)

report$save()

