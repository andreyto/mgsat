library(ggplot2)
#library(knitr)
library(pander)

panderOptions("table.style","rmarkdown")
panderOptions("table.split.table",180)
panderOptions("table.alignment.default","left")
panderOptions("evals.messages",F)
#panderOptions("graph.fontsize",10)
evalsOptions("cache",F)
evalsOptions("cache.mode","environment")
#evalsOptions("output",c("result"))
evalsOptions("output",c("all"))
evalsOptions("graph.output","png")
evalsOptions("graph.unify",F)
evalsOptions("width",640)
evalsOptions("height",480)
evalsOptions("res",75)
evalsOptions("hi.res",T)
evalsOptions("hi.res.width",800)

make.global <- function(var) {
  assign(deparse(substitute(var)),var,envir=globalenv()) 
}

mget.stack<-function(x,ifnotfound) {
  pars = lapply(rev(sys.parents()),sys.frame)
  #print(pars)
  #print(sys.parent())
  #print(sys.frames())
  for(fr in pars) {
    y = tryCatch(get(x,envir=fr,inherits=F),error=function(x) NULL)
    #print(paste("y=",y,"ifnotfound=",ifnotfound))
    if(!is.null(y)) {
      #print(paste("y=",y,"ifnotfound=",ifnotfound))
      return (y)
    }
  }
  return (ifnotfound)
}

get.default.section<-function() {
  x = new.env(parent=emptyenv())
  x$path = list(0)
  return (x)
}

clone.report.section<-function(x) {
  return (as.environment(as.list(x, all.names=TRUE)))
}

get.report.section<-function(default=get.default.section()) {
  y = mget.stack("report.section",ifnotfound=0)
  if(!identical(y,0)) return (y)
  return (default)
}

push.report.section<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section(default=NULL)
    if(is.null(x)) {
      x = get.default.section()
    }
    else {
      x = clone.report.section(x)
    }
  }
  x$path[[length(x$path)+1]] = 0
  return (x)
}

pop.report.section<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section(default=NULL)
    if(is.null(x)) {
      x = get.default.section()
    }
  }
  l = length(x$path)
  if(l>1)
    x$path = x$path[1:l-1]
  return (x)
}

incr.report.section<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section(default=NULL)
    if(is.null(x)) {
      x = get.default.section()
    }
  }
  last = length(x$path)
  if(last>0) {x$path[[last]] = x$path[[last]] + 1}
  else {x$path[[1]] = 1}
  return (x)
}

format.report.section<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()
  }
  return (paste("\\(",paste(x$path,sep="",collapse="."),"\\)",sep=""))
}



pandoc.as.printed.return <- function(x,attrs="") {
  x = capture.output(print(x))
  paste0('\n', repChar('`', 7), 
         ifelse(attrs == '', '', sprintf('{%s}', attrs)), 
         '\n', paste(x, 
                     collapse = '\n'), 
         '\n', repChar('`', 7), '\n')
}

PandocAT <- setRefClass('PandocAT', contains = "Pandoc", 
                        fields = list('sections' = 'list')
                        )

## private service method - should be called whenever an element is
## appended to the .self$body
PandocAT$methods(priv.append.section = function() {
  .self$sections[[length(.self$sections)+1]] = get.report.section()$path
})

PandocAT$methods(priv.format.caption = function(caption,section=NULL) {
  if(!is.null(caption)) {
    return (paste(format.report.section(section),caption))
  }
  else {
    return (caption)
  }
})

PandocAT$methods(add = function(x,new.paragraph=T,caption=NULL,...) {
  
  timer           <- proc.time()
  par_fr <- parent.frame()
  if(!identical(.GlobalEnv, par_fr)) {
    env             <-  list2env(as.list(par_fr, all.names=TRUE),parent=parent.frame(2))
  }
  else {
    env = NULL
  }
  res             <- evals(deparse(match.call()[[2]]),env=env,...)
  if(new.paragraph) {
    .self$add.p("")
  }

  if(!is.null(caption)) {
    .self$add.p(.self$priv.format.caption(caption))
  }
  .self$body      <- c(.self$body,res)
  .self$priv.append.section()
  .self$proc.time <- .self$proc.time + as.numeric(proc.time() - timer)[3]
  
})

PandocAT$methods(add.list = function(x,...) {
  return(.self$add(as.list(x),...))
})

PandocAT$methods(get.section = function(...) {
  return(get.report.section(...))
})

PandocAT$methods(incr.section = function(...) {
  return(incr.report.section(...))
})

PandocAT$methods(push.section = function(...) {
  return(push.report.section(...))
})

PandocAT$methods(pop.section = function(...) {
  return(pop.report.section(...))
})

PandocAT$methods(add.table = function(x,show.row.names=F,
                                      echo=T,
                                      caption=NULL,
                                      split.tables=180,style="rmarkdown",...) {
  caption = .self$priv.format.caption(caption)
  if(is.null(x) || nrow(x)==0) {
    if(!is.null(caption)) {
      .self$add.p(caption)
    }
    return(.self$add.p("Empty dataset"))
  }
  if(!show.row.names) {
    rownames(x) <- NULL
  }
  tbl_p <- pandoc.table.return(x,split.tables=split.tables,style=style,caption=caption,...)
  if(echo) {
    print(tbl_p)
  }
  return(.self$add.p(tbl_p))
})

PandocAT$methods(add.vector = function(x,name=NULL,
                                       show.row.names=T,
                                       caption=NULL,
                                       wrap.vals=T,
                                       ...) {
  if(is.null(x) || length(x)==0) {
    if(!is.null(caption)) {
      .self$add.p(.self$priv.format.caption(caption))
    }
    return(.self$add.p("Empty dataset"))
  }
  if(wrap.vals) {
    x = wrap(x)
  }
  y = data.frame(x=x)
  if(!is.null(name)) {
    names(y) <- c(name)
  }
  if(is.null(names(x))) {
    show.row.names = F
  }
  if(show.row.names) {
    row.names(y) <- names(x)
  }
  return(.self$add.table(y,caption=caption,show.row.names=show.row.names,...))
})

PandocAT$methods(add.p = function(x,rule=F,echo=T,...) {
  if(rule) {
    .self$priv.append.section()
    .self$add.paragraph(pandoc.horizontal.rule.return())
  }
  if(echo) {
    cat(format.report.section(),x,"\n")
  }
  .self$priv.append.section()
  return(.self$add.paragraph(x,...))
})

PandocAT$methods(add.descr = function(x,...) {
  .self$add.p(pandoc.strong.return(x),...)
})

PandocAT$methods(add.package.citation = function(x,...) {
  .self$add.p(capture.output(print(citation(x),style="text")))
})

PandocAT$methods(add.printed = function(x,caption=NULL,echo=T,...) {
  if(!is.null(caption)) {
    .self$add.p(.self$priv.format.caption(caption))
  }
  return(.self$add.p(pandoc.as.printed.return(x,...),echo=echo))
})

PandocAT$methods(add.header = function(x,level=NULL,section.action="incr",echo=T,...) {
  
  report.section = switch(section.action,
                 incr=incr.report.section(),
                 push=incr.report.section(),
                 keep=get.report.section())
  if (is.null(level)) {
    level = min(6,length(report.section))
  }
  x = paste(format.report.section(report.section),x)
  .self$add.p(pandoc.header.return(x,level=level,...),echo=echo)
  
  if(section.action=="push") {
    report.section = push.report.section()
  }
  return (report.section)
  
})

PandocAT$methods(save = function(f) {
  
  if (missing(f))
    f <- tempfile('pander-', getwd())
  fp    <- sprintf('%s.md', f)
  timer <- proc.time()
  
  ## create pandoc file
  cat(pandoc.title.return(.self$title, .self$author, .self$date), file = fp)
  f_sections = .self$sections
  make.global(f_sections)
  ##sort by section lexicographically, using a stable sort
  sect_ord = sort.list(
    unlist(lapply(.self$sections,function(x) paste(x,sep="",collapse="."))),
    method="shell")
  lapply(.self$body[sect_ord], function(x) cat(paste(pandoc.return(x$result), collapse = '\n'), file = fp, append = TRUE))
    
})

test_report.sections<-function() {
  
  report <- PandocAT$new("noone@email.com","test")
  
  get.pandoc.section.3<-function() {
    print("get.pandoc.section.3")
    #print(incr.report.section())
    report.section = report$add.header("get.pandoc.section.3",section.action="push")
    report.section = report$add.header("get.pandoc.section.3.1")
    report$add.table(data.frame(A=c("a","b"),B=c(1,2)),caption="Table")
    return (report.section)
  }
  
  get.pandoc.section.2<-function() {
    print("get.pandoc.section.2")  
    #report.section = incr.report.section()
    print(get.report.section())
    report.section = report$add.header("get.pandoc.section.2",section.action="push")
    get.pandoc.section.3()
    #report.section = report$pop.section()
    get.pandoc.section.3()
    report.section = report$pop.section()
    get.pandoc.section.3()
  }
  
  get.pandoc.section.1<-function() {
    print("get.pandoc.section.1")
    report.section = report$get.section()
    report$add.header("get.pandoc.section.1")
    #report.section = incr.report.section()
    #print(pandoc.section)
    print(get.report.section())
    report$add.header("get.pandoc.section.1")
    report$add(plot(x <- sort(rnorm(47))),caption="Figure")
    get.pandoc.section.2()
    #report.section = push.report.section()
    print(get.report.section())
    report$add.header("get.pandoc.section.4")
  }
  
  get.pandoc.section.1()
  
  report$save("test_report")
  Pandoc.convert("test_report.md",format="html",open=F)
  
}


my_p <- function(x, wrap = panderOptions('p.wrap'), sep = panderOptions('p.sep'), copula = panderOptions('p.copula'), limit = Inf, keep.trailing.zeros = panderOptions('keep.trailing.zeros')){
  
  stopifnot(is.vector(x))
  stopifnot(all(sapply(list(wrap, sep, copula), function(x) is.character(x) && length(x) == 1)))
  nm <- names(x)
  attributes(x) <- NULL
  x.len <- length(x)
  if (x.len == 0)
    return('')
  stopifnot(x.len <= limit)
  
  ## prettify numbers
  if (is.numeric(x)) {
    
    x <- round(x, panderOptions('round'))
    x <- format(x, trim = TRUE, digits = panderOptions('digits'), decimal.mark = panderOptions('decimal.mark'))
    
    ## optionally remove trailing zeros
    if (!keep.trailing.zeros)
      x <- sub('(?:(\\..*[^0])0+|\\.0+)$', '\\1', x)
    
  }
  
  if (x.len == 1)
    wrap(x, wrap)
  else if (x.len == 2)
    paste(wrap(x, wrap), collapse = copula)
  else
    paste0(paste(wrap(head(x, -1), wrap), collapse = sep), copula, wrap(tail(x, 1), wrap))
}

#library(R.utils)
#reassignInPackage("p", pkgName="pander", my_p)


