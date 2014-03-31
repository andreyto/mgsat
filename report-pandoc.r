library(ggplot2)
#library(knitr)
library(pander)

panderOptions("table.style","rmarkdown")
panderOptions("table.split.table",180)
panderOptions("table.alignment.default","left")
panderOptions("evals.messages",F)

evalsOptions("cache",F)
evalsOptions("cache.mode","environment")
#evalsOptions("output",c("result"))
evalsOptions("output",c("all"))
evalsOptions("graph.unify",F)
evalsOptions("res",75)
evalsOptions("hi.res",T)

make.global <- function(var) {
  assign(deparse(substitute(var)),var,envir=globalenv()) 
}

pandoc.as.printed.return <- function(x,attrs="") {
  x = capture.output(print(x))
  paste0('\n', repChar('`', 7), 
         ifelse(attrs == '', '', sprintf('{%s}', attrs)), 
         '\n', paste(x, 
                     collapse = '\n'), 
         '\n', repChar('`', 7), '\n')
}

PandocAT <- setRefClass('PandocAT', contains = "Pandoc", fields = list('tag' = 'character'))


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
    .self$add.p(caption)
  }
  .self$body      <- c(.self$body, res)
  .self$proc.time <- .self$proc.time + as.numeric(proc.time() - timer)[3]
  
})

PandocAT$methods(add.table = function(x,show.row.names=F,
                                      echo=T,
                                      caption=NULL,
                                      split.tables=180,style="rmarkdown",...) {
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
  return(.self$add.paragraph(tbl_p))
})

PandocAT$methods(add.vector = function(x,name=NULL,
                                       show.row.names=T,
                                       caption=NULL,
                                       ...) {
  if(is.null(x) || length(x)==0) {
    if(!is.null(caption)) {
      .self$add.p(caption)
    }
    return(.self$add.p("Empty dataset"))
  }
  y = data.frame(x=x)
  if(!is.null(name)) {
    names(y) <- c(name)
  }
  if(show.row.names) {
    row.names(y) <- names(x)
  }
  return(.self$add.table(y,caption=caption,show.row.names=show.row.names,...))
})

PandocAT$methods(add.p = function(x,rule=F,echo=T,...) {
  if(rule) {
    .self$add.paragraph(pandoc.horizontal.rule.return())
  }
  #if(!is.null(.self$tag) && length(.self$tag)!=0 && .self$tag != "") {
  #  x = paste("tag",.self$tag,x,collapse=":")
  #}
  if(echo) {
    cat("tag:",.self$tag,x,"\n")
  }
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
    .self$add.p(caption)
  }
  return(.self$add.p(pandoc.as.printed.return(x,...),echo=echo))
})

PandocAT$methods(add.header = function(x,echo=T,...) {
  return(.self$add.p(pandoc.header.return(x,...),echo=echo))
})

PandocAT$methods(set.tag = function(tag,echo=T,...) {
  .self$tag = tag
})

PandocAT$methods(save = function(f) {
  
  if (missing(f))
    f <- tempfile('pander-', getwd())
  fp    <- sprintf('%s.md', f)
  timer <- proc.time()
  
  ## create pandoc file
  cat(pandoc.title.return(.self$title, .self$author, .self$date), file = fp)
  lapply(.self$body, function(x) cat(paste(pandoc.return(x$result), collapse = '\n'), file = fp, append = TRUE))
    
})


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


