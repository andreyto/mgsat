library(ggplot2)
#library(knitr)
library(pander)

panderOptions("round",4)
panderOptions("table.style","rmarkdown") #"grid"
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
## This is ggplot2 function that changes
## text base size within a current theme for entire session
#theme_set(theme_gray(base_size = 20))

make.global <- function(var) {
  assign(deparse(substitute(var)),var,envir=globalenv()) 
}

# print a named list as a string of named function arguments
arg.list.as.str<-function(x,collapse=",") {
  paste("[",
        paste(capture.output(str(x,no.list=T,comp.str="",give.attr=F,give.head=F)),collapse=collapse),
        "]",
        sep=""
  )
}

#' Adopted from phyloseq code
#' Computes text size of axis label based on the number of
#' labels.
#' Maybe R strwidth can be used even with ggplot2?
calc.text.size <- function(n, mins=0.5, maxs=4, B=6, D=100){
  # empirically selected size-value calculator.
  s <- B * exp(-n/D)
  # enforce a floor.
  s <- ifelse(s > mins, s, mins)
  # enforce a max
  s <- ifelse(s < maxs, s, maxs)
  return(s)
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

format.report.section.as.file<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()
  }
  return (paste(x$path,sep="",collapse="."))
}


pandoc.as.printed.return <- function(x,attrs="") {
  x = capture.output(print(x))
  paste0('\n', repChar('`', 7), 
         ifelse(attrs == '', '', sprintf('{%s}', attrs)), 
         '\n', paste(x, 
                     collapse = '\n'), 
         '\n', repChar('`', 7), '\n')
}

pandoc.special.symb = "`*_{}()#+!~"

pandoc.escape.special <- function(x) {
  gsub(paste('([',pandoc.special.symb,'])',sep=''),"\\\\\\1",
       format(x,digits=panderOptions("digits")))
}

pandoc.link.verbatim.return <- function(url,text=NULL) {
  if(is.null(text)) {
    text = url
  }
  pandoc.link.return(url,pandoc.verbatim.return(text))
}

## make string x a (more or less) valid file name
str.to.file.name <- function(x,max.length=0) {
  x = gsub('[^-[:alnum:]._]+','.',x)
  if(max.length>0) {
    x = substring(x,1,max.length)
  }
  x
}

PandocAT <- setRefClass('PandocAT', contains = "Pandoc", 
                        fields = list(
                          'sections' = 'list',
                          'incremental.save' = 'logical',
                          'out.file.md' = 'character',
                          'out.formats' = 'character',
                          'portable.html' = 'logical',
                          'object.index' = 'list',
                          'data.dir' = 'character'
                        )
)

PandocAT$methods(initialize = function(
  author = "Anonymous",
  title = "Analysis",
  out.file.md = "report.md",
  out.formats = c("html"),
  incremental.save = F,
  portable.html=T,
  ...
) {
  #.self$author=author
  #.self$title=title
  .self$out.file.md=out.file.md
  .self$out.formats=out.formats
  .self$incremental.save=incremental.save
  .self$portable.html=portable.html
  .self$object.index=list(table=1,figure=1)
  .self$data.dir = "output"
  unlink(.self$data.dir,recursive=T,force=T)
  dir.create(.self$data.dir, showWarnings = FALSE, recursive = TRUE)
  
  graph.dir = evalsOptions("graph.dir")
  unlink(graph.dir,recursive=T,force=T)
  dir.create(graph.dir, showWarnings = FALSE, recursive = TRUE)
  
  callSuper(author=author,title=title,...)
})

## private service method - should be called whenever an element is
## appended to the .self$body
PandocAT$methods(priv.append.section = function() {
  .self$sections[[length(.self$sections)+1]] = get.report.section()$path
})

## private service method - should be called whenever an element 
## that needs its own index in the report is
## appended to the .self$body
PandocAT$methods(priv.append.index = function(type) {
  val = .self$object.index[[type]]
  if(is.null(val)) {
    val = .self$object.index[[type]] = 1
  }
  .self$object.index[[type]] = val + 1
  return(val)
})

PandocAT$methods(priv.format.caption = function(caption,section=NULL,type=NULL) {
  if(!is.null(caption)) {  
    if(is.null(type)) {
      type = ""
      ind = ""
    }
    else {
      ind = .self$priv.append.index(type)
    }
    if(nzchar(ind)) {
      anchor.name = sprintf('%s.%s',type,ind)
      anchor = sprintf('<a name="%s"></a>',anchor.name)
      name = sprintf("%s%s",
                     anchor,
                     pandoc.link.return(sprintf("#%s",anchor.name),
                                        sprintf("%s %s.",capitalize(type),ind)))
    }
    else {
      name = ""
    }
    caption = paste(format.report.section(section),name,caption)
    if(substr(caption,nchar(caption),nchar(caption)+1)!=".") {
      caption = paste(caption,".",sep="")
    }
  }
  
  return (caption)
})

PandocAT$methods(add = function(x,new.paragraph=T,
                                caption=NULL,
                                show.image.links=T,
                                caption.type=NULL,
                                ...) {
  
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
  
  is.image = F
  caption.res = ""
  for (r in res) {
    if(r$type=="image") {
      if(show.image.links) {
        rr = r$result
        caption.res = paste(caption.res,
                            sprintf("Image file: %s.",
                                    pandoc.link.verbatim.return(as.character(rr)))
        )
        hres.ref = attr(rr,"href")
        if(!is.null(hres.ref)) {
          caption.res = paste(caption.res,
                              sprintf("High resolution image file: %s.",
                                      pandoc.link.verbatim.return(hres.ref))
          )
        }
      }
      is.image = T
    }
  }
  if(is.null(caption.type)) {
    if(is.image) {
      caption.type = "figure"
    }
  }
  if(!is.null(caption)) {
    caption = .self$priv.format.caption(caption,type=caption.type)
  }
  if(nzchar(caption.res)) {
    caption = paste(caption,caption.res)
  }
  if(!is.null(caption)) {
    .self$add.p(caption)
  }
  .self$body      <- c(.self$body,res)
  .self$priv.append.section()
  .self$proc.time <- .self$proc.time + as.numeric(proc.time() - timer)[3]
  return(res)
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

PandocAT$methods(add.table = function(x,
                                      show.row.names=F,
                                      echo=T,
                                      caption=NULL,
                                      wrap.vals=T,
                                      wrap.caption=T,
                                      split.tables=180,
                                      style="rmarkdown",
                                      export.to.file=T,
                                      show.first.rows=200,
                                      show.first.cols=200,
                                      ...) {
  if (wrap.caption && !is.null(caption)) {
    caption = pandoc.escape.special(caption)
  }
  
  caption = .self$priv.format.caption(caption,type="table")
  
  if(is.null(x) || nrow(x)==0) {
    if(!is.null(caption)) {
      .self$add.p(caption)
    }
    return(.self$add.p("Empty dataset"))
  }
  if(show.first.rows > 0) {
    if(show.first.rows >= nrow(x)) {
      show.first.rows = 0
    }
  }
  if(show.first.rows > 0) {
    caption = paste(caption,sprintf("Showing only %s first rows.",show.first.rows))
  }
  if(show.first.cols > 0) {
    if(show.first.cols >= ncol(x)) {
      show.first.cols = 0
    }
  }
  if(show.first.cols > 0) {
    caption = paste(caption,sprintf("Showing only %s first columns.",show.first.cols))
  }
  
  if(export.to.file) {
    file.name = .self$write.table.file(x,
                                       name.base=paste(str.to.file.name(caption,20),".tsv",sep=""),
                                       descr=NULL,
                                       row.names=show.row.names,
                                       row.names.header=T)    
    caption = paste(caption,
                    "Full dataset is also saved in a delimited text file",
                    pandoc.link.return(file.name,pandoc.verbatim.return(file.name))
    )
  }
  
  if(show.first.rows > 0) {
    x = x[1:show.first.rows,,drop=F]
  }
  if(show.first.cols > 0) {
    x = x[,1:show.first.cols,drop=F]
  }
  
  if(!show.row.names) {
    rownames(x) <- NULL
  }
  
  if(wrap.vals) {
    rn = rownames(x)
    if(!is.null(rn)) {
      rn = pandoc.escape.special(rn)
    }
    if(is.matrix(x)) {
      x = as.data.frame(x)
    }
    x = sapply(x,pandoc.escape.special,USE.NAMES=F,simplify=T)
    if(!is.matrix(x)) {
      x = t(as.matrix(x))
    }
    rownames(x) = rn
  }
  
  .self$add.p(caption)
  tbl_p <- pandoc.table.return(x,split.tables=split.tables,style=style,caption=NULL,...)
  if(echo) {
    print(tbl_p)
  }
  return(.self$add.p(tbl_p))
})

PandocAT$methods(add.vector = function(x,name=NULL,
                                       show.row.names=T,
                                       caption=NULL,
                                       ...) {
  if(is.null(x) || length(x)==0) {
    if(!is.null(caption)) {
      .self$add.p(.self$priv.format.caption(caption))
    }
    return(.self$add.p("Empty dataset"))
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
    ##headers will shift to the left above level 5 in HTML output
    level = min(5,length(report.section$path))
  }
  x = paste(format.report.section(report.section),x)
  .self$add.p(pandoc.header.return(x,level=level,...),echo=echo)
  
  if(section.action=="push") {
    #w/o argument it clones
    report.section = push.report.section()
  }
  
  if(.self$incremental.save) {
    .self$save()
  }
  
  return (report.section)
  
})

PandocAT$methods(make.file.name = function(name.base="",make.unique=T) {
  if(length(name.base)==0) {
    name.base = ""
  }
  if(name.base=="" && !make.unique) {
    stop("Need either non-empty name.base or make.unique=T")
  }
  fn.start = format.report.section.as.file()
  if(make.unique) {
    fn = tempfile(paste(fn.start,"-",sep=""),tmpdir=.self$data.dir,fileext=name.base)
  }
  else {
    fn = file.path(.self$data.dir,paste(fn.start,name.base,sep="-"))
  }
  return(fn)
})

PandocAT$methods(write.table.file = function(data,
                                             name.base,
                                             make.unique=T,
                                             descr=NULL,
                                             row.names=F,
                                             row.names.header=T) {
  ## if we write row.names, Excel shifts header row to the left when loading
  if(row.names && row.names.header) {
    data = cbind(rownames=rownames(data),data)
    row.names=F
  }
  fn = .self$make.file.name(name.base,make.unique=make.unique)
  write.table(data,
              fn,
              sep="\t",
              row.names = row.names)
  if (!is.null(descr)) {
    .self$add.descr(paste("Wrote",descr,"to file",
                          pandoc.link.return(fn,fn)))
  }
  return(fn)
})

PandocAT$methods(save = function(out.file.md.loc,out.formats.loc,portable.html.loc,sort.by.sections=F) {
  
  if (missing(out.file.md.loc)) {
    out.file.md.loc = .self$out.file.md
    if(missing(out.file.md.loc)) {
      out.file.md.loc = "report.md"
    }
  }
  
  if (missing(out.formats.loc)) {
    out.formats.loc = .self$out.formats
    if(missing(out.formats.loc)) {
      out.formats.loc = c("html")
    }
  }
  
  if (missing(portable.html.loc)) {
    portable.html.loc = .self$portable.html
    if(missing(portable.html.loc)) {
      portable.html.loc = T
    }
  }
  
  fp    <- out.file.md.loc
  timer <- proc.time()
  
  ## create pandoc file
  cat(pandoc.title.return(.self$title, .self$author, .self$date), file = fp)
  f_sections = .self$sections
  f_body = .self$body
  if(sort.by.sections) {
    ##sort by section lexicographically, using a stable sort
    sect_ord = sort.list(
      unlist(lapply(.self$sections,function(x) paste(x,sep="",collapse="."))),
      method="shell")
    f_body = f_body[sect_ord]
  }
  lapply(f_body, function(x) cat(paste(pander.return(x$result), collapse = '\n'), 
                                 file = fp, append = TRUE))
  
  for(out.format in out.formats.loc) {
    ## It would be nice to add `options="-s -S"` to support
    ## Pandoc's subscript and suprscript extensions, but
    ## this will entirely replace internal default options and
    ## break TOC etc
    Pandoc.convert(fp,format=out.format,open=F,footer=F,
                   portable.html=portable.html.loc)
  }
  
})

test_report.sections<-function() {
  
  report <- PandocAT$new("noone@email.com","developer")
  
  get.pandoc.section.3<-function() {
    print("get.pandoc.section.3")
    #print(incr.report.section())
    report.section = report$add.header("get.pandoc.section.3",section.action="push")
    report.section = report$add.header("get.pandoc.section.3.1")
    report$add.table(data.frame(A=c("a","b"),B=c(1,2)),caption="Table")
    report$add.descr(paste("File name with extra output is ",report$make.file.name("data.csv")))
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
  
  report$save("test_sections")
  Pandoc.convert("test_sections.md",format="html",open=F)
  
}

test_report.tables<-function() {
  
  report <- PandocAT$new("noone@email.com","developer")
  
  y = data.frame(a=c(1,2),b=c(2,4),c=c("zzz","mmmm"))
  
  report$add.table(y,caption="Test table 1",style="multiline")
  
  y = data.frame(a=c(1))
  
  report$add.table(y,caption="Test table 2",style="grid")
  
  x = c(x="_a_",y="b",z="c")
  
  report$add.vector(x,name="Clade",caption="Test vector 1",
                    show.row.names=T,wrap.vals=F,style="grid")
  report$add.vector(x,name="Clade",caption="Test vector 1",
                    show.row.names=T,wrap.vals=T,style="grid")
  report$add.vector(x,name="Clade",caption="Test vector 1",
                    show.row.names=F,wrap.vals=T,style="grid")
  
  report$save("test_tables")
  
  Pandoc.convert("test_tables.md",format="html",open=F)
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


