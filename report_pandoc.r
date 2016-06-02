library(ggplot2)
#library(knitr)
library(pander)

panderOptions("round",4)
panderOptions("table.style","rmarkdown") #"grid"
panderOptions("table.split.table",Inf)
panderOptions("table.alignment.default","left")
panderOptions("evals.messages",F)
#panderOptions("graph.fontsize",10)
evalsOptions("cache",F)
evalsOptions("cache.mode","environment")
#evalsOptions("output",c("result"))
evalsOptions("output",c("all"))
evalsOptions("graph.output","svg") #"png"
evalsOptions("graph.unify",F)
evalsOptions("width",800)
evalsOptions("height",640)
evalsOptions("res",75)
evalsOptions("hi.res",T)
evalsOptions("hi.res.width",1200)
evalsOptions("graph.env",F)
evalsOptions("graph.recordplot",F)
evalsOptions("graph.RDS",F)
##using graph.name option causes seemingly unconnected
##errors starting with warnings like:
##`No pander method for "ggplot", reverting to default`
##described in the URL but still not solved apparently
#https://github.com/Rapporter/rapport/issues/98
#evalsOptions("graph.name","plot-%d-%n")
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

are.automatic.rownames <- function(df) {
  all(rownames(df) == paste(seq(nrow(df))))
}

tempfile.unix <- function(...) {
  x = tempfile(...)
  gsub("\\","/",x,fixed=T)
}

html.path <- function(...) {
  paste(...,sep="/")
}

#' Write citations for a vector of package names into file in BibTex format
#' TODO can just use bibtex::write.bib
#' TODO after writing (or before if BibTex allows that type), optionally replace 
#' Manual type with TechReport that Zotero understands in BibTex (converts to Report). 
#' Otherwise Zotero
#' replaces Manual with Book. The replacement parameter should be a list of
#' to:from tuples. In the path BibTex -> Zotero -> RIS -> Endnote Web Page
#' gets converted to Journal Article still. It seems that in certain styles in
#' EndNote (ACS), the only way too show URL is to set type to Web Page. Otherwise
#' it is not clear at all that packages are CRAN packages.
citation.to.file <- function(package,file.name,append=F,...) {
  cit = unlist(sapply(package,function(p) toBibtex(citation(p,...))))
  write(cit,file.name,append=append)
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

new_section_path_el <- function(num=0,sub=F,has.header=F) {
  list(num=num,sub=sub,has.header=has.header)
}

extract.path.nums.section.path <- function(x) {
  sapply(x,function(y) y$num)
}

extract.path.nums.report.section <- function(x) {
  extract.path.nums.section.path(x$path)
}

extract.path.subs.section.path <- function(x) {
  sapply(x,function(y) y$sub)
}

extract.path.subs.report.section <- function(x) {
  extract.path.subs.section.path(x$path)
}

get.sub.level.section.path <- function(x) {
  sum(extract.path.subs.section.path(x))
}

cut.to.bottom.sub.section.path <- function(x) {
  subs = extract.path.subs.section.path(x)
  pos = length(subs) - match(T,rev(subs)) + 1
  if(is.na(pos)) {
    return (list())
  }
  else {
    if(pos>1)
      return (x[1:pos-1])
    else
      return (list())
  }
}

get.default.section<-function() {
  x = new.env(parent=emptyenv())
  x$path = list(new_section_path_el())
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

## when x is NULL, this function creates a clone
push.report.section<-function(x=NULL,sub=F,has.header=F) {
  if(is.null(x)) {
    x = get.report.section(default=NULL)
    if(is.null(x)) {
      x = get.default.section()
    }
    else {
      x = clone.report.section(x)
    }
  }
  x$path[[length(x$path)+1]] = new_section_path_el(num=0,sub=sub,has.header=has.header)
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
  if(last>0) {x$path[[last]]$num = x$path[[last]]$num + 1}
  else {x$path[[1]]$num = 1}
  return (x)
}

incr.report.section.if.zero <- function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()
  }
  last = length(x$path)
  if(last>0 && x$path[[last]]$num == 0) {x$path[[last]]$num = x$path[[last]]$num + 1}
  return (x)
}

format.section.path<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()$path
  }
  num = extract.path.nums.section.path(x)
  return (paste("\\(",paste(num,sep="",collapse="."),"\\)",sep=""))
}

format.report.section<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()
  }
  format.section.path(x$path)
}

format.section.path.as.file<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()$path
  }
  num = extract.path.nums.section.path(x)
  return (paste(num,sep="",collapse="."))
}

format.report.section.as.file<-function(x=NULL) {
  if(is.null(x)) {
    x = get.report.section()
  }
  format.section.path.as.file(x$path)
}

pandoc.as.printed.return <- function(x,attrs="") {
  x = capture.output(print(x))
  paste0('\n', repChar('`', 7), 
         ifelse(attrs == '', '', sprintf('{%s}', attrs)), 
         '\n', paste(x, 
                     collapse = '\n'), 
         '\n', repChar('`', 7), '\n')
}

pandoc.special.symb = "\\-\\[\\]`*_{}()#+!~"

pandoc.escape.special <- function(x) {
  x = gsub(paste('([',pandoc.special.symb,'])',sep=''),"\\\\\\1",
           format(x,digits=panderOptions("digits")),perl=T)
  ##the only way to have backstick is table cells is to use special symbol
  gsub('[|]','&#124;',x)
}

pandoc.link.verbatim.return <- function(url,text=NULL) {
  if(is.null(text)) {
    text = url
  }
  pandoc.link.return(url,pandoc.verbatim.return(text))
}

pandoc.anchor.return <- function(anchor,text) {
  anchor.tag = sprintf('<a name="%s"></a>',anchor)
  ret = sprintf("%s%s",
                anchor.tag,
                pandoc.link.verbatim.return(sprintf("#%s",anchor),
                                            text))
  return(ret)
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
                          'out.file' = 'character',
                          'out.formats' = 'character',
                          'portable.html' = 'logical',
                          'object.index' = 'list',
                          'data.dir' = 'character',
                          'widget.dir' = 'character',
                          'widget.deps.dir' = 'character'
                        )
)

PandocAT$methods(initialize = function(
  author = "Anonymous",
  title = "Analysis",
  out.file = "report",
  out.formats = c("html"),
  incremental.save = F,
  portable.html=T,
  ...
) {
  #.self$author=author
  #.self$title=title
  .self$out.file=out.file
  .self$out.formats=out.formats
  .self$incremental.save=incremental.save
  .self$portable.html=portable.html
  .self$object.index=list(table=1,figure=1)
  .self$data.dir = "data"
  unlink(.self$data.dir,recursive=T,force=T)
  dir.create(.self$data.dir, showWarnings = FALSE, recursive = TRUE)
  
  graph.dir = evalsOptions("graph.dir")
  unlink(graph.dir,recursive=T,force=T)
  dir.create(graph.dir, showWarnings = FALSE, recursive = TRUE)
  
  .self$widget.dir = "." #normalizePath(graph.dir,winslash = "/")
  .self$widget.deps.dir = "widget_deps" #html.path(.self$widget.dir,"widget_deps")
  dir.create(.self$widget.deps.dir, showWarnings = FALSE, recursive = TRUE)
  
  callSuper(author=author,title=title,...)
  
})

## private service method - should be called whenever an element is
## appended to the .self$body
PandocAT$methods(priv.append.section = function() {
  incr.report.section.if.zero()
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

PandocAT$methods(format.caption = function(caption,section=NULL,type=NULL) {
  if(!is.null(caption)) {  
    if(is.null(type)) {
      type = ""
      ind = ""
    }
    else {
      ind = .self$priv.append.index(type)
    }
    if(nzchar(ind)) {
      anchor = sprintf('%s.%s',type,ind)
      anchor.name = sprintf("%s %s.",Hmisc::capitalize(type),ind)
      name = pandoc.anchor.return(anchor,anchor.name)
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


PandocAT$methods(add.widget = function(x,new.paragraph=T,
                                       caption=NULL,
                                       show.image.links=T,
                                       width = 800,
                                       height = 800,
                                       data.export = NULL,
                                       data.export.descr = NULL,
                                       show.inline = T,
                                       ...) {
  
  require(htmlwidgets)
  
  if(new.paragraph) {
    .self$add.p("")
  }
  
  name.base=paste(str.to.file.name(caption,20),".html",sep="")
  
  fn = .self$make.file.name(name.base,dir=.self$widget.dir,make.unique=T)
  
  saveWidget(x,fn,selfcontained = F,libdir=.self$widget.deps.dir)
  
  caption.res = sprintf("Click to see HTML widget file in full window: %s",
                        pandoc.link.verbatim.return(fn))
  
  caption.type = "widget"
  
  if(!is.null(caption)) {
    caption = .self$format.caption(caption,type=caption.type)
  }
  
  caption = paste(caption,caption.res)
  
  if(!is.null(data.export)) {
    caption = sprintf("%s. Dataset is also saved here: %s", 
                      caption,
                      pandoc.link.verbatim.return(data.export))
    if(!is.null(data.export.descr)) {
      caption = sprintf("%s, %s",caption,data.export.descr)
    }
  }
  
  if(!is.null(caption)) {
    .self$add.p(caption)
  }
  if(show.inline) {
    if(is.null(width)) {
      width = evalsOptions("width")
    }
    
    if(is.null(height)) {
      height = evalsOptions("height")
    }
    iframe.tpl = '<iframe style="max-width=100%" 
        src="fn" 
        sandbox="allow-same-origin allow-scripts" 
        width="100%" 
        height="%s" 
        scrolling="no" 
        seamless="seamless" 
        frameBorder="0"></iframe>'
    iframe.tpl = '<iframe src="%s" width="%s" height="%s"> </iframe>'
    report$add(sprintf(iframe.tpl,
                       fn,
                       width,
                       height))
  }
  return(x)
})


PandocAT$methods(add = function(x,new.paragraph=T,
                                caption=NULL,
                                show.image.links=T,
                                caption.type=NULL,
                                graph.output = pander::evalsOptions("graph.output"),
                                hi.res = pander::evalsOptions("hi.res"),
                                ...) {
  
  timer           <- proc.time()
  par_fr <- parent.frame()
  if(!identical(.GlobalEnv, par_fr)) {
    env             <-  list2env(as.list(par_fr, all.names=TRUE),parent=parent.frame(2))
  }
  else {
    env = NULL
  }
  ## work around pander bug in v.0.6.0 where hi.res is created as a broken symlink plots/normal.res
  ## instead of just normal.res
  if(graph.output == 'svg') {
    hi.res = F
  }
  res             <- pander::evals(deparse(match.call()[[2]]),env=env,
                                   graph.output = graph.output,
                                   hi.res = hi.res,
                                   ...)
  if(new.paragraph) {
    .self$add.p("")
  }
  
  is.image = F
  caption.res = ""
  for (r in res) {
    if(any(r$type=="image")) {
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
    caption = .self$format.caption(caption,type=caption.type)
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

PandocAT$methods(reset.section = function(...) {
  return(NULL)
})

PandocAT$methods(default.section = function(...) {
  return(get.default.section(...))
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

PandocAT$methods(add.file = function(x,
                                      caption=NULL,
                                      wrap.caption=T,
                                      skip.if.empty=F,
                                      ...) {
  if (wrap.caption && !is.null(caption)) {
    caption = pandoc.escape.special(caption)
  }
  
  caption = .self$format.caption(caption,type="dataset")
  
  if(is.null(x)) {
    if(!skip.if.empty) {
      if(!is.null(caption)) {
        .self$add.p(caption)
      }
      return(.self$add.p("Empty dataset"))
    }
    else {
      return(.self)
    }
  }
    caption = paste(caption,
                    "Dataset is saved in a file (click to download)",
                    pandoc.link.verbatim.return(x)
    )

  return(.self$add.p(caption))
})


PandocAT$methods(add.table = function(x,
                                      show.row.names=is.matrix(x),
                                      echo=T,
                                      caption=NULL,
                                      wrap.vals=T,
                                      wrap.caption=T,
                                      split.tables=Inf,
                                      style="rmarkdown",
                                      export.to.file=T,
                                      show.first.rows=200,
                                      show.first.cols=200,
                                      skip.if.empty=F,
                                      ...) {
  if (wrap.caption && !is.null(caption)) {
    caption = pandoc.escape.special(caption)
  }
  
  caption = .self$format.caption(caption,type="table")
  
  if(is.null(x) || nrow(x)==0) {
    if(!skip.if.empty) {
      if(!is.null(caption)) {
        .self$add.p(caption)
      }
      return(.self$add.p("Empty dataset"))
    }
    else {
      return(.self)
    }
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
                    "Full dataset is also saved in a delimited text file (click to download and open e.g. in Excel)",
                    pandoc.link.return(file.name,pandoc.verbatim.return(file.name))
    )
  }
  
  if(show.first.rows > 0) {
    if(inherits(x,"data.table")) x = x[1:show.first.rows]
    else x = x[1:show.first.rows,,drop=F]
  }
  if(show.first.cols > 0) {
    if(inherits(x,"data.table")) x = x[,1:show.first.cols,with=F]
    else x = x[,1:show.first.cols,drop=F]
  }
  
  ## With data.table, I am getting this message:
  ## `data.table inherits from data.frame (from v1.5) but this data.table does not`
  ## when calling `rn = rownames()` below. Converting to data.frame here to get rid of it.
  
  if(inherits(x,"data.table")) x = as.data.frame(x)
  
  if(!show.row.names) {
    rownames(x) <- NULL
  }
  
  if(wrap.vals) {
    rn = rownames(x)
    if(show.row.names && !are.automatic.rownames(x)) {
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
    colnames(x) = pandoc.escape.special(colnames(x))
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
      .self$add.p(.self$format.caption(caption))
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
  .self$add.p(x,...)
})

PandocAT$methods(add.package.citation = function(x,...) {
  .self$add.p(capture.output(print(citation(x),style="text")))
})

PandocAT$methods(add.printed = function(x,caption=NULL,echo=T,...) {
  if(!is.null(caption)) {
    .self$add.p(.self$format.caption(caption))
  }
  return(.self$add.p(pandoc.as.printed.return(x,...),echo=echo))
})

PandocAT$methods(add.header = function(x,level=NULL,report.section=NULL,section.action="incr",echo=T,sub=F,...) {
  if(sub) {
    section.action = "push"
  }
  do.clone = is.null(report.section)
  report.section = switch(section.action,
                          incr=incr.report.section(report.section),
                          push=incr.report.section(report.section),
                          keep=get.report.section(report.section))
  num = extract.path.nums.report.section(report.section)
  if (is.null(level)) {
    ##headers will shift to the left above level 5 in HTML output
    level = min(5,length(num))
  }
  x = paste(format.report.section(report.section),x)
  ##newlines currently break header formatting, remove them
  x = gsub("\n"," ",x)
  .self$add.p(pandoc.header.return(x,level=level,...),echo=echo)
  
  if(section.action=="push") {
    rep.sec.push = NULL
    if(!do.clone) rep.sec.push = report.section
    #w/o argument it clones
    report.section = push.report.section(rep.sec.push,sub=sub,has.header=T)
  }
  
  if(.self$incremental.save) {
    .self$save()
  }
  
  return (report.section)
  
})

PandocAT$methods(make.file.name = function(name.base="",
                                           make.unique=T,
                                           dir=NULL,
                                           section.path=NULL) {
  if(is.null(dir)) {
    dir = .self$data.dir
  }
  if(length(name.base)==0) {
    name.base = ""
  }
  if(name.base=="" && !make.unique) {
    stop("Need either non-empty name.base or make.unique=T")
  }
  fn.start = format.section.path.as.file(section.path)
  if(make.unique) {
    fn = tempfile.unix(paste(fn.start,"-",sep=""),tmpdir=dir,fileext=name.base)
  }
  else {
    fn = file.path(dir,paste(fn.start,name.base,sep="-"),fsep="/")
  }
  return(fn)
})

## Save data as a delimited text file
## ... are optional arguments to write.table
PandocAT$methods(write.table.file = function(data,
                                             name.base,
                                             make.unique=T,
                                             descr=NULL,
                                             row.names=F,
                                             row.names.header=T,
                                             ...) {
  ## if we write row.names, Excel shifts header row to the left when loading
  if(row.names && row.names.header) {
    data = cbind(rownames=rownames(data),data)
    row.names=F
  }
  fn = .self$make.file.name(name.base,make.unique=make.unique)
  write.table(data,
              fn,
              sep="\t",
              row.names = row.names,
              ...)
  if (!is.null(descr)) {
    .self$add.descr(paste("Wrote",descr,"to file",
                          pandoc.link.return(fn,fn)))
  }
  return(fn)
})

PandocAT$methods(save = function(out.file.loc,out.formats.loc,portable.html.loc,sort.by.sections=F) {
  
  if (missing(out.file.loc)) {
    out.file.loc = .self$out.file
    if(missing(out.file.loc)) {
      out.file.loc = "report"
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
  
  fp    <- out.file.loc
  timer <- proc.time()
  
  fp.all = list()
  
  f_sections = .self$sections
  f_body = .self$body
  
  if(sort.by.sections) {
    ##sort by section lexicographically, using a stable sort
    sect_ord = sort.list(
      unlist(lapply(.self$sections,format.section.path)),
      method="shell")
    f_body = f_body[sect_ord]
  }
  
  write.el <- function(el,fp) {
    el.str = pander_return(el$result)
    cat(paste(el.str, collapse = '\n'), 
        file = fp, append = TRUE)    
  }
  
  for(i.el in seq_along(f_body)) {
    section = f_sections[[i.el]]
    el = f_body[[i.el]]
    
    section.par = cut.to.bottom.sub.section.path(section)
    
    #print(paste("Full section:",paste(section,collapse=" ")))
    #print(paste("Par section:",paste(section.par,collapse=" ")))
    
    if(length(section.par) > 0) {
      sub.path = section.par
    }
    else {
      sub.path = NULL
    }
    fp.sub = make.file.name(name.base=fp,
                            make.unique=F,
                            dir=".",
                            section.path=sub.path)
    fp.sub.md = paste(fp.sub,".Rmd",sep="") #".md"
    #print(paste("fp.sub=",fp.sub))
    
    if(is.null(fp.all[[fp.sub.md]])) {
      cat(pandoc.title.return(.self$title, .self$author, .self$date), file = fp.sub.md)
    }    
    
    if(i.el>1) {
      sub.level.prev = get.sub.level.section.path(f_sections[[i.el-1]])
      sub.level = get.sub.level.section.path(section)
      if(sub.level > sub.level.prev) {
        cat(pandoc.link.verbatim.return(paste(fp.sub,".html",sep=""),"Subreport"), #".html"
            file = fp.sub.md.prev, append = TRUE)
        if(section[[length(section)]]$has.header) {
          write.el(f_body[[i.el-1]],fp.sub.md)
        }
      }
    }    
    
    write.el(el,fp.sub.md)
    fp.all[[fp.sub.md]] = 1
    fp.sub.md.prev = fp.sub.md
  }
  
  for(out.format in out.formats.loc) {
    for(fp.sub.md in names(fp.all)) {
      ## It would be nice to add `options="-s -S"` to support
      ## Pandoc's subscript and suprscript extensions, but
      ## this will entirely replace internal default options and
      ## break TOC etc
      cat(sprintf("Pandoc converting markdown file %s to %s format\n",fp.sub.md,out.format))
      Pandoc.convert(fp.sub.md,format=out.format,open=F,footer=F,
                     portable.html=portable.html.loc)
    }
  }
  
})


tmp.save <- function(report,out.file.loc,out.formats.loc,portable.html.loc,sort.by.sections=F) {
  .self = report$.self
  
  if (missing(out.file.loc)) {
    out.file.loc = .self$out.file
    if(missing(out.file.loc)) {
      out.file.loc = "report"
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
  
  fp    <- out.file.loc
  timer <- proc.time()
  
  fp.all = list()
  
  f_sections = .self$sections
  f_body = .self$body
  
  if(sort.by.sections) {
    ##sort by section lexicographically, using a stable sort
    sect_ord = sort.list(
      unlist(lapply(.self$sections,format.section.path)),
      method="shell")
    f_body = f_body[sect_ord]
  }
  
  write.el <- function(el,fp) {
    el.str = pander_return(el$result)
    el.str = gsub("_[","_\\[",el.str,fixed=T)
    cat(paste(el.str, collapse = '\n'), 
        file = fp, append = TRUE)    
  }
  
  for(i.el in seq_along(f_body)) {
    section = f_sections[[i.el]]
    el = f_body[[i.el]]
    
    section.par = cut.to.bottom.sub.section.path(section)
    
    #print(paste("Full section:",paste(section,collapse=" ")))
    #print(paste("Par section:",paste(section.par,collapse=" ")))
    
    if(length(section.par) > 0) {
      sub.path = section.par
    }
    else {
      sub.path = NULL
    }
    fp.sub = report$make.file.name(name.base=fp,
                                   make.unique=F,
                                   dir=".",
                                   section.path=sub.path)
    fp.sub.md = paste(fp.sub,".md",sep="")
    #print(paste("fp.sub=",fp.sub))
    
    if(is.null(fp.all[[fp.sub.md]])) {
      cat(pandoc.title.return(.self$title, .self$author, .self$date), file = fp.sub.md)
    }    
    
    if(i.el>1) {
      sub.level.prev = get.sub.level.section.path(f_sections[[i.el-1]])
      sub.level = get.sub.level.section.path(section)
      if(sub.level > sub.level.prev) {
        cat(pandoc.link.verbatim.return(paste(fp.sub,".html",sep=""),"Subreport"), 
            file = fp.sub.md.prev, append = TRUE)
        if(section[[length(section)]]$has.header) {
          write.el(f_body[[i.el-1]],fp.sub.md)
        }
      }
    }    
    
    write.el(el,fp.sub.md)
    fp.all[[fp.sub.md]] = 1
    fp.sub.md.prev = fp.sub.md
  }
  
  for(out.format in out.formats.loc) {
    for(fp.sub.md in names(fp.all)) {
      ## It would be nice to add `options="-s -S"` to support
      ## Pandoc's subscript and suprscript extensions, but
      ## this will entirely replace internal default options and
      ## break TOC etc
      cat(sprintf("Pandoc converting markdown file %s to %s format\n",fp.sub.md,out.format))
      Pandoc.convert(fp.sub.md,format=out.format,open=F,footer=F,
                     portable.html=portable.html.loc)
    }
  }
  
}


test_report.sections<-function() {
  
  report <- PandocAT$new("noone@email.com","developer")
  
  get.pandoc.section.3<-function() {
    print("get.pandoc.section.3")
    #print(incr.report.section())
    report.section = report$add.header("get.pandoc.section.3",section.action="push",sub=T)
    report$add.descr("Plots are shown with relation to various combinations of meta 
                   data variables and in different graphical representations. Lots of plots here.")
    
    #report$add.header("Iterating over all combinations of grouping variables")
    
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


