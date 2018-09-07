
set_trace_options<-function(try.debug=T) {
  #tell our custom tryCatchAndWarn not to catch anything
  options(try.debug=try.debug)
  #setup more verbose error reporting
  #set warn to 2 to convert warinings to errors; 0 - to keep going
  options(warn = 0, keep.source = TRUE, error = 
            quote({ 
              cat("Environment:\n", file=stderr()); 
              
              # TODO: setup option for dumping to a file (?)
              # Set `to.file` argument to write this to a file for post-mortem debugging    
              dump.frames();  # writes to last.dump
              
              #
              # Debugging in R
              #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/index.shtml
              #
              # Post-mortem debugging
              #   http://www.stats.uwo.ca/faculty/murdoch/software/debuggingR/pmd.shtml
              #
              # Relation functions:
              #   dump.frames
              #   recover
              # >>limitedLabels  (formatting of the dump with source/line numbers)
              #   sys.frame (and associated)
              #   traceback
              #   geterrmessage
              #
              # Output based on the debugger function definition.
              
              n <- length(last.dump)
              calls <- names(last.dump)
              cat(paste("  ", 1L:n, ": ", calls, sep = ""), sep = "\n", file=stderr())
              cat("\n", file=stderr())
              
              if (!interactive()) {
                q()
              }
            }))
  #options(error=recover)
}



## bind variables in the current scope to the global environment
## can be used in debugging
## if NULL is passed, the parents (calling) environment 
## is cloned and assigned to a global variable passed in name
## except – list of names not to make global (in case there are already globals with these names that should be preserved).
## Example:
## Somewhere in the calling stack of your code, write:
## make.global()
## stop(“DEBUG”)
## Then, all local variables before the stop() will become global. Note: sometimes something slips inside R interpreter, and make.global() does nothing, 
## or takes a few seconds of idle shell to take effect. Make.global(name=“dbg”) is more robust. It creates new environment, and copies local variables into it.
## You access them like dbg$var_name
## Another use pattern:
## You have some local variable var_name, then:
## make.global(var_name)
## The above only makes var_name global. Note: these are copies as per usual R semantics, not reference.
## Beware to start new session for production (non-debugging runs), otherwise you can get tripped by the global variables that you created during debugging.
make.global <- function(var=NULL,name=NULL,except=c("report.section")) {
  if(is.null(var)) {
    if(is.null(name)) {
      name="global"
    }
    p.e = parent.frame()
    if(name=="global") {
      t.e = globalenv()
      except = c("...",except)
      for(n in ls(p.e, all.names=TRUE)) {
        if(! (n %in% except )) assign(n, get(n, p.e), t.e)
      }
      return()
    }
    else {
      var = as.environment(as.list(p.e, all.names=TRUE))
    }
  }
  if(is.null(name)) {
    name = deparse(substitute(var))
  }
  assign(name,var,envir=globalenv()) 
}


take_first<-function(x,n) {
  return (x[1:min(length(x),n)])
}

str_blank <- function(x) {
  return (nchar(gsub("\\s","",x))==0)
}

## Take named list and replace values with data frame columns
## if original values match data frame column names; otherwise
## keep original values
interpret.args.in.df <- function(args,data) {
  args = plyr::compact(args)
  names = colnames(data)
  args = lapply(args,function(x) if(length(x) && is.character(x) && (x %in% names)) data[,x] else x)
  args
}


## convert accented characters to regular ASCII equivalents
str_to_ascii <- function(x) iconv(x,to="ASCII//TRANSLIT")

str_dedent <- function(x,n_header=1,ignore_length_of_blank_lines=T) {
  lines_full = str_split(x,"\n")[[1]]
  n_lines = length(lines_full)
  n_header = min(n_header,n_lines)
  lines_header = lines_full[1:n_header]
  if(n_header >= n_lines) {
    return (x)
  }
  else {
    lines_full = lines_full[(n_header+1):n_lines]
    len_full = unlist(str_length(lines_full))
    lines_trimmed = str_trim(lines_full,"left")
    len_trimmed = unlist(str_length(lines_trimmed))
    len_diff = len_full - len_trimmed
    if(ignore_length_of_blank_lines) {
      len_diff = len_diff[len_trimmed>0]
    }
    len_dedent = min(len_diff)
    lines_padded = unlist(str_pad(lines_trimmed,len_full - len_dedent,side="left"))
    return (str_c(c(lines_header,lines_padded),collapse = "\n"))
  }
  #str_c(lines_padded,collapse = "\n")
  stop("Should not get here")
}

## Given character vector that has values as numbers with optional
## character suffix, create a factor with levels ordered numerically.
## Optionally, lef-pad and right-align the levels (and data values)
## so that lexicographic string sort would order them numerically.
## Returns either ordered (default) or unordered factor.
order.factor.by.numeric.prefix <- function(x,pad=NULL,ordered=T) {
  library(stringr)
  x = factor(x)
  xlev = levels(x)
  patt = "(.*[0-9])([^0-9]*$)"
  xlevpref = str_match(xlev,patt)[,-1]
  xlevnum = as.numeric(xlevpref[,1])
  xlevpref = xlevpref[order(xlevnum),]
  x_reformat <- function(x,pad=' ') {
    if(!is.null(pad)) {
      x = paste(format(as.numeric(x[,1]),
                       preserve.width="common",drop0trailing=T,trim=F),
                x[,2],
                sep="")
      if(pad!=' ') {
        xtr = str_trim(x,side = "left")
        x = str_pad(xtr,str_length(x),side="left",pad=pad)
      }
    }
    else {
      x = paste(x[,1],x[,2],sep="")
    }
    x
  }
  xlev = x_reformat(xlevpref,pad=pad)
  if(!is.null(pad)) {
    
    xpref = str_match(as.character(x),patt)[,-1]
    x = x_reformat(xpref,pad=pad)
  }
  factor(x,xlev,ordered = ordered)
}

file.dep.updated <- function(dep,targ) {
  dep = unlist(dep)
  targ = unlist(targ)
  if(!all(file.exists(targ))) {
    return (T)
  }
  if(!all(file.exists(dep))) {
    return (T)
  }
  for(dt in file.mtime(dep)) for(tt in file.mtime(targ)) {
    if(dt > tt) return (T)
  }
  return (F)
}

label.size.points <- function(x,what="width",resolution=72) {
  require(graphics)
  if(what=="width") {
    y = strwidth(x,units = "inches")
  }
  else if(what=="height") {
    y = strheight(x,units = "inches")
  }
  else {
    stop(sprintf("Uknown `what` parameter value: %s",what))
  }
  return (max(y)*resolution)
}

## replace oldnames with newnames in vector allnames
## to be used when allnames = names(data.frame)
replace.col.names<-function(allnames,oldnames,newnames,do.checks=T) {
  if (do.checks) {
    not.found = oldnames[!(oldnames %in% allnames)]
    if (length(not.found) > 0) {
      stop(sprintf("These old names were not found: %s",paste(not.found,collapse=", ")))
    }
    if(length(oldnames)!=length(newnames)) {
      stop(sprintf("There are %s old names and %s new names",lenght(oldnames),length(newnames)))
    }
  }
  allnames[match(oldnames,allnames)] = newnames
  return(allnames)
}

sub.col.names <- function(x,pattern,replacement,...) {
  cn = colnames(x)
  cn.new = sub(pattern,replacement,cn,...)
  colnames(x) = cn.new
  x
}

all_non_null <- function(...) {
  Filter(Negate(is.null), list(...))
}

first_non_null <- function(...) {
  x = all_non_null(...)
  if(length(x)>0) x[[1]]
  else NULL
}

## Update fields in named list x with fields in named list y
## The semantics as in Python dict.update()
update.list <- function(x,y) {
  x[names(y)] = y
  x
}

are.identical <- function(x,y) {
  if(is.character(x)) x = get(x)
  if(is.character(y)) y = get(y)
  identical(x,y)
}

## Return a copy of data frames or matrix with selected columns removed
drop.columns <- function(df,drop_names=c()) {
  mask_col_ignore = colnames(df) %in% drop_names
  return ( df[,!mask_col_ignore] )
}

## Return row names as factor with levels sorted in the original order of rows
rownames.as.factor <- function(x,ordered=F) {
  rn = rownames(x)
  factor(rn,levels=rn,ordered=ordered)
}

## Return column names as factor with levels sorted in the original order of columns
colnames.as.factor <- function(x,ordered=F) {
  rn = colnames(x)
  factor(rn,levels=rn,ordered=ordered)
}

## Wrapper for functions from various packages that recode values in a vector
recode.values <- function(x,...,flavor=c("car","DescTools")) {
  flavor = flavor[[1]]
  if(flavor == "car") {
    car::Recode(x,...)
  }
  else if(flavor == "DescTools") {
    DescTools::Recode(x,...)
  }
  else {
    stop(sprintf("Unknown flavor value: %s",flavor))
  }
}

#' Call quantcut and return its result as ordered factor.
#' As in gtools::quantcut(), q can be either the desired number of
#' quantiles, or a vector of quantiles (probabilities).
#' If q_are_breaks==T, the q is interpreted as a vector with breaks of length n_groups + 1,
#' and `base::cut` is called.
quantcut.ordered <- function(x,na.rm=T,q=4,q_are_breaks=F,...) {
  ##contrary to the docs, na.rm is ignored by quantcut causing stop, do it manually here
  if(na.rm) {
    x.nna = x[!is.na(x),drop=F]
    ##this line works around the error "'breaks' are not unique" when there are
    ##fewer unique values in x than the q. Note that this might get slow if x
    ##is a huge continuous vector
    if(length(q)==1) {
      q = min(q,length(unique(x.nna))+1)
    }
    if(q_are_breaks) {
      qq.nna = cut(x.nna,breaks = q,ordered_result = T,...)
    }
    else {
      qq.nna = gtools::quantcut(x.nna,q=q,...)
    }
    y = factor(rep(NA,length(x)),levels=c(levels(qq.nna),NA),ordered = T)
    y[!is.na(x)] = qq.nna
  }
  else {
    if(length(q)==1) {
      q = min(q,length(unique(x))+1)
    }
    if(q_are_breaks) {
      y = cut(x,breaks = q,ordered_result = T,...)
    }
    else {
      y = gtools::quantcut(x,q=q,...)
    }
  }
  as.ordered(y)
}

index_as_left_padded_str <- function(x,n_zeros=NULL) {
  if(is.null(n_zeros)) {
    n_zeros = ceiling(log10(length(x)))
  }
  frmt = paste0("%0",n_zeros,"i")
  sprintf(frmt,x)
}

#' Geometric mean of a 1-D vector
#' copied from https://stackoverflow.com/a/25555105
gm_mean = function(x, na.rm=TRUE, zero.propagate = FALSE){
  if(any(x < 0, na.rm = TRUE)){
    return(NaN)
  }
  if(zero.propagate){
    if(any(x == 0, na.rm = TRUE)){
      return(0)
    }
    exp(mean(log(x), na.rm = na.rm))
  } else {
    exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
  }
}

#' Geometric median (spatial median or L1-median) of a matrix (rows are observations)
gm_median = function(x,...) {
  Gmedian::Weiszfeld(x,...)$median
}

missing.join.keys.report <- function(x,y,by,do.report=T,name.x="x",name.y="y",report.only.keys=F) {
  x_y = dplyr::anti_join(x,y,by=by)
  if(report.only.keys) {
    if(!is.null(names(by))) {
      by_x = names(by)
    }
    else {
      by_x = by
    }
    x_y = dplyr::count_(x_y,by_x,sort=T)
  }
  if(do.report) {
    report$add.table(x_y,wrap.vals = T,skip.if.empty=T,
                     caption=sprintf("Keys present in '%s' but missing from '%s'",name.x,name.y))
  }
  return(x_y)
}

## Wrapper around dcast.data.table that allows to specify
## value.var="..." to mean "all variables not used in formula",
## reorder value columns by RHS of the formula and use RHS as
## a prefix of the output columnsrather than suffix as default.
dcast.ext.data.table <- function(data,formula,sep = "_", 
                                 value.var = "...",
                                 value.is.suffix = T,
                                 order.by.rhs = T,
                                 ...) {
  library(data.table)
  library(stringr)
  library(dplyr)
  if(any(value.var=="...")) {
    value.var = setdiff(colnames(data),all.vars(formula))
  }
  if(value.is.suffix || order.by.rhs) {
    use.sep = "@#%$^&*!" #something not likely to be in the value data
  }
  else {
    use.sep = sep
  }
  x = dcast.data.table(data = data,formula = formula,value.var = value.var, sep = use.sep,...)
  if(value.is.suffix || order.by.rhs) {
    cn = colnames(x)
    mask.val.col = grepl(use.sep,cn,fixed = T)
    cn.val.col = cn[mask.val.col]
    cn.other.col = cn[!mask.val.col]
    ##split by our unique separator and swap columns
    cn.val.col.spl = as.data.frame(str_split_fixed(cn.val.col,fixed(use.sep),2)[,c(2,1)])
    cn.val.col.spl$name = cn.val.col
    if(order.by.rhs) {
      cn.val.col.spl = dplyr::arrange_(cn.val.col.spl,colnames(tmp$cn.val.col.spl))
      setcolorder(x,c(cn.other.col,cn.val.col.spl$name))
    }
    if(!value.is.suffix){
      cn.val.col.spl = cn.val.col.spl[,c(2,1,3)]
    }
    cn.val.col.spl$new_name = paste(cn.val.col.spl[,1],cn.val.col.spl[,2],sep = sep)
    setnames(x,cn.val.col.spl$name,cn.val.col.spl$new_name)
  }
  x
}

## return x converted to class of `other`.
## This calls the constructor of other, which must exist,
## and passes x along with the remaining aurguments
call.ctor <- function(other,x,...) {
  do.call(class(other)[[1]],
          c(list(x),list(...)))
}

## copied from rMSA, GPL
mgsat.mafft <- function(x, param="--auto") {
  
  ## get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    file.remove(Sys.glob(paste(temp_file, ".*", sep="")))
    setwd(dir)
  })
  setwd(wd)
  
  infile <- paste(temp_file, ".in", sep="")
  outfile <- paste(temp_file, ".aln", sep="")
  reader <- if(is(x, "RNAStringSet")) Biostrings::readRNAMultipleAlignment
  else if(is(x, "DNAStringSet")) Biostrings::readDNAMultipleAlignment
  else if(is(x, "AAStringSet")) Biostrings::readAAMultipleAlignment
  else stop("Unknown sequence type!")
  
  
  Biostrings::writeXStringSet(x, infile, append=FALSE, format="fasta")
  
  system(paste("mafft", param, "--clustalout",
               infile,">", outfile))
  
  reader(outfile, format="clustal")
}

## copied from rMSA, GPL
mgsat.muscle <- function(x, param="") {
  
  ## get temp files and change working directory
  wd <- tempdir()
  dir <- getwd()
  temp_file <- basename(tempfile(tmpdir = wd))
  on.exit({
    file.remove(Sys.glob(paste(temp_file, ".*", sep=""))) 
    setwd(dir)
  })
  setwd(wd)
  
  infile <- paste(temp_file, ".in", sep="")
  outfile <- paste(temp_file, ".aln", sep="")
  reader <- if(is(x, "RNAStringSet")) Biostrings::readRNAMultipleAlignment
  else if(is(x, "DNAStringSet")) Biostrings::readDNAMultipleAlignment
  else if(is(x, "AAStringSet")) Biostrings::readAAMultipleAlignment
  else stop("Unknown sequence type!")
  
  
  Biostrings::writeXStringSet(x, infile, append=FALSE, format="fasta")
  
  system(paste("muscle", param, "-clwstrict",
               "-in", infile, "-out", outfile))
  
  ### muscle is missing a blank line in the output!
  rows <- scan(outfile, what = "", sep = "\n", strip.white = FALSE, 
               quiet = TRUE, blank.lines.skip = FALSE)
  if(rows[[3L]] != "") {
    rows <- c(rows[1:2], "", rows[3:length(rows)])
  }
  cat(rows, file=outfile, sep="\n")
  
  
  ### FIXME: Sequences need to be reordered!
  reader(outfile, format="clustal")
}


mgsat.msa <- function(x,msa.method="Muscle",...) {
  library(Biostrings) #otherwise cannot use rownames(y)
  if(msa.method=="Muscle") {
    y = mgsat.muscle(x,...)
  }
  else if(msa.method=="MAFFT") {
    library(data.table)
    nmap = data.table(orig=names(x))[,tmp:=as.character(.I)]
    names(x) = nmap$tmp
    y = mgsat.mafft(x,...)
    stopifnot(names(y)==nmap$tmp)
    rownames(y) = nmap$orig
  }
  as(y,"MultipleAlignment") #drop derived class in case it is there
}

## Return either shallow or deep copy of MultipleAlignment object with masked rows and/or
## columns dropped as requested. Shallow means that a new MSA object us returned that was
## created by assignment operator and its masked modified according to drop.rows and drop.columns
## arguments. Deep means that new object has underlying sequences modified accordingly and mask
## is empty. type.out ="XStringSet" always implies `deep`.
masked.copy.ali <- function(ali, drop.rows = T, drop.columns = T, deep = F, type.out = c("MultipleAlignment","XStringSet")) {
  ##this to make sure that we can call the AAMultipleAlignment constructor as class(ali) if
  ##the ali object is msaAAMultipleAlignment as returned by msa::msa()
  ali = as(ali,"MultipleAlignment")
  library(Biostrings)
  type.out = type.out[[1]]
  ##What conversions below do:
  ##assignment creates a shallow copy that has independent masking attribute
  ##unmasked() would return AAStringSet with masks stripped
  ##that one could be converted back to AAMultipleAlignment to be able
  ##to set only needed masks on a new object.
  ##
  ##as(ali,"XStringSeq") returns a proper type stringset with masked
  ##regions cut out, which is then converted back to AAMultipleAlignment.
  ##Not clear how many deep copies are actually made in the process.
  if(!drop.rows && !drop.columns && type.out=="XStringSet") {
    ret = unmasked(ali)
  }
  else {
    ret = ali
    if(!drop.rows) {
      rowmask(ret) = NULL
    }
    if(!drop.columns) {
      colmask(ret) = NULL
    }
    if(deep || type.out == "XStringSet") {
      #as(as.character(ret),class(ret)[[1]])
      ##This is a generic way to reconstruct object with the same type as ali.
      ret = as(ret,"XStringSet")
      if(type.out == "MultipleAlignment"){
        ret = call.ctor(ali,ret)
      }
    }
  }
  return (ret)
}

get.mask.of.unmasked.rows.ali <- function(ali) {
  bm = rep(T,nrow(ali))
  bm[as.integer(rowmask(ali))] = F
  bm
}

filter.ali <- function(ali,rn,drop.rows = T, drop.columns = T, deep = F, type.out = c("MultipleAlignment","XStringSet")) {
  row.sel = rownames(ali) %in% rn
  rowmask(ali,append="union") = IRanges(start=!row.sel)
  masked.copy.ali(ali,drop.rows = drop.rows, drop.columns = drop.columns, deep = deep, type.out = type.out)
}


## maskGaps can give an error on alignments w/o any gaps -
## this function supresses such errors and return the
## original alignment unaltered
mask.gaps.ali <- function(ali,...) {
  ok = T
  msg = NULL
  x = tryCatch(maskGaps(ali,...),
               error = function(e) { 
                 ok <<- F
                 msg <<- e$message
               }
  )
  if(!ok) {
    stopifnot(grepl(".*subscript out of bounds.*",msg))
    x = ali
  }
  x
}

ungapped.seq.ali <- function(ali,drop.rows = T,drop.columns = T) {
  ss = masked.copy.ali(ali,
                       drop.rows = drop.rows,
                       drop.columns = drop.columns,
                       type.out = "XStringSet")
  call.ctor(ss,gsub("-","",ss,fixed = T),use.names=T)
}

ungapped.export.ali <- function(ali,file.name,drop.rows = T,drop.columns = T,name = NULL) {
  if(!is.null(name)) {
    rownames(ali) = paste(rownames(ali),name)
  }
  seq = ungapped.seq.ali(ali,drop.rows = drop.rows,drop.columns = drop.columns)
  writeXStringSet(seq,file.name,format = "fasta")
}

remove.end.stop.seq <- function(x,stop.seq=c("*","X"),chop.all=F) {
  ind = (as.character(subseq(x,-1)) %in% stop.seq)
  if(!chop.all) {
    subseq(x[ind],-1) = NULL
  }
  else {
    if(any(ind)) {
      subseq(x,-1) = NULL
    }
  }
  x
}

##If you do not want to count here stops (`*`) at the end of AA sequences,
##remove them prior to this call with remove.end.stops.seq()
nongenic.frequency.seq <- function(x,degen.symb=NULL,out.format=c("rowSums","matrix")) {
  out.format = out.format[[1]]
  if(is.null(degen.symb)) {
    if(inherits(x,"AAStringSet")) {
      degen.symb = "BJZX*+." #see ?AAString for meaning
    }
    else stop("Not implemented yet")
  }
  degen.symb = unlist(strsplit(degen.symb,""))
  res = letterFrequency(x,degen.symb)
  if(out.format=="rowSums") {
    res = rowSums(res)
  }
  res
}

##If threshold==0, take a simple majority vote at each position
get.consensus.from.ali <- function(ali,drop.rows=T,drop.columns=F,threshold=0.5,...) {
  library(Biostrings)
  seq = masked.copy.ali(ali,drop.rows = drop.rows,drop.columns = drop.columns, type.out = "XStringSet")
  if(threshold==0) {
    library(stringr)
    conmat = consensusMatrix(seq,...)
    str_c(rownames(conmat)[max.col(t(conmat))],collapse = "")
  }
  else {
    ## consensusString gives an error on MultipleAlignment object even if it has no mask set
    ## need to give it XStringSet
    consensusString(seq,threshold=threshold,...)
  }
}

## Reorder rows of MultipleSequenceAlignment.
## order.names is a sequence of row names to put first.
## The order of the remaining rows can be optionally randomized.
## Existings masks will be retained.
reorder.rows.ali <- function(ali,order.names,random.rest=T) {
  rmask = rowmask(ali)
  cmask = colmask(ali)
  seq = unmasked(ali)
  rmask.bool = rep(F,length(seq))
  names(rmask.bool) = names(seq)
  rmask.bool[as.integer(rmask)] = T
  headseq = seq[order.names]
  headrmask.bool = rmask.bool[order.names]
  rest.bool = !(names(seq) %in% order.names)
  tailseq = seq[rest.bool]
  tailrmask.bool = rmask.bool[rest.bool]
  if(random.rest) {
    ind.rnd = sample.int(length(tailseq))
    tailseq = tailseq[ind.rnd]
    tailrmask.bool = tailrmask.bool[ind.rnd]
  }
  seq = c(
    headseq,
    tailseq
  )
  rmask.bool = c(
    headrmask.bool,
    tailrmask.bool
  )
  call.ctor(ali,seq,
            rowmask=as(rmask.bool,"NormalIRanges"),
            colmask=cmask)
}

## see vignette phangorn-specials
add.gap.level.phyDat <- function(phdat) {
  contr = attr(phdat,"contrast")
  lev = attr(phdat,"levels")
  lev = c(lev,"-")
  contr = cbind(contr,`-`=0)
  contr["-",] = 0
  contr["-","-"] = 1
  contr["?","-"] = 1
  phyDat(phdat,type="USER",contrast=contr)
}

alignment.as.phyDat <- function(ali,method=c("file","memory"),with.gaps=F) {
  library(phangorn)
  library(Biostrings)
  method = method[[1]]
  if(method == "memory") {
    ##directly; as.matrix() removes masked positions
    mat = Biostrings::as.matrix(ali)
    ##phangorn::as.phyDat is broken for AA unless type argument is
    ##provided, we just use this internal function
    if(inherits(ali,"AAMultipleAlignment")) {
      phdat = phangorn:::phyDat.AA(mat)
    }
    else {
      phdat = phangorn::as.phyDat(mat)
    }
  }
  else if(method == "file") {
    ali.file = tempfile("tmp.phydat_convert.",tmpdir=getwd(),fileext=".fasta")
    writeXStringSet(as(ali,"XStringSet"),ali.file,format = "fasta")
    type = "DNA"
    if(inherits(ali,"AAMultipleAlignment")) {
      type = "AA"
    }
    phdat = phangorn::read.phyDat(ali.file, format="fasta", type=type)
    unlink(ali.file)
  }
  else stop(sprintf("Unknown method: %s",method))
  if(with.gaps) {
    phdat = add.gap.level.phyDat(phdat)
  }
  phdat
}

get.seq.from.phyDat <- function(x,to_upper=T) {
  #library(matrixStats)
  library(Biobase)
  library(stringr)
  cn = attr(x,"levels")
  index = attr(x,"index")
  seq = list()
  ## looping like this we maintain all attributes of x,
  ## as opposed to using lapply and return new list
  for(i in seq_along(x)) {
    y = x[[i]]
    colnames(y) = cn
    rn = cn[max.col(y)]
    ##gap is a row of equal values if gap was not a level in phyDat
    #rn[rowMin(y) == rowMax(y)] = "-"
    rownames(y) = rn
    x[[i]] = y
    seq = c(seq,str_c(rn[index],collapse = ""))
  }
  seq=unlist(seq)
  if(to_upper) seq = str_to_upper(seq)
  names(seq) = names(x)
  list(phdat=x,seq=seq)
}

reconstruct.ancestors <- function(ali,tree) {
  phdat = alignment.as.phyDat(ali,method = "file",with.gaps = T)
  library(phangorn)
  library(ape)
  tree.file = "pratchet.tree"
  if(T || !file.exists(tree.file)) {
    tree.opt = pratchet(phdat,start=tree,maxit = 100,k=5)
    write.tree(tree,tree.file,tree.names = T)
  }
  else {
    tree.opt = read.tree(tree.file,tree.names = T)
  }
  
  tree.opt = acctran(tree.opt, phdat)
  anc.mpr = ancestral.pars(tree.opt, phdat, "ACCTRAN") #"MPR"
  res = get.seq.from.phyDat(anc.mpr)
  res$tree = tree.opt
  res
}

##TODO: use optional tree$node.label (see ?phylo)
get.tree.node.label <- function(tree,node.ind) {
  if(node.ind<=length(tree$tip.label)) {
    node.lab = tree$tip.label[node.ind]
  }
  else {
    node.lab = as.character(node.ind)
  }
  node.lab
}

##TODO: use optional tree$node.label (see ?phylo)
get.tree.node.ind <- function(tree,node.label) {
  ind = seq(length(tree$tip.label)+tree$Nnode)
  names(ind) = c(tree$tip.label,as.character(ind[(length(tree$tip.label)+1):length(ind)]))
  ind[node.label]
}

get.tree.node.medoid <- function(tree,format=c("index","label")) {
  library(ape)
  library(Biobase)
  format = format[[1]]
  d.nd = as.matrix(ape::dist.nodes(tree))
  ## note that passing node number as string to ape::root segfaults R
  node.med = as.numeric(rownames(d.nd)[which.min(rowSums(d.nd))][1])
  if(format=="label") {
    node.med = get.tree.node.label(tree,node.med)
  }
  node.med
}

get.tree.mrca <- function(tree,tips,format=c("index","label")) {
  library(ape)
  node.med = ape::getMRCA(tree,tips)
  if(format=="label") {
    node.med = get.tree.node.label(tree,node.med)
  }
  node.med
}

get.center.of.tree <- function(tree,nodes=NULL,format=c("index","label")) {
  library(ape)
  library(Biobase)
  format = format[[1]]
  if(is.null(nodes)) nodes = seq(length(tree$tip.label))
  else nodes = get.tree.node.ind(tree,nodes)
  d.nd = as.matrix(ape::dist.nodes(tree))
  ## only consider distances to selected nodes
  d.nd = d.nd[,nodes]  
  node.med = as.numeric(rownames(d.nd)[which.min(rowMedians(d.nd))][1])
  if(format=="label") {
    node.med = get.tree.node.label(tree,node.med)
  }
  node.med
}


balance.alignment.rows.by.group <- function(ali,ali.attr,group.attr.name,group.target) {
  ali = masked.copy.ali(ali,drop.rows = T,drop.columns = T,deep = T)
  ali.attr = filter.seq.attr.by.ali(ali,ali.attr,drop.rows=T)
  ali.attr = ali.attr[balanced.sample(ali.attr[[group.attr.name]],
                                      target=group.target,
                                      drop.smaller = F)]
  ali = filter.ali(ali,ali.attr$rn)
  list(ali=ali,ali.attr=ali.attr)
}

filter.seq.attr.by.ali <- function(ali,seq.attr,drop.rows=T) {
  setkey(seq.attr,rn)
  rn = rownames(ali)
  if(drop.rows) rn = rn[get.mask.of.unmasked.rows.ali(ali)]
  seq.attr = seq.attr[rn,nomatch=0]
  stopifnot(nrow(seq.attr)==length(rn))
  seq.attr
}

cluster.with.dist.cutoff <- function(d,h) {
  cl.res = hclust(d)
  gr = cutree(cl.res,h = h)
  stopifnot(all(labels(d)==names(gr)))
  list(group=gr,clust=cl.res)
}

alignment.report <- function (x,
                              caption=NULL,
                              export.to.file=T,
                              ungapped.export.to.file=F,
                              ungapped.drop.rows=T,
                              ungapped.drop.columns=T,
                              show.inline=F,
                              show.first.rows=200,
                              show.first.cols=200,
                              skip.if.empty=F,
                              add.widget.args=list(),
                              name=NULL,
                              ...) {
  if(export.to.file) {
    library(Biostrings)
    name.base=paste(str.to.file.name(caption,20),".fasta",sep="")
    ali.file = report$make.file.name(name.base,make.unique=T)
    s = as(x,"XStringSet")
    if(!is.null(name)) {
      names(s) = paste(names(s),name)
    }
    writeXStringSet(s,ali.file,format = "fasta")
    ali.file.descr = "In FASTA format"
  }
  else {
    ali.file=NULL
  }
  if(ungapped.export.to.file) {
    name.base=paste(str.to.file.name(caption,20),".ungapped.fasta",sep="")
    ungapped.file = report$make.file.name(name.base,make.unique=T)
    ungapped.export.ali(x,ungapped.file,drop.rows=ungapped.drop.rows,drop.columns=ungapped.drop.columns,name=name)
    ungapped.file.descr = "Ungapped sequences in FASTA format."
  }
  else {
    ungapped.file=NULL
  }
  
  library(rbiojsmsa)
  do.call(report$add.widget,
          c(
            list(rbiojsmsa(x,name=name,...),
                 caption = caption,
                 data.export=ali.file,
                 data.export.descr=ali.file.descr,
                 show.inline=show.inline
            ),
            add.widget.args
          )
  )
  if(!is.null(ungapped.file)) {
    report$add.file(ungapped.file,caption=paste(caption,ungapped.file.descr,sep=". "))
  }
}




count_matr_from_df<-function(dat,col_ignore=c()) {
  mask_col_ignore = names(dat) %in% col_ignore
  as.matrix(dat[!mask_col_ignore])
}

split_count_df<-function(dat,col_ignore=c(),rownames_col=NULL) {
  dat = as.data.frame(dat)
  if(is.null(rownames_col)) {
    rn = rownames(dat)
  }
  else {
    rn = dat[[rownames_col]]
    if(!rownames_col %in% col_ignore) {
      col_ignore = c(col_ignore,rownames_col)
    }
  }
  if(is.null(rn)) {
    stop("Row IDs must be defined")
  }
  mask_col_ignore = names(dat) %in% col_ignore
  if(!all(mask_col_ignore)) {
    m = as.matrix(dat[,!mask_col_ignore])
  }
  else {
    m = NULL
  }
  attr = dat[mask_col_ignore]
  rownames(attr) = rn
  rownames(m) = rn
  list(count=m,attr=attr)
}

## operation opposite to split_count_df()
join_count_df<-function(m_a) {
  
  if(!is.null(m_a$count)) {
    count = as.data.frame(m_a$count)
    if(!is.null(m_a$attr)) {
      cbind(count,m_a$attr)
    }
    else {
      count
    }
  }
  else {
    m_a$attr
  }
  
}

new.m_a <- function(count=NULL,attr=NULL,attr_feat=NULL,validate=F) {
  if(is.null(attr)) {
    if(!is.null(count)) {
      attr = data.frame(SampleID=rownames(count))
    }
  }
  if(validate) {
    if(!(is.null(count) || is.null(attr))) {
      stopifnot(all(rownames(count)==rownames(attr)))
    }
    if(!(is.null(count) || is.null(attr_feat))) {
      stopifnot(all(colnames(count)==rownames(attr_feat)))
    }
  }
  list(count=count,attr=attr,attr_feat=attr_feat)
}

as.m_a.dds <- function(dds,do.rlog.norm,...) {
  library(DESeq2)
  if(!do.rlog.norm) {
    count = t(as.matrix(assay(dds)))
  }
  else {
    count = rlog.dds(dds,...)
  }
  attr = as.data.frame(colData(dds))
  attr_feat = as.data.frame(rowData(dds))
  m_a = new.m_a(count=count,attr=attr,attr_feat = attr_feat)
}

as.dds.m_a <- function(m_a,formula.rhs,force.lib.size=T,round.to.int=T) {
  
  library(DESeq2)
  if(round.to.int) {
    m_a$count = round(m_a$count)
  }
  ## DESeq2 does not work on ordered factors; here we convert all orderd
  ## factors used in the formula into unordered but with the same order of levels, 
  ## in a copy of the colData
  design = as.formula(paste("~",formula.rhs))
  colData = data.table::copy(m_a$attr) #in case this is not a DF but DT, works for DF too
  for(v in all.vars(design)) {
    if(is.ordered(colData[[v]])) {
      colData[[v]] = factor(colData[[v]],ordered = F,levels=levels(colData[[v]]))
    }
  }
  dds <- DESeqDataSetFromMatrix(countData = t(m_a$count),
                                colData = colData,
                                design = design)  
  if(force.lib.size) {
    ## from phyloseq vignette at 
    ## http://www.bioconductor.org/packages/release/bioc/vignettes/phyloseq/inst/doc/phyloseq-mixture-models.html
    gm_mean = function(x, na.rm=TRUE){
      exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
    }
    geoMeans = apply(counts(dds), 1, gm_mean)
    dds = estimateSizeFactors(dds, geoMeans = geoMeans)
  }
  
  return(dds)
}

norm.prop <- function(x, ...) {
  UseMethod('norm.prop', x)
}

norm.prop.default <- function(x.f) {
  x.f/sum(x.f)
}

norm.prop.matrix<-function(x.f,mar=1) {
  if(mar==1)
    x.f/rowSums(x.f)
  else if(mar==2)
    x.f/colSums(x.f)
  else
    stop("mar can be only 1 or 2")
}

norm.prop.m_a<-function(m_a,mar=1) {
  m_a$count = norm.prop(m_a$count,mar=mar)
  return (m_a)
}

norm.prop.data.frame<-function(dat,col_ignore=c(),mar=1) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  matnorm<-norm.prop(x$count,mar=mar)
  datnorm<-as.data.frame(matnorm) 
  res = cbind(x$attr,datnorm)
  #print(apply(res[!(names(res) %in% col_ignore)],1,sum) == 1)
  stopifnot(all(abs(apply(res[!(names(res) %in% col_ignore)],1,sum) - 1)<1e-6))
  res
}

norm.meta.data<-function(dat,col_ignore=c(),norm.func=NULL,...) {
  mask_col_ignore = names(dat) %in% col_ignore
  x<-split_count_df(dat,col_ignore)
  if( is.null(norm.func) ) {
    norm.func = ihs
  }
  matnorm<-norm.func(x$count,...)
  datnorm<-as.data.frame(matnorm) 
  cbind(x$attr,datnorm)
}

# CLR and ALR transforms are copied from SpiecEasi package
# https://github.com/zdk123/SpiecEasi
# If data is non-normalized count OTU/data table with samples on rows
# and features/OTUs in columns, then the transform is applied as
# clr(data)
# By default, it does not add any offset of before applying the transform, and acts
# on rows, returning matrix in the same order as input. Adding offset=1 might
# be a good idea when you are dealing with integer count data where counts are 
# generally much higher than one, and a count of one is close to noise level.
# Here is why: the implementation will just ignore all features that are
# close to zero within the tolerance. Features with count one, on the other hand,
# will contribute zero to the nominator of the mean, but zero to the divider when the mean is taken.
# Thus, at lower depth of sequencing, you will get more features which are exactly zero and do
# not increase the divider, resulting in overal higher level of all non-zero features after
# transformation. Adding offset one removes this effect, and should be OK under an 
# assumption that one count should be considered random noise that should never affect
# the analysis results. Another option is to use offset=0.5 (uninformative prior)
# The defaul base=2 in
# order to produce fold change between columns:
# (m[,k]/m[,l] == 2**(clr.mgsat(m)[,k]-clr.mgsat(m)[,l]))
# Note that in the original SpiecEasi implementation the transform should
# be applied as
# t(clr(data+1, 1, base=2))
# because it uses apply(mar=1) which transposes the result, and default base=exp(1)
# See `compositions` package for other Aitchison transforms
# Note: CLR as implemented here might be quite sensitive to the "depth of sequencing" 
# ("collection effort"). At larger depth and uneven abundance distribution, you get 
# much longer tail with small but non-zero values. Depending on the domain of you data
# (e.g. integer counts or already proportions) and your offset and tol values, samples
# with higher depth of sequencing may get very negative mean(log), dominated by the long tail with
# low values

#' Centered log-ratio functions
#' @export
norm.clr <- function(x, ...) {
  UseMethod('norm.clr', x)
}

#' @method clr default
#' @export
norm.clr.default <- function(x.f, offset=0, base=2, tol=.Machine$double.eps) {
  ## this is invariant to a constant multiplier (well, not quite becase of the
  ## offset), so there is no need to combine it with normalization to simple
  ## proportions
  #x.f = ifelse(x.f>0,x.f,offset)
  x.f = x.f + offset
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  ifelse(nzero, LOG - mean(LOG)/mean(nzero), 0.0)
}

#' @method clr matrix
#' @export
norm.clr.matrix <- function(x.f, mar=1, ...) {
  y = plyr::aaply(x.f, mar, norm.clr, ...)
  if(mar==2) {
    y = t(y)
  }
  return(y)
}

#' @method clr data.frame
#' @export
norm.clr.data.frame <- function(x.f, ...) {
  as.data.frame(norm.clr(as.matrix(x.f), ...))
}

#' Additive log-ratio functions
#' @export
norm.alr <- function(x, ...) {
  UseMethod('norm.alr', x)
}


norm.alr.default <- function(x.f, divcomp=1, offset=1, base=2, remove.divcomp=TRUE,
                             tol=.Machine$double.eps) {
  if(is.character(divcomp)) {
    divcomp = match(divcomp,colnames(x.f))
  }
  #x.f = ifelse(x.f>0,x.f,offset)
  x.f = x.f + offset
  nzero <- (x.f >= tol)
  LOG <- log(ifelse(nzero, x.f, 1), base)
  x.alr <- ifelse(nzero, LOG - LOG[divcomp], 0.0)
  if (remove.divcomp) x.alr[-divcomp]
  else x.alr
}

norm.alr.matrix <- function(x.f, mar=1, divcomp=1, remove.divcomp=T, ...) {
  ##TODO: optimize by doing matrix operations here instead of plyr::aaply
  if(is.character(divcomp)) {
    name.divcomp = divcomp
    divcomp = match(divcomp,colnames(x.f))
    if(length(divcomp)!=1) {
      stop(sprintf("One and one only feature name should match the divider component name %s",name.divcomp))
    }
  }
  y = plyr::aaply(x.f, mar, norm.alr, divcomp=divcomp, remove.divcomp=F, ...)
  if(mar==2) {
    y = t(y)
  }
  if(remove.divcomp) {
    y = y[,-divcomp]
  }
  return(y)
}

norm.alr.data.frame <- function(x.f, ...) {
  as.data.frame(norm.alr(as.matrix(x.f), ...))
}

#' Fold-ratio functions
#' @export
norm.fr <- function(x, ...) {
  UseMethod('norm.fr', x)
}


norm.fr.default <- function(x.f, divcomp=1, offset=1, remove.divcomp=TRUE,
                            tol=.Machine$double.eps,offset.is.prop=F) {
  if(is.character(divcomp)) {
    divcomp = match(divcomp,colnames(x.f))
  }
  if(offset.is.prop) {
    s = sum(x.f)
    if(abs(s) < tol) {
      if (remove.divcomp) x.f = x.f[-divcomp]
      return (x.f)
    }
    x.f = x.f / s
  }
  x.f = x.f + offset
  x.fr <- x.f/x.f[divcomp]
  if (remove.divcomp) x.fr[-divcomp]
  else x.fr
}

norm.fr.matrix <- function(x.f, mar=1, divcomp=1, remove.divcomp=T, ...) {
  ##TODO: optimize by doing matrix operations here instead of plyr::aaply
  if(is.character(divcomp)) {
    name.divcomp = divcomp
    divcomp = match(divcomp,colnames(x.f))
    if(length(divcomp)!=1) {
      stop(sprintf("One and one only feature name should match the divider component name %s",name.divcomp))
    }
  }
  y = plyr::aaply(x.f, mar, norm.fr, divcomp=divcomp, remove.divcomp=F, ...)
  if(mar==2) {
    y = t(y)
  }
  if(remove.divcomp) {
    y = y[,-divcomp]
  }
  return(y)
}

norm.fr.data.frame <- function(x.f, ...) {
  as.data.frame(norm.fr(as.matrix(x.f), ...))
}


#' @method clr on results of DESeq2 Rlog transform. x.f is not used
norm.clr.rlog.dds <- function(x.f,dds,...) {
  plyr::aaply(rlog.dds(dds,...),1,function(x) (x - mean(x)))
}

#' @method results of DESeq2 Rlog transform. x.f is not used
norm.rlog.dds <- function(x.f,dds,...) {
  rlog.dds(dds,...)
}

#' @method extract Rlog transformed values from DESeq2 object as sample-row matrix
#' Note the discussion of the proper value for blind in the DESeq2::rlog help page:
#' "blind=FALSE should be used for transforming data for downstream analysis, 
#' where the full use of the design information should be made".
rlog.dds <- function(dds,blind=T,fitType="local",...) {
  t(assay(rlog(dds,blind=blind,fitType=fitType,...)))
}


## IHS (inverse hyperbolic sign) transform
## asinh(x) == log(x+(x**2+1)**0.5)
## Increasing theta reduces the influence of unit addition
## under the square root and make the function behave more
## like (A + log(x)/B). Generally, for variance stabilization,
## it is assumed that x >> 1. Thus, if you are applying it to
## simple proportions, set theta to some large number (~100).
ihs <- function(x,theta=1) {
  asinh(theta*x)/theta
}

#Function for symbolic derivative of arbitrary order
#Copied from help page for 'deriv()'
#Use as:
#dd.expr = DD(expression(log(x+(x**2+1)**0.5)),"x",2)
#To create a function from the 'expression' output:
#f = function(x) eval(dd.expr)
DD <- function(expr, name, order = 1) {
  if(order < 1) stop("'order' must be >= 1")
  if(order == 1) D(expr, name)
  else DD(D(expr, name), name, order - 1)
}

## Return a function that computes the derivative of expr
make.DD <- function(expr,name,order = 1) {
  dn.expr = DD(expr,name,order=order)
  return (function(x) {eval(dn.expr)})
}

ihs.d1 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",1)
ihs.d2 = make.DD(expression(log(x+(x**2+1)**0.5)),"x",2)

norm.ihs = ihs

## Boxcox transform
boxcox <- function(x,lambda1,lambda2=0) {
  #print(paste("l1=",lambda1,"l2=",lambda2))
  if (lambda1 != 0) {
    ((x+lambda2)**lambda1-1)/lambda1
  }
  else {
    log(x+lambda2)
  }
}

## Copied from geoR:: to avoid having to install tcltk
"boxcoxfit" <-
  function(object, xmat, lambda, lambda2 = NULL, add.to.data = 0,...)
  {
    call.fc <- match.call()
    data <- object + add.to.data
    if(is.null(lambda2) && any(data <= 0))
      stop("Transformation requires positive data")
    ##
    data <- as.vector(data)
    n <- length(data)
    if(missing(xmat)) xmat <- rep(1, n)
    xmat <- as.matrix(xmat)
    if(any(xmat[,1] != 1)) xmat <- cbind(1, xmat)
    ## do not reverse order of the next two lines:
    xmat <- xmat[!is.na(data),,drop=FALSE]
    data <- data[!is.na(data)]
    n <- length(data)
    ##
    beta.size <- ncol(xmat)
    if(nrow(xmat) != length(data))
      stop("xmat and data have incompatible lengths")
    ##  lik.method <- match.arg(lik.method, choices = c("ML", "RML"))
    lik.method <- "ML"
    ##
    if(all(data > 0)) absmin <- 0
    else absmin <- abs(min(data)) + 0.00001 * diff(range(data))
    if(!is.null(lambda2)){
      if(missing(lambda)) lambda.ini <- seq(-2, 2, by=0.2)
      else lambda.ini <- lambda
      lambda2.ini <- 0
      if(isTRUE(lambda2)) lambda2.ini <- absmin
      else if(mode(lambda2) == "numeric") lambda2.ini <- lambda2
      lambdas.ini <- as.matrix(expand.grid(lambda.ini, lambda2.ini))
      ##
      if(length(as.matrix(lambdas.ini)) > 2){
        lamlik <- apply(lambdas.ini, 1, .negloglik.boxcox, data=data + absmin,
                        xmat=xmat, lik.method=lik.method)
        lambdas.ini <- lambdas.ini[which(lamlik == min(lamlik)),]
      }
      lambdas.ini <- unname(drop(lambdas.ini))
      lik.lambda <- optim(par=lambdas.ini, fn = .negloglik.boxcox,
                          method="L-BFGS-B",
                          #hessian = TRUE, 
                          lower = c(-Inf, absmin), 
                          data = data, xmat = xmat, lik.method = lik.method)
    }
    else{
      lik.lambda <- optimize(.negloglik.boxcox, interval = c(-5, 5), data = data,
                             xmat = xmat, lik.method = lik.method)
      lik.lambda <- list(par = lik.lambda$minimum, value = lik.lambda$objective,
                         convergence = 0, message = "function optimize used")
    }
    ##
    ##  hess <- sqrt(diag(solve(as.matrix(lik.lambda$hessian))))
    lambda.fit <- lik.lambda$par
    if(length(lambda.fit) == 1) lambda.fit <- c(lambda.fit, 0)
    data <- data + lambda.fit[2]
    ##
    #  if(abs(lambda.fit[1]) < 0.0001) yt <- log(data)
    if(isTRUE(all.equal(unname(lambda.fit[1]),0))) yt <- log(data)
    else yt <- ((data^lambda.fit[1]) - 1)/lambda.fit[1]
    beta <- solve(crossprod(xmat), crossprod(xmat, yt))
    mu <- drop(xmat %*% beta)
    sigmasq <- sum((yt - mu)^2)/n
    if(lik.method == "ML")
      loglik <- drop((-(n/2) * (log(2*pi) + log(sigmasq) + 1)) + (lambda.fit[1]-1) * sum(log(data)))
    ##  if(lik.method == "RML")
    ##    loglik <- drop(-lik.lambda$value - (n/2)*log(2*pi) - (n-beta.size)*(log(n) - 1))
    ##
    temp <- 1 + lambda.fit[1] * mu
    fitted.y <- ((temp^((1/lambda.fit[1]) - 2)) *
                   (temp^2 + ((1-lambda.fit[1])/2) * sigmasq))
    variance.y <-  (temp^((2/lambda.fit[1]) - 2)) * sigmasq
    if(beta.size == 1){
      fitted.y <- unique(fitted.y)
      variance.y <- unique(fitted.y)
    }
    ##
    beta <- drop(beta)
    if(length(beta) > 1)
      names(beta) <- paste("beta", 0:(beta.size-1), sep="")
    if(length(lik.lambda$par) == 1) lambda.fit <- lambda.fit[1]
    if(length(lik.lambda$par) == 2) names(lambda.fit) <- c("lambda", "lambda2")
    res <- list(lambda = lambda.fit, beta.normal = drop(beta),
                sigmasq.normal = sigmasq, 
                loglik = loglik, optim.results = lik.lambda)
    ## res$hessian <- c(lambda = hess) 
    res$call <- call.fc
    oldClass(res) <- "boxcoxfit"
    return(res)
  }

## Copied from geoR:: to avoid having to install tcltk
".negloglik.boxcox" <-
  function(lambda.val, data, xmat, lik.method = "ML")
  {
    if(length(lambda.val) == 2){
      data <- data + lambda.val[2]
      lambda <- lambda.val[1]
    }
    else lambda <- lambda.val
    lambda <- unname(lambda)
    n <- length(data)
    beta.size <- ncol(xmat)
    if(isTRUE(all.equal(unname(lambda), 0))) yt <- log(data)
    else yt <- ((data^lambda) - 1)/lambda
    beta <- solve(crossprod(xmat), crossprod(xmat, yt))
    ss <- sum((drop(yt) - drop(xmat %*% beta))^2)
    if(lik.method == "ML")
      neglik <- (n/2) * log(ss) - ((lambda - 1) * sum(log(data)))
    if(lik.method == "RML"){
      xx <- crossprod(xmat)
      if(length(as.vector(xx)) == 1)
        choldet <- 0.5 * log(xx)
      else
        choldet <- sum(log(diag(chol(xx))))
      neglik <- ((n-beta.size)/2) * log(ss) + choldet -
        ((lambda - 1) * sum(log(data)))
    }
    if(mode(neglik) != "numeric") neglik <- Inf
    return(drop(neglik))
  }


## Fit boxcox transform to the data and transform the data
boxcox.transform.vec <- function(x) {
  l.2 = T
  if(!all(x>0)) {
    b = boxcoxfit(x,lambda2=T)
  }
  else {
    b = boxcoxfit(x)
    l.2 = F
  }
  #b = try(boxcoxfit(x,lambda2=T),TRUE)
  #for some reason the above fails if all(b>0)
  #we just repeat it with lambda2 NULL
  #if (inherits(b, "try-error")) {
  #  b = boxcoxfit(x)
  #  l.2 = T
  #  print("err")
  #  print(b)
  #}
  lambda1 = NA
  lambda2 = 0
  if (l.2) { 
    lambda1 = b$lambda["lambda"] 
    lambda2 = b$lambda["lambda2"] 
  }
  else {
    lambda1 = b$lambda
  }
  return (list(x=boxcox(x,lambda1=lambda1,lambda2=lambda2),
               boxcox=b))
}

norm.boxcox <- function(x, ...) {
  UseMethod('norm.boxcox', x)
}

norm.boxcox.default <- function(x.f) {
  boxcox.transform.vec(x.f)$x
}

norm.boxcox.matrix <- function(x.f, mar=2, ...) {
  y = plyr::aaply(x.f, mar, norm.boxcox, ...)
  if(mar==2) {
    y = t(y)
  }
  return(y)
}

norm.ihs.prop <- function(x,theta=1,mar=1) {
  norm.ihs(norm.prop(x,mar=mar),theta=theta)
}

norm.boxcox.prop <- function(x,mar.prop=1,mar.boxcox=2) {
  norm.boxcox(norm.prop(x,mar=mar.prop),mar=mar.boxcox)
}


## normalize raw count data according to one of the
## methods defined above.

norm.count <- function(x, ...) {
  UseMethod('norm.count', x)
}

norm.count.matrix <- function(count,method,drop.features=c("other"),method.args=list()) {
  if( ! (method %in% c("norm.ident","ident")) ) {
    count = do.call(method,c(
      list(count),
      method.args
    )
    )
  }
  if(!is.null(drop.features)) {
    count = count[,!colnames(count) %in% drop.features,drop=F]
  }
  return (count)
}

norm.count.m_a <- function(m_a,...) {
  m_a.norm = m_a
  m_a.norm$count = norm.count(m_a.norm$count,...)
  return(m_a.norm)
}

## Jensen-Shannon distance
## Adapted from http://enterotype.embl.de/enterotypes.html
dist.js <- function(x, offset=.Machine$double.eps, ...) {
  x=t(x)
  KLD <- function(x,y) sum(x *log(x/y))
  JSD<- function(x,y) sqrt(0.5 * KLD(x, (x+y)/2) + 0.5 * KLD(y, (x+y)/2))
  matrixColSize <- length(colnames(x))
  matrixRowSize <- length(rownames(x))
  colnames <- colnames(x)
  resultsMatrix <- matrix(0, matrixColSize, matrixColSize)
  
  x = ifelse(x==0,offset,x)
  
  for(i in 1:matrixColSize) {
    for(j in 1:matrixColSize) { 
      resultsMatrix[i,j]=JSD(as.vector(x[,i]),
                             as.vector(x[,j]))
    }
  }
  colnames -> colnames(resultsMatrix) -> rownames(resultsMatrix)
  as.dist(resultsMatrix)->resultsMatrix
  attr(resultsMatrix, "method") <- "dist"
  return(resultsMatrix) 
}

names.to.ind <- function(x,ind) {
  if(!is.numeric(ind)) {
    ind = as.character(ind)
    ind.names = seq(nrow(x))
    names(ind.names) = rownames(x)
    ind = ind.names[ind]
  }  
  ind
}

## Take sample x feature matrix and id.ref that selects 
## "reference" samples in that matrix, and return distance
## to those samples from the remaining samples.
## id.ref can be either a scalar, in which case a distance to
## that single sample from all remaining samples is computed,
## or a vector of length nrow(x), with each element pointing to
## some (possibly repeating elements of x.
## id.ref elements are excluded from the returned vector.
## TODO: currently a square matrix all-against-all is computed and then
## subset is returned. It should be done instead only for requested
## combinations.
dist.to.reference <- function(x,id.ref,dist.method,drop.ref=T,drop.na=T) {
  if(dist.method=="js") {
    ##Jensen-Shannon
    dist.o = dist.js(x)
  }
  else if(dist.method=="hell") {
    ##Hellinger
    dist.o = vegdist(sqrt(x),method = "euclidean")/sqrt(2)
  }
  else {
    dist.o = vegdist(x,method = dist.method)
  }
  dist.o = as.matrix(dist.o)
  id.ref = names.to.ind(x,id.ref) 
  ## index will be a two column matrix
  ind = as.matrix(data.frame(seq(nrow(x)),id.ref))
  res = dist.o[ind]
  names(res) = rownames(x)
  if(drop.ref) {
    res = res[-id.ref[!is.na(id.ref)]]
  }
  if(drop.na) {
    res = res[!is.na(res)]
  }
  res
}

## Similar to dist.to.reference but returns a matrix of profile-profile subtractions (contrasts)
contrast.to.reference <- function(x,id.ref,drop.ref=T) {
  id.ref = names.to.ind(x,id.ref)
  if(length(id.ref)==1 && nrow(x)>1) {
    id.ref = rep(id.ref,nrow(x))
  }
  x.ref = x[id.ref,]
  res = x - x.ref
  rownames(res) = rownames(x)
  if(drop.ref) {
    res = res[-id.ref[!is.na(id.ref)],]
  }
  res
}

rownames.set.m_a <- function(m_a,x,set.attr.field=T) {
  x = as.character(x)
  rownames(m_a$count) = x
  rownames(m_a$attr) = x
  if(set.attr.field) {
    m_a$attr$SampleID = x
  }
  m_a
}

## subset method that will use the same subset argument on all data objects in m_a
subset.m_a <- function(m_a,subset=NULL,select.count=NULL,select.attr=NULL,na.index.is.false=T) {
  
  set.index.na.to.false <- function(ind) {
    if(all(sapply(ind,is.logical))) {
      ind[is.na(ind)] = F
    }
    else {
      ind = ind[!is.na(ind)]
    }
    ind
  }
  
  if(is.null(select.count)) select.count = T
  if(is.null(select.attr)) select.attr = T
  if(is.null(subset)) subset = T
  
  if(na.index.is.false) {
    select.count = set.index.na.to.false(select.count)
    select.attr = set.index.na.to.false(select.attr)
    subset = set.index.na.to.false(subset)
  }
  
  m_a$count = m_a$count[subset,select.count,drop=F]
  m_a$attr = m_a$attr[subset,select.attr,drop=F]
  
  ## remove no longer used levels from all factors
  for(colnam in colnames(m_a$attr)) {
    if(is.factor(m_a$attr[,colnam])) {
      m_a$attr[,colnam] = factor(m_a$attr[,colnam])
    }
  }
  
  if(!is.null(m_a$attr_feat)) {
    m_a$attr_feat = m_a$attr_feat[select.count,,drop=F]
    for(colnam in colnames(m_a$attr_feat)) {
      if(is.factor(m_a$attr_feat[,colnam])) {
        m_a$attr_feat[,colnam] = factor(m_a$attr_feat[,colnam])
      }
    }
  }
  
  return(m_a)
}

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

long_data.to.m_a <- function(dat,
                             col.attr,
                             col.count,
                             fun.aggregate = function(x) mean(x,na.rm = T),
                             fill=0,
                             value.var = "Value",
                             row.names=NULL,
                             ...
) {
  library(data.table)
  
  dat.cast = dcast.data.table(as.data.table(dat),sprintf("%s~%s",paste(col.attr,collapse="+"),paste(col.count,collapse="+")),
                              fun.aggregate = fun.aggregate,fill=fill,value.var = value.var,...)
  
  dat.cast = setDF(dat.cast)
  if(!is.null(row.names)) {
    rownames(dat.cast) = eval(parse(text=row.names),dat.cast)
    dat.cast$SampleID = rownames(dat.cast)
  }
  else {
    rownames(dat.cast) = paste("S",1:nrow(dat.cast),sep="")
    dat.cast$SampleID = rownames(dat.cast)
  }
  split_count_df(dat.cast,col_ignore = unique(c("SampleID",col.attr)))
  
}

quant.mask <- function(x,prob,drop.zero=T) {
  if(drop.zero) {
    x.q = x[x>0]
  }
  else {
    x.q = x
  }
  x.cut = quantile(x.q, probs=prob)
  return (x>=x.cut)
}

## If the other_cnt column is already present, it will be incremented with counts of features
## relegated to the "other" in this call; otherwise, the new column with this name will be
## created.
## Count columns will be sorted in decreasing order of the column mean frequencies.
## Various column exclusion criteria are applied cumulatively. Row-wise fractions used by some of
## the exclusion criteria are computed on the *original* set of columns.
## All columns that are excluded are summed into the "other" category, which is placed last 
## regardless of its counts.
## drop.names can contain the name of the "other" category, in which case it will be dropped as well.
## n.top (if > 0) will return only top (by mean frequency) features *after* other exclusions, and
## not counting "other", which will be also returned.
count.filter.m_a <- function(m_a,
                             min_max_frac=0.0,
                             min_max=0,
                             min_mean=0,
                             min_mean_frac=0.0,
                             min_quant_mean_frac=0.0,
                             min_incidence_frac=0.0,
                             min_quant_incidence_frac=0.0,
                             min_row_sum=0,
                             max_row_sum=.Machine$integer.max,
                             other_cnt="other",
                             keep.names=NULL,
                             drop.except.names=NULL,
                             drop.names=NULL,
                             n.top=0,
                             drop.zero=T,
                             drop.unclassified=F) {
  ##Note that filtering columns by a median value would not be a good idea - if we have a slightly
  ##unbalanced dataset where one group is 60% of samples and has zero presence in some column,
  ##and another group is 40% and has a large presence, then median filter will always throw this
  ##column away.
  #m.call = match.call()
  ##make.global(m.call)
  #stop("DEBUG")
  x = m_a
  
  ## use this to drop rows and to find the count in "other" category at the end
  row_cnt = rowSums(x$count)
  
  ## index of rows to drop
  row_sel = row_cnt >= min_row_sum & row_cnt < max_row_sum
  
  ## drop rows from data and attribute matrices and derived matrices/vectors
  cnt = x$count[row_sel,,drop=F]
  row_cnt = row_cnt[row_sel,drop=F]
  attr = x$attr[row_sel,,drop=F]
  
  ## get proportions
  cnt_norm = norm.prop(cnt)
  
  ## reorder all columns
  ind_col_ord = order(colSums(cnt_norm),decreasing=T)
  cnt = cnt[,ind_col_ord,drop=F]
  cnt_norm = cnt_norm[,ind_col_ord,drop=F]
  
  ## bitmask of columns to keep; start with keep all
  mask_col_sel = rep(T,ncol(cnt))
  
  if(is.function(drop.except.names)) {
    drop.except.names = do.call(drop.except.names,list(cnt,cnt_norm,attr))
  }
  
  if(is.function(drop.names)) {
    drop.names = do.call(drop.names,list(cnt,cnt_norm,attr))
  }
  
  if(is.function(keep.names)) {
    keep.names = do.call(keep.names,list(cnt,cnt_norm,attr))
  }
  
  ## if drop.except.names defined, only keep names provided where
  if(!is.null(drop.except.names)) {
    mask_col_sel = mask_col_sel & (colnames(cnt) %in% drop.except.names)
  }
  
  if(drop.unclassified) {
    drop.names = c(drop.names,colnames(cnt)[grepl("^Unclassified.*",colnames(cnt),ignore.case=T)])
  }
  
  ## if drop.names defined, drop all names from it; this means drop.names
  ## overrides drop.except.names
  if(!is.null(drop.names)) {
    mask_col_sel = mask_col_sel & !(colnames(cnt) %in% drop.names)
  }  
  
  ## dropping columns with various criteria, each criteria is applied column-wise
  ## (does not depend on other columns). Building mask of columns to *keep*.
  mask_col_sel = mask_col_sel & apply(cnt_norm,2,max) >= min_max_frac
  #cnt = cnt[,ind_col_sel,drop=F]
  #cnt_norm = cnt_norm[,ind_col_sel,drop=F]
  mask_col_sel = mask_col_sel & (apply(cnt_norm,2,mean) >= min_mean_frac)
  #cnt = cnt[,ind_col_sel,drop=F]
  #cnt_norm = cnt_norm[,ind_col_sel,drop=F]
  if(drop.zero) {
    mask_col_sel = mask_col_sel & (!apply(cnt==0,2,all))
    #cnt = cnt[,!apply(cnt==0,2,all),drop=F]
  }
  if(min_quant_mean_frac > 0) {
    frac.sums = colSums(cnt_norm)
    mask_col_sel = mask_col_sel & quant.mask(frac.sums,min_quant_mean_frac,drop.zero=drop.zero)
  }
  mask_col_sel = mask_col_sel & (apply(cnt,2,max) >= min_max)
  #cnt = cnt[,apply(cnt,2,max) >= min_max,drop=F]
  mask_col_sel = mask_col_sel & (apply(cnt,2,mean) >= min_mean)
  #cnt = cnt[,apply(cnt,2,mean) >= min_mean,drop=F]
  mask_col_sel = mask_col_sel & (apply(cnt>0,2,mean) >= min_incidence_frac)
  
  if(min_quant_incidence_frac > 0) {
    frac.inc = apply(cnt>0,2,mean)
    mask_col_sel = mask_col_sel & quant.mask(frac.inc,min_quant_incidence_frac,drop.zero=drop.zero)
  }
  
  ## if keep.names defined, keep all names listed in it; this means keep.names
  ## overrides everything else
  if(!is.null(keep.names)) {
    mask_col_sel = mask_col_sel | (colnames(cnt) %in% keep.names)
  } 
  
  ## drop all columns accumulated so far; not updating cnt_norm because not needed anymore
  cnt = cnt[,mask_col_sel,drop=F]
  cnt_norm = NULL
  
  ## cut to top by abundance columns
  if(n.top>0) {
    cnt = cnt[,1:min(n.top,ncol(cnt)),drop=F]
  }
  
  if(is.null(drop.names) || (! (other_cnt %in% drop.names) ) || (other_cnt %in% keep.names)) {
    
    cnt_col_other = as.matrix(row_cnt - rowSums(cnt))
    
    if (!all(abs(cnt_col_other[,1,drop=F])<=sqrt(.Machine$double.eps)) && !is.null(other_cnt)) {
      colnames(cnt_col_other) = c(other_cnt)
      if (other_cnt %in% colnames(cnt)) {
        cnt[,other_cnt] = cnt[,other_cnt] + cnt_col_other[,other_cnt]
      }
      else {
        cnt = cbind(cnt,cnt_col_other)
      }
    }
  }
  
  x$attr = attr
  x$count = cnt
  return(x)
}

count.filter<-function(dat,
                       col_ignore=c(),
                       ...) {
  ##Note that filtering columns by a median value would not be a good idea - if we have a slightly
  ##unbalanced dataset where one group is 60% of samples and has zero presence in some column,
  ##and another group is 40% and has a large presence, then median filter will always throw this
  ##column away.
  #m.call = match.call()
  ##make.global(m.call)
  #stop("DEBUG")
  x<-split_count_df(dat,col_ignore)
  x = count.filter.m_a(x,...)
  return (join_count_df(x))
}

sample.rows<-function(x,size,replace=F,prob=NULL) {
  x[sample(dim(x)[1],size,replace,prob),]
}

list_to_df<-function(x,col_names=NULL,row_names=NULL) {
  y = as.data.frame(do.call(rbind, x),row.names=row_names)
  if (!is.null(col_names)) {
    names(y) = col_names
  }
  y
}

tryCatchAndWarn<-function(expr,catch.warnings=F,catch.errors=T) {
  if(getOption("try.debug",F)) {
    message("In debug mode, not catching anything in tryCatchAndWarn")
    return (eval.parent(expr))
  }
  ok=T
  #somehow if both warning and error arguments are used
  #in a single call to tryCatch, then error handler is not
  #called - probably only the handler for whatever happens first
  #is called
  
  if(catch.warnings) {
    warnHand = function(w) {
      warning(paste("Warning caught in tryCatch; returning NULL: ",w,"\n"))
      ok<<-F
    }
  }
  
  if(catch.errors) {
    errHand = function(e) {
      warning(paste("Error caught in tryCatch; returning NULL: ",e,"\n"))
      ok<<-F
    }
  }
  
  if(catch.warnings && catch.errors) {
    res = tryCatch(tryCatch(expr,warning=warnHand),error=errHand)
  }
  else if(catch.warnings && !catch.errors) {
    res = tryCatch(expr,warning=warnHand)
  }
  else if(!catch.warnings && catch.errors) {
    res = tryCatch(expr,error=errHand)
  }
  else {
    res = tryCatch(expr)
  }
  
  #options(warn=warn_saved)
  if (!ok) {
    #print("Error occured in test, returning NA")
    res = NULL
  }
  return(res)
}

sorted.freq<-function(count.df,col_ignore=c(),col_group="group") {
  freq = norm.prop(count.df,col_ignore=col_ignore)
  freq_m = count_matr_from_df(freq,col_ignore=col_ignore)
  group = count.df[,col_group]
  freq.mean = aggregate(freq_m,list(group),mean,simplify=TRUE)
  row.names(freq.mean) = freq.mean[,"Group.1"]
  freq.mean$Group.1 = NULL
  freq.mean = t(as.matrix(freq.mean))
  return(freq.mean)
}


kelvin.rnames.to.ids<-function(row.names) {
  col_names = c("id_repl","group")
  ids = list_to_df(str_split(row.names,"_",2),row_names=row.names(data),col_names=col_names)
}

read.kelvin.summary.matr<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="Sample")
  stopifnot(all(data$Total==rowSums(subset(data,select=c(-Total)))))
  data$Total = NULL
  ids=kelvin.rnames.to.ids(row.names(data))
  data$id_repl = ids$id_repl
  data$group = ids$group
  data
}

power.kelvin<-function() {
  group_sel = c("Stool","Saliva")
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls")
  
  data = data_all[data_all$group %in% group_sel,]
  
  ##x = Data.filter(data$x, "sample", 1000, 10)
  data = count.filter(data,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)
  cnt_m = count_matr_from_df(data,col_ignore=attr_names)
  
  
  
  x1 = cnt_m[data$group == group_sel[1],]
  x2 = cnt_m[data$group == group_sel[2],]
  
  n_samp = 15
  
  n_iter = 10
  
  alpha = 0.05
  
  p_vals = c()
  for (i_iter in seq(n_iter)) {
    mygroup = list(sample.rows(x1,n_samp),sample.rows(x2,n_samp))
    p_vals = c(p_vals,Xmcupo.sevsample(mygroup,dim(x1)[2])[[2]])
  }
  power = mean(p_vals<=alpha,na.rm = T)
  
  print(power)
  list(x1,x2)
}


read.koren<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="Taxon")
  data = as.data.frame(t(data))
  rnames = row.names(data)
  data$time = as.ordered(as.integer(substr(rnames,3,3)))
  data$id_repl = as.factor(as.integer(substr(rnames,4,300)))
  data
}

wilcox.test.multi <- function(mat,group=NULL,
                              type=c("unpaired","paired","onesample"),
                              impl=c("wilcox.exact","wilcox.coin")) {
  type = type[1]
  impl = impl[1]
  if(is.null(group)) {
    stopifnot(type=="onesample")
  }
  else {
    group <- factor(group)
    stopifnot(nlevels(group) == 2) 
    group.lev = levels(group)
  }
  if(impl == "wilcox.exact") {
    packages = "exactRankTests"
  }
  else if(impl == "wilcox.coin") {
    packages = "coin"
    data = cbind(data.frame(.group=group),mat)    
  }
  pvals = foreach(resp.var=colnames(mat),.combine=c,.packages=packages) %do% {
    if(impl == "wilcox.exact") {
      resp = mat[,resp.var]
      pval = wilcox.exact(x=resp[group==group.lev[2]],
                          y=resp[group==group.lev[1]],
                          paired=(type=="paired"))$p.value
    }
    else if(impl == "wilcox.coin") {
      form = as.formula(paste(resp.var,".group",sep="~"))
      
      if(type == "unpaired") {
        res.test = wilcox_test(form,
                               data=data)
      }
      else if(type == "paired") {
        res.test = wilcoxsign_test(form,
                                   data=data)
      }
      pval = pvalue(res.test)
    }
    pval
  }
  names(pvals) = colnames(mat)
  return (pvals)
}

wilcox.test.multi.fast <- function(mat,group=NULL,type="unpaired",pval=T,only.pval=T,pvalues.impl="genesel") {
  tr.mat=t(mat)
  if(is.null(group)) {
    stopifnot(type=="onesample")
    group = rep(1,ncol(tr.mat))
  }
  res = RankingWilcoxonAT(tr.mat,group,type=type,pvalues=pval,pvalues.impl="genesel")
  res = toplist(res,nrow(tr.mat))
  res[res$index,] = res
  if(only.pval) {
    return (res$pval)
  }
  else {
    return (res)
  }
}

#test.method argument is for passing function definition to the new
#parallel process on WIndows under Snow. Set it to wilcox.test.multi.fast
#when calling 'boot'
booted.wilcox.test.multi.fast <- function(n,mat,group,indices,test.method) {
  library(GeneSelector)
  ##replace=TRUE to select more samples than in original dataset
  ##TODO: apply strata to keep group count ratio
  ind.n = sample(indices, n, replace = TRUE, prob = NULL)
  d <- mat[ind.n,] # allows boot to select n samples
  g <- group[ind.n]
  #print(summary(g))
  test.method(d,g)
}


# run one test from boot indices
# pass exact=F because boot convertes warnings "cannot compute exact p-value with ties"
#into errors
booted.wilcox.test <- function(n,is.paired,exact,data,indices) {
  d <- data[indices[1:n],] # allows boot to select n sample
  with(d,wilcox.test(freq.x,freq.y,paired=is.paired,exact=exact))$p.value
}


booted.adonis.test.koren <- function(n,data,indices) {
  #nlev = length(levels(data$time))
  d <- data[indices,]
  if(!is.null(n)) {
    d <- sample.rows(d,n*length(levels(d$time)),replace=F) #allows boot to select n sample
  }
  l = count_matr_from_df(d,col_ignore=data_factors)
  ad.res = adonis(l~time,data=d,permutations=1000,method="bray")
  #print (ad.res)
  ad.res$aov.tab$"Pr(>F)"[1]
}

power.koren<-function() {
  
  data = read.koren("GG_OTU_table_L6.txt")
  data_factors = c("time","id_repl")
  #data = count.filter(data,col_ignore=data_factors,min_max_frac=0.1,min_row_sum=0)
  data_col = melt(data,id.vars=data_factors,variable.name="feature",value.name="freq")
  data_time = merge(data_col[data_col$time==1,],data_col[data_col$time==2,],by=c("id_repl","feature"))
  data_time_summary = plyr::ddply(data_time,c("feature"),summarize,
                                  mean.x=mean(freq.x), sd.x=sd(freq.x), n.x = length(freq.x),
                                  mean.y=mean(freq.y), sd.y=sd(freq.y), n.y = length(freq.y),
                                  mean_paired=mean(freq.x-freq.y),
                                  sd_paired=sd(freq.x-freq.y),n_paired=length(freq.x))
  data_time_summary$effect_t_paired = data_time_summary$mean_paired/data_time_summary$sd_paired
  #get pooled sd (original group ratios)
  data_time_summary$sd = with(data_time_summary,sqrt(((n.x-1)*sd.x**2+(n.y-1)*sd.y**2)/(n.x+n.y)))
  data_time_summary$effect_t = with(data_time_summary,(mean.x-mean.y)/sd)
  
  
  
  #feature_sel = "Root;k__Bacteria;p__Actinobacteria"
  alpha = 0.05
  power.cut = 0.8
  n = 50 # equal sample sizes per group
  eff.div = 1
  is.paired = T
  
  mv.test.res = booted.adonis.test.koren(NULL,data,seq(nrow(data)))
  
  print (paste("adonis p-value: ",mv.test.res))
  
  power.mv.res = mean(boot(data=data,statistic=booted.adonis.test.koren,R=100,n=n,strata=data$time)$t<=alpha)
  
  print(paste("adonis power: ",power.mv.res))
  
  n.tests = length(levels(data_time$feature))
  alpha.bn = alpha/n.tests
  feature = c()
  p.value = c()
  power = c()
  for (feature_sel in levels(data_time$feature)) {
    
    data_sel = data_time[data_time$feature==feature_sel,]
    
    #test.res = with(data_sel,t.test(freq.x,freq.y,paired=is.paired))
    test.res = with(data_sel,wilcox.test(freq.x,freq.y,paired=is.paired,exact=F))
    
    feature = c(feature,feature_sel)
    p.value = c(p.value,test.res$p.value)
    
    data_summ_sel = data_time_summary[data_time_summary$feature==feature_sel,]
    if(is.paired) {
      test.type = "paired"
    }
    else {
      test.type = "two.sample"
    }
    
    sd_power = with(data_summ_sel,sqrt((n-1)/n*(sd.x**2+sd.y**2)/2))
    #power.res = with(data_summ_sel,power.t.test(n=n,delta=(mean.x-mean.y)/eff.div,sd=sd_power,sig.level=alpha.bn,type=test.type))$power
    #power.res = with(data_summ_sel,mean(replicate(1000, wilcox.test(rnorm(n,mean.x,sd_power), rnorm(n,mean.y,sd_power), paired=is.paired,exact=F)$p.value)<=alpha.bn))
    power.res = mean(boot(data=data_sel,statistic=booted.wilcox.test,R=100,n=n,is.paired=is.paired,exact=F)$t<=alpha.bn)
    
    power = c(power,power.res)
  }
  
  p.value.bh = p.adjust(p.value,method="BH")
  p.value.bn = p.value*n.tests
  test.res = data.frame(feature,p.value,p.value.bh,p.value.bn,power)
  
  #we use BN correction to pick those features which we consider to be significant and
  #for which we measure the power
  test.res.sig = test.res[test.res$p.value.bn<=alpha,]
  test.res.sig$power.ok = test.res.sig$power >= power.cut
  
  test.res.sig = merge(test.res.sig,data_time_summary,by=c("feature"))
  
  #print (test.res.sig)
  
  print(paste("Mean power of the test: ",mean(test.res.sig$power)))
  
}

## compute Cramer V (also known as Cramer Phi) from the
## results of chisq.test or g.test of independence (on contingency table)
## of goodness of fit test
## see also function vcd::assocstats or lsr::cramersV
## For test of independence:
##V = sqrt(X^2 / [nobs * (min(ncols, nrows) - 1)])
## For goodnes of fit:
##V = sqrt(X^2 / [nobs * (nrows - 1)])
cramers.v <- function(x) {
  obs = x$observed
  k = NA
  if(length(dim(obs)) == 2 && all(dim(obs)>=2)) {
    k = min(nrow(obs),ncol(obs))
  }
  else {
    k = length(obs)
  }
  cv = sqrt(x$statistic /
              (sum(obs) * (k - 1)))
  return(as.numeric(cv))
}


Xmcupo.effectsize.par <-
  function(par.groups,reads.groups){
    par.groups <- par.groups
    assert.non.zero.dm.par(par.groups)    
    for(i.group in seq(length(par.groups))) {
      par.groups[[i.group]]$reads = reads.groups[[i.group]]
    }
    
    N.group <- length(par.groups)
    Kc <- length(par.groups[[1]]$pi)
    N.total <- do.call(sum, lapply(par.groups, function(x){sum(x$reads)}))
    
    group.parameter.estimated <- lapply(par.groups, function(x){c(x$pi, x$theta, x$reads)})
    Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K=Kc)
    
    if(Kc >= N.group){
      pi.groups <- diag(rep(1, N.group))
    }else{
      stop("The number of taxa must be greater than the number of groups.")
    }
    
    gp.max <- lapply(as.list(1:N.group), function(x,gp.1, pi.groups){
      gpp <- gp.1[[x]]
      ret <- c(pi.groups[x,], gpp$theta, gpp$reads)
      return(ret)
    }, gp.1=par.groups, pi.groups=pi.groups)
    
    Xmcupo.gp <- Xmcupo.statistics(gp.max, K=N.group)
    CramerV<- sqrt(Xmcupo[1]/(N.total*min(Kc-1, N.group-1)))
    Mod.CramerV <- sqrt(Xmcupo[1]/(Xmcupo.gp[1]*min(Kc-1, N.group-1)))
    
    result <- c(Xmcupo[1],CramerV, Mod.CramerV)
    names(result) <- c("Chi-Squared", "Cramer Phi", "Modified-Cramer Phi")
    return(result)
  }

Xmcupo.sevsample.par <-
  function(par.groups, reads.groups, K){
    if(missing(par.groups) || missing(reads.groups) || missing(K))
      stop("par.groups, reads.groups or K is missing.")
    assert.non.zero.dm.par(par.groups)
    n.groups <- length(par.groups)
    index <- as.matrix(seq(1:n.groups))
    group.parameter.estimated <- list()
    
    for(x in index){
      par <- par.groups[[x]]
      nreads.data <- as.matrix(reads.groups[[x]])
      
      group.parameter.estimated[[x]] <- c(par$pi, par$theta, t(nreads.data))
    }
    
    Xmcupo <- Xmcupo.statistics(group.parameter.estimated, K)
    p.value <- 1-pchisq(q=Xmcupo, df=(n.groups-1)*(K-1), ncp=0, lower.tail=TRUE)
    
    sevRAD.mean.test.upo <- list(Xmcupo, p.value)
    names(sevRAD.mean.test.upo) <- c("Xmcupo statistics", "p value")
    
    return(sevRAD.mean.test.upo)
  }


read.nistal<-function(file_name) {
  data = read.csv(file=file_name,head=TRUE,sep="\t",row.names="id")
  data$group = as.factor(data$group)  
  data
}

power.nistal<-function() {
  
  data_all = read.nistal("children.txt")
  groups = as.factor(c("Celiac","Control"))
  attr_names = c("group")
  all_features = colnames(data_all)[!(colnames(data_all) %in% attr_names)]
  test_features = as.factor(c("Prevotella.sp.","Streptococcus.sp."))
  #test_features = as.factor(c("Comamonas.sp."))
  
  n.samp.grp = c(180,60)
  n.seq.range = seq(1500,2000,by=1)
  #n.seq.range <- rep(50,100)
  #n.samp.grp = c(8,5)
  
  dm.par.orig = list()
  
  #set number of sequences here so that we can estimate
  #effect size and DM test statistics outside of sample generation loop
  n.seq = list()
  
  for (i.group in seq(length(groups))) {
    group = groups[i.group]
    data = data_all[data_all$group == group,]
    cnt_m = count_matr_from_df(data,col_ignore=attr_names)
    dm.par.orig[[i.group]] <- DM.MoM(cnt_m)
    n.seq[[i.group]] <- sample(n.seq.range, size=n.samp.grp[i.group]) 
  }
  
  res.power = power.dirmult.range(
    dm.par.orig=dm.par.orig,
    n.seq=n.seq,
    groups=groups,
    all_features=all_features,
    test_features=test_features,
    n.samp.grp=n.samp.grp
  )    
  write.csv(res.power,"res.power.csv")
}

dirmult.comp.gamma <- function(pi,theta) {
  # standard formula for gamma
  pi * (1-theta)/theta
}

assert.non.zero.dm.par <- function(dm.par) {
  all.pi = plyr::laply(dm.par,function(l){l$pi})
  stopifnot(all(colSums(all.pi)>0))
}

batch.jobs.file.dir <- function(id) {
  tempfile(pattern = paste(id,".",sep=""), tmpdir = getwd(), fileext = ".batch_reg")
}

batch.jobs.chunked.ids <- function(reg,n.chunks) {
  chunk(getJobIds(reg), n.chunks=n.chunks, shuffle=TRUE)
}

power.dirmult.range<-function(
  dm.par.orig,
  n.seq,
  groups,
  all_features,
  test_features,
  n.samp.grp,
  alpha=0.05,
  effect.size.range = c(0, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0),
  n.rep = 100,
  n.adonis.perm = 400
  
) {
  
  assert.non.zero.dm.par(dm.par.orig)
  
  n_features = length(all_features)
  
  test.ad.res.power.eta = c()
  r2.ad.eta = c()
  test.wil.res.power.sel.eta = c()
  mod.cramer.phi = c()
  test.Xmcupo.res = c()
  
  reg = makeRegistry("onesamp",
                     file.dir=batch.jobs.file.dir("onesamp"),
                     packages=packages,
                     skip=F)
  
  par.list = foreach (eta = effect.size.range) %do% {
    # copying the DM parameters of the alternative hypothesis from one of the parent datasets.
    dm.par <- dm.par.orig
    # interpolate the proportion vectors of the two parent datasets
    dm.par[[2]]$pi <- dm.par[[1]]$pi * (1-eta) + dm.par.orig[[2]]$pi  *(eta)
    dm.par[[2]]$theta <- dm.par[[1]]$theta * (1-eta) + dm.par.orig[[2]]$theta  *(eta)  
    # standard formula for gamma
    dm.par[[2]]$gamma <- dirmult.comp.gamma(dm.par[[2]]$pi,dm.par[[2]]$theta)
    
    mod.cramer.phi = Xmcupo.effectsize.par(dm.par,n.seq)["Modified-Cramer Phi"]
    test.Xmcupo.res = Xmcupo.sevsample.par(dm.par,n.seq,n_features)$"p value"
    
    list(dm.par=dm.par,eta=eta,mod.cramer.phi=mod.cramer.phi,test.Xmcupo.res=test.Xmcupo.res)
  }    
  one.sample.action = function(par,i.rep,all_features,groups,n.seq,n.samp.grp,n.adonis.perm) {
    dm.par = par$dm.par
    ### Generate a random vector of number of reads per sample
    dm.counts = matrix(nrow=0,ncol=length(all_features),dimnames=list(NULL,all_features))
    dm.group = c()
    for (i.group in seq(length(groups))) {
      dm.counts <- rbind(dm.counts,Dirichlet.multinomial(n.seq[[i.group]], dm.par[[i.group]]$gamma))
      dm.group <- c(dm.group,rep(groups[i.group],n.samp.grp[i.group]))
    }
    dm.freq = norm.prop(dm.counts)
    dm.group = as.factor(dm.group)
    ad.res = adonis(dm.freq~dm.group,permutations=n.adonis.perm,method="bray")
    #print (ad.res)
    
    test.wil.res = foreach (feature = all_features,.combine=c) %do% {
      #should use groups[1] and groups[2], but rep(factor,...) loses labels
      dm.freq.1 = dm.freq[dm.group==1,feature]
      dm.freq.2 = dm.freq[dm.group==2,feature]
      wilcox.test(dm.freq.1,dm.freq.2,paired=FALSE,exact=F)$p.value
    }
    
    return(list(ad.res=ad.res,test.wil.res=test.wil.res,
                eta=par$eta,
                mod.cramer.phi=par$mod.cramer.phi,
                test.Xmcupo.res=par$test.Xmcupo.res
    ))
  }
  
  batchExpandGrid(reg,one.sample.action,
                  par=par.list,i.rep=seq(n.rep),
                  more.args=list(all_features,groups,n.seq,n.samp.grp,n.adonis.perm))
  
  submitJobs(reg, ids=batch.jobs.chunked.ids(reg,n.chunks=10))
  
  #reduceResultsDataFrame loses column names, so we use plyr::ldply on a list:
  res.all = plyr::ldply(reduceResultsList(reg,fun=function(job,res,all_features,test_features,alpha) {
    test.ad.res = res$ad.res$aov.tab$"Pr(>F)"[1]
    r2.ad = res$ad.res$aov.tab$R2[1]
    test.wil.res = res$test.wil.res
    names(test.wil.res) = all_features
    test.wil.res = t(test.wil.res)
    test.ad.res.power = mean(test.ad.res<=alpha)
    for (i.col in seq(ncol(test.wil.res))) {
      test.wil.res[,i.col] = p.adjust(test.wil.res[,i.col],method="BH")
    }
    test.wil.res.power = rowMeans(test.wil.res<=alpha)
    test.wil.res.power.sel = test.wil.res.power[all_features %in% test_features]
    cbind(eta=res$eta,mod.cramer.phi=res$mod.cramer.phi,r2.ad,test.Xmcupo.res=res$test.Xmcupo.res,test.ad.res,test.wil.res)
  },
  test_features=test_features,
  all_features=all_features,
  alpha=alpha)
  )
  res = plyr::ddply(res.all,c("eta"),
                    function(x,alpha) { 
                      first = 
                        summarize(x[,-which(names(x) %in% all_features)],
                                  mod.cramer.phi = mean(mod.cramer.phi),
                                  test.Xmcupo.res = mean(test.Xmcupo.res),
                                  r2.ad = mean(r2.ad),
                                  test.ad.res.power = mean(test.ad.res<=alpha))
                      #use dfrm[col_ind] notation because dfrm[,col_ind] will
                      #return a vector if length(col_ind) == 1
                      second = 
                        as.data.frame(colMeans(x[which(names(x) %in% test_features)]<=alpha))
                      cbind(first,t(second))
                    },
                    alpha = alpha
  )
  res = res[order(res$eta),]
  col.names.features = plyr::laply(names(res)[names(res) %in% test_features],function(x){paste(x,"Mann-Whitney U power")})
  colnames(res) = c(c("Parameter scaling coeff.", "Modified Cramer Phi", 
                      "Gen. Wald-type p-value", "Adonis R2", "Adonis power"),
                    col.names.features)
  
  return(t(res))
}


sorted.freq.kelvin<-function() {
  group_sel = c("Stool","Saliva")
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls")
  
  data = data_all[data_all$group %in% group_sel,]
  data = count.filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
  
  freq.mean = sorted.freq(data,col_ignore=attr_names)
  
  write.csv(as.matrix(sort(freq.mean[,"Saliva"],decreasing=TRUE)),"freq.top.saliva.csv")
  
  return(freq.mean)
}

sorted.freq.nistal<-function() {
  
  attr_names = c("group")
  data_all = read.nistal("children.txt")
  data = data_all
  #data = count.filter(data,col_ignore=attr_names,min_max_frac=0.2,min_row_sum=2000)  
  #return(data_all)
  
  freq.mean = sorted.freq(data,col_ignore=attr_names)
  
  write.csv(as.matrix(sort(freq.mean[,"Control"],decreasing=TRUE)),"freq.top.control.csv")
  
  return(freq.mean)
}

#freq.mean = sorted.freq.kelvin()
#freq.mean = sorted.freq.nistal()
#res = power.nistal()
#x = res[[1]]
#y = res[[2]]
#z = res[[3]]

read.humann.summary<-function(file_name) {
  #skip first 4 rows with diversity indices
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  data = data[5:nrow(data),]
  data$ID = NULL
  row.names(data) = data$NAME
  data$NAME = NULL
  data = t(data)
  row.names(data) = unlist(lapply(strsplit(row.names(data), "\\."),function(x) x[1]))
  return (as.data.frame(data))
}

annot.levels <- function(annot.type) {
  if( annot.type == "humann" ) {
    return (c("level.3"))
  }
  levs = switch(annot.type,
                cog=c("level.2"),
                c("level.2","level.3")
  )
  levs = c(levs,c("function.","gene"))
  return (levs)
}

read.mgrast.summary<-function(file_name,file_name.id.map=NULL) {
  #if quote!=NULL, resulting dataset is cropped in mid-file with a warning about
  #EOF inside a quoted string. Quotes probably do not match in pairs.
  data = read.delim(file_name, header=T,stringsAsFactors=T,quote=NULL)
  #there are more unique id values than function values, which means
  #a mapping id -> function is 1->many. We create field "gene" that gives
  #a descriptive name to the id
  data$gene = as.factor(paste(data$function.,data$id,sep="."))
  if (!is.null(file_name.id.map)) {
    id.map = read.delim(file_name.id.map, header=T,stringsAsFactors=T,row.names="mgrast_metagenome_id")
  }
  return (merge(data,id.map,by.x="metagenome",by.y="row.names"))
}

## Somehow ggplot silently drops labels that have "/" in them; spaces get converted to dots
## in other situations. This replaces problematic symbols. This has to be applied uniformly
## to all files that are loaded and contain taxonomy names. Of course, if you write it out,
## you will get a mismatch in any downstream analysis by other tools. This is therefore a hack.
## Use with care!
sanitize.taxa.names <- function(x) {
  gsub("[/ ]","_",x)
}

read.mothur.otu.shared <- function(file_name,sanitize=T) {
  require(data.table)
  #data = read.delim(file_name, header=T,stringsAsFactors=T)
  #when read.delim is used, X == NA column comes because there is an extra delimiter at the end of line
  #data$X = NULL
  data = fread(file_name, header=T, sep="\t", stringsAsFactors=T, data.table=F)
  last.col = ncol(data)
  if(all(is.na(data[,last.col]))) {
    data[,last.col] = NULL
  }
  data$label = NULL
  data$numOtus = NULL
  row.names(data) = data$Group
  data$Group = NULL
  if(sanitize) {
    names(data) = sanitize.taxa.names(names(data))
  }
  return (as.matrix(data))
}

read.mothur.cons.taxonomy <- function(file_name,sanitize=T,taxa.level="otu",taxa.levels.mix=0) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  row.names(data) = data$OTU
  data$Taxa = plyr::laply(strsplit(as.character(data$Taxonomy),"\\([0-9]*\\);"),
                          function(x,taxa.level) {
                            if(!is.taxa.level.otu(taxa.level)) {
                              ##only look that deep in lineage
                              x = x[1:as.integer(taxa.level)]
                            }
                            ##pmatch returns the index of the first match, so the expression below
                            ##returns the prefix of lineage till (excluding) the first "unclassified" element
                            before_tail = pmatch("unclassified",x,nomatch=length(x)+1,dup=T)-1
                            if(before_tail>0) {
                              no.tail = x[1:before_tail]
                              last = no.tail[length(no.tail)]
                              if(length(no.tail) < (length(x) - taxa.levels.mix)) {
                                last = paste("Unclassified",last,sep="_")
                              }
                            }
                            else {
                              last = "Unclassified"
                            }
                            last
                          },
                          taxa.level=taxa.level
  )
  if(sanitize) {
    data$Taxa = sanitize.taxa.names(data$Taxa)
  }
  return (data)
}

is.taxa.level.otu <- function(taxa.level) {
  as.character(taxa.level) == "otu"
}

read.mothur.otu.with.taxa <- function(otu.shared.file,cons.taxonomy.file,sanitize=T,taxa.level="otu",
                                      count.basis="seq", otu.count.filter.options=NULL,taxa.levels.mix=0) {
  otu.df = read.mothur.otu.shared(otu.shared.file,sanitize=sanitize)
  taxa.df = read.mothur.cons.taxonomy(cons.taxonomy.file,sanitize=sanitize,taxa.level=taxa.level,taxa.levels.mix=taxa.levels.mix)
  ##make.global(name="otu.fr")
  stopifnot(ncol(otu.df) == nrow(taxa.df))
  stopifnot(all(colnames(otu.df) %in% taxa.df$OTU))
  otu.df = new.m_a(count=otu.df)
  otu.df = report.count.filter.m_a(otu.df,count.filter.options=otu.count.filter.options)$count
  otu.name.ind = match(colnames(otu.df),taxa.df$OTU)
  ## this will drop "other" category
  otu.name.mask = !is.na(otu.name.ind)
  otu.name.ind = otu.name.ind[otu.name.mask]
  taxa.df = taxa.df[otu.name.ind,]
  otu.df = otu.df[,otu.name.mask]
  stopifnot(all(colnames(otu.df) == taxa.df$OTU))
  if(count.basis=="otu") {
    otu.df = ifelse(otu.df > 0,1,0)
  }
  if(!is.taxa.level.otu(taxa.level)) {
    x = aggregate(t(otu.df),list(row.ids=taxa.df$Taxa),sum)
    rownames(x) = x$row.ids
    x$row.ids = NULL
    otu.df = t(x)
  }
  else {
    colnames(otu.df) = paste(taxa.df$Taxa,taxa.df$OTU,sep=".")
  }
  ## Order columns by taxa name
  ##otu.df = otu.df[,order(colnames(otu.df))]
  return (otu.df)
}

read.mothur.otu.with.taxa.m_a <- function(...) {
  count = read.mothur.otu.with.taxa(...)
  return(new.m_a(count=count))
}

make.mothur.taxa.summary.feature.names <- function(taxa.summary) {
  taxon = taxa.summary$taxon
  names(taxon) = taxa.summary$rankID
  ind_unclass = which(taxon == "unclassified")
  ind_unclass.ini = ind_unclass
  while(length(ind_unclass)>0) {
    taxon[ind_unclass] = taxon[unlist(strsplit(names(taxon)[ind_unclass], "\\.[0-9]+$"))]
    ind_unclass = ind_unclass[taxon[ind_unclass]=="unclassified"]
  }
  ind_unclass.ini = ind_unclass.ini[taxon[ind_unclass.ini] != "unknown"]
  taxon = as.character(taxon)
  taxon[ind_unclass.ini] = paste("Unclassified",taxon[ind_unclass.ini],sep="_")
  return(taxon)
}

read.mothur.taxa.summary <- function(file_name,sanitize=T) {
  data = read.delim(file_name, header=T,stringsAsFactors=T)
  #X == NA column comes because there is an extra delimiter at the end of line
  data$X = NULL
  if(sanitize) {
    data$taxon = sanitize.taxa.names(data$taxon)
  }
  data$feature = as.factor(make.mothur.taxa.summary.feature.names(data))
  return (data)
}


multi.mothur.to.abund.m_a <- function(data,level) {
  data.level = data[data$taxlevel==level,,drop=F]
  attr = c("taxlevel","rankID","taxon","daughterlevels","total","feature")
  x = split_count_df(data.level,col_ignore=attr)
  row.names(x$count) = x$attr$feature
  x$count = t(x$count)
  x$attr.feat = x$attr
  x$attr = data.frame(SampleID=rownames(x$count))
  return (x)
}

multi.mothur.to.abund.df <- function(data,level) {
  x = multi.mothur.to.abund.m_a(data,level)
  return (as.data.frame(x$count))
}

mgrast.to.abund.df <- function(data,level) {
  x = data[,c("project_sample_id",level,"abundance")]
  x = plyr::ddply(x,c("project_sample_id",level),summarize,abundance=sum(abundance))
  x = dcast(x,list("project_sample_id",level),value.var="abundance",fill=0.)
  row.names(x) = x$project_sample_id
  x$project_sample_id = NULL
  return (x)
}


#Merge counts data frame with meta data frame.
#Merge in this method is always done on row.names. If you need to
#merge on a component of row.names(counts), first
#get row.names(counts), split them into components, merge with 
#meta data frame, and then use the resulting frame as a new meta data
##TODO: remove intermediate conversion into data.frame
merge.counts.with.meta <- function(x,y,suffixes=c("","meta")) {
  ## if metadata is not provided, just make a very simple one right here
  if(is.null(y)) {
    y = data.frame(SampleID=rownames(x))
    rownames(y) = y$SampleID
    y$Group = y$SampleID
  }
  mrg = merge(x,y,by="row.names",suffixes=suffixes)
  #Row.names column is generated by merge() when by="row.names"
  #the assignment below also serves as assertion that count records were
  #not duplicated as a result of merge
  row.names(mrg) = mrg$Row.names
  mrg$Row.names = NULL
  all.names = colnames(mrg)
  attr.names = all.names[!(all.names %in% colnames(x))]
  ret = split_count_df(mrg,col_ignore=attr.names)
  return (ret)
}

aggregate.by.meta.data.m_a <- function(m_a,
                                       group_col,
                                       count_aggr=sum,
                                       attr_aggr=NULL,
                                       group_col_result_name="SampleID",
                                       colwise=T) {
  x = m_a
  
  groups = list()
  groups[[group_col_result_name]] = x$attr[,group_col]
  
  #We need to drop resulting grouping name from the input,
  #otherwise we will get duplicated names (e.g. two SampleID columns. 
  #The assumption here is that if such name is supplied for the grouping
  #output column, then the existing column with the same name is not needed.
  if(group_col_result_name %in% names(x$attr)) {
    x$attr = drop.columns(x$attr,c(group_col_result_name))
  }
  ##TODO: make the default attr_aggr to drop all columns
  ## with more that one value in any group. 
  if(is.null(attr_aggr)) {
    attr_aggr = function(y) {y[1]}
  }
  x$attr = aggregate(x$attr,groups,attr_aggr)
  row.names(x$attr) = x$attr[,group_col_result_name]
  x$attr = droplevels(x$attr)
  if(!is.null(x$count)) {
    count_names = colnames(x$count)
    if(colwise) {
      x$count = aggregate(x$count,groups,count_aggr)  
      row.names(x$count) = x$count[,group_col_result_name]
      x$count = drop.columns(x$count,c(group_col_result_name))
    }
    else {
      dt = data.table(count=x$count,group=groups[[group_col_result_name]])
      dt = dt[,data.table(do.call(count_aggr,list(.SD))),by=group]
      rn = dt[,group]
      dt[,group:=NULL]
      setnames(dt,count_names)
      x$count = setDF(dt)
      row.names(x$count) = rn
    }
    return (merge.counts.with.meta(x$count,x$attr))
  }
  else {
    return (x)
  }
  
}


aggregate.by.meta.data <- function(meta_data,
                                   col_ignore=c(),
                                   ...) {
  m_a = split_count_df(meta_data,col_ignore=col_ignore)
  m_a = aggregate.by.meta.data.m_a(m_a,...)
  return (join_count_df(m_a))
}




melt.abund.meta <- function(data,id.vars,attr.names,value.name="abundance") {
  feature.names = get.feature.names(data,attr.names)
  data$.record.id=rownames(data)
  if(!is.null(data$.record.id)) {
    id.vars = c(".record.id",id.vars)
  }
  return (reshape2::melt(data,id.vars=id.vars,measure.vars=feature.names,variable.name="feature",value.name=value.name))
}

sort.factor.by.total <- function(factor.val,sort.val,ordered=F,decreasing=TRUE) {
  o = aggregate(sort.val, list(factor.val), sum)
  o = o[order(o$x, decreasing=decreasing),1]  
  
  return (factor(factor.val, levels = o, ordered = ordered))
}

order.levels <- function(lev,keys) {
  lev[order(keys,decreasing=T)]
}

#' split columns names with factors into two subsets to have about equal number of levels
#' for each dimension in ggplot facet grid plot
split.by.total.levels.data.frame <- function(x) {
  ##transpose and column subsetting below is to get first row as a vector,
  ##otherwise the input to cumprod is a data.frame with a single row, and cumprod then works column-wise,
  ##doing nothing
  cuprod = cumprod(t(plyr::colwise(function(y) length(levels(factor(y))))(x))[,1])
  cuprod.max = cuprod[length(cuprod)]
  ind.split = which(cuprod >= (cuprod.max/2))[1]
  colnam = colnames(x)
  if(ind.split==1) {
    if(length(colnam)<2) {
      first_part = c()
    }
    else {
      first_part = colnam[1]
    }
  }
  else{
    first_part = colnam[1:(ind.split-1)]
  }
  second_part = colnam[!colnam %in% first_part]
  return(list(first_part,second_part))
}

ggplot.hue.colors <- function(n,l=60,c=200) {
  hues = seq(15, 375, length=n+1)
  grDevices::hcl(h=hues, l=l, c=c)[1:n]
}

#' copied from Heatplus RainbowPastel
mgsat.rainbow.pastel <- function (n, blanche = 200, ...) 
{
  cv = rainbow(n, ...)
  rgbcv = col2rgb(cv)
  rgbcv = pmin(rgbcv + blanche, 255)
  rgb(rgbcv[1, ], rgbcv[2, ], rgbcv[3, ], maxColorValue = 255)
}

#' exclude.qual.colors - these colors will be dropped from qualitative palettes. Current
#' default drops yellow because it yellow points get washed out on screens
#' Other suggested exclusions for Brewer palettes:
#' Dark2 -> "#66A61E" (second green color)
generate.colors.mgsat <- function(x,value=c("colors","palette"),family=c("brewer","ggplot"),
                                  brewer.pal.name="Set1",
                                  exclude.qual.colors=c("#FFFF33")) {
  family = family[[1]]
  if(family=="brewer") library(RColorBrewer)
  else if(family=="ggplot") library(ggplot2)
  
  value = value[[1]]
  if(is.character(x)) {
    x = factor(x)
  }
  else if(is.logical(x)) {
    x = factor(x)
  }
  if(is.numeric(x)) {
    require(lattice)
    at = do.breaks(range(x), 30)
    ret = level.colors(x, at = at,col.regions = topo.colors)
  }
  else if(is.factor(x)) {
    lev = levels(x)
    if(!is.ordered(x)) {
      if(family=="brewer") {
        pal.info = brewer.pal.info[brewer.pal.name,]
        n.color.orig = pal.info$maxcolors
        n.excl = length(exclude.qual.colors)
        n.request = n.excl + length(lev)
        if(n.request>n.color.orig) {
          n.request = n.color.orig
        }
        palette = brewer.pal(n.request, brewer.pal.name)
      }
      else if(family=="ggplot") {
        palette = ggplot.hue.colors(n.request) 
      }
      palette = palette[!(palette %in% exclude.qual.colors)]
      #get.palette = colorRampPalette(palette)
      #palette = get.palette(max(length(features),n.color.orig))
      palette = rep_len(palette,length(lev))
    }
    else {
      palette = ggplot.hue.colors(length(lev))
    }
    names(palette) = lev
    if(value=="colors") {
      ret = palette[x]
      names(ret) = NULL
    }
    else {
      ret = palette
    }
  }
  return(ret)
}


pairs.scatter.plot <- function(m_a,group=NULL,tooltip=NULL,color=NULL,...) {
  library(pairsD3)
  count = m_a$count
  attr = m_a$attr
  if(!is.null(group)) {
    group = eval(parse(text=group),attr)
  }
  if(!is.null(tooltip)) {
    tooltip = eval(parse(text=tooltip),attr)
  }
  if(is.null(color) && !is.null(group)) {
    color = generate.colors.mgsat(group,value="palette",brewer.pal.name = "Set1")
  }
  pairsD3(x=count,
          group=group,
          tooltip=tooltip,
          big=T,
          col=color,
          ...)
}

pairs.scatter.plot.rbokeh <- function(m_a,
                                      width.cell=200,
                                      height.cell=200,
                                      tools=c("pan", "wheel_zoom", "box_zoom", "box_select", "reset"),
                                      lod_threshold=10) {
  require(rbokeh)
  var.names = colnames(m_a$count)
  data = join_count_df(m_a)
  nms <- expand.grid(var.names, rev(var.names), stringsAsFactors = FALSE)
  n.nms = length(var.names)
  nms$yaxis <- rep(c(TRUE, rep(FALSE, n.nms-1)), n.nms)
  nms$xaxis <- c(rep(FALSE, (n.nms-1)*n.nms), rep(TRUE, n.nms))
  #nms$h <- height.cell
  #nms$w <- width.cell
  #nms$h[nms$xaxis] <- height.cell
  #nms$w[nms$yaxis] <- width.cell
  splom_list <- vector("list", n.nms*n.nms)
  color = "Species"
  for(i in seq_along(splom_list)) {
    
    if(i==1) legend = T
    else legend = F
    
    pl = figure(width = width.cell, height = height.cell,
                tools = tools, min_border = 2, lod_threshold = lod_threshold)
    pl = eval(parse(text=sprintf('ly_points(pl,%s, %s, data = data,
              color = %s, hover = Species, size = Petal.Length*2, legend = %s)',nms$Var1[i],nms$Var2[i],color,legend)))
    pl = pl %>% x_axis(visible = nms$xaxis[i]) %>% y_axis(visible = nms$yaxis[i])
    pl = pl %>% theme_axis(c("x", "y"),major_label_orientation = "vertical")
    splom_list[[i]] = pl
  }
  
  grid_plot(splom_list, nrow = n.nms, ncol = n.nms, same_axes = c(T,T), link_data = TRUE, simplify_axes = T) %>% 
    tool_save() %>% tool_resize()
  #saveWidget(paired.scatter.plot(split_count_df(iris,col_ignore = "Species")),"tmp.html",selfcontained = F)
}

dev.editable.svg <- function(file, width=800, height=600,bg="white",...) {
  ## this graphics device saves text as text, unlike the default grDevices::svg that
  ## converts text into curves. InkScape or Adobe tools can be used to edit SVG files.
  library(RSvgDevice)
  devSVG(file = file, width = width/72, height = height/72, bg = bg,...)
}


require(scales) # trans_new() is in the scales library
signed_sqrt_trans = function() trans_new("signed_sqrt", function(x) sign(x)*sqrt(abs(x)), function(x) sign(x)*sqrt(abs(x)))

plot.abund.meta <- function(m_a,
                            ci=NULL,
                            id.vars=c(),
                            value.name="Abundance",
                            file_name=NULL,
                            ggp.comp=NULL,
                            facet_grid.margins=FALSE,
                            features.order=NULL,
                            feature.names.conv=NULL,
                            geom="bar",
                            line.show.points=T,
                            line.legend.thickness=rel(3),
                            n.top=20,
                            id.var.dodge=NULL,
                            flip.coords=T,
                            sqrt.scale=F,
                            stat_summary.fun.y="mean",
                            make.summary.table=T,
                            facet_wrap_ncol=3,
                            legend.title=NULL,
                            record.label=NULL,
                            theme_font_size = 0.8,
                            hide.ticks.x=F,
                            hide.axis.text.x=F,
                            show.samp.n=T,
                            axis.text.rel=1.,
                            axis.text.rel.x=axis.text.rel,
                            axis.text.rel.y=axis.text.rel,
                            dodge.width=0.7) {
  
  if(is.null(id.var.dodge)) {
    id.vars.facet = id.vars
  }
  else {
    id.vars.facet = id.vars[id.vars != id.var.dodge]
  }
  if(flip.coords && sqrt.scale) {
    stop("SQRT coordinate transformation cannot be used in flipped horizontal and vertical coords")
  }
  
  #TODO: change code in this method to use m_a directly
  data = join_count_df(m_a)
  attr.names = names(m_a$attr)
  dat = melt.abund.meta(data,id.vars=c(id.vars,record.label),attr.names=attr.names,value.name=value.name)
  if(!is.null(ci)) {
    ci.m = dcast(melt(ci,varnames=c(".record.id","feature","var")),.record.id+feature~var)
    n.rec.dat = nrow(dat)
    dat = plyr::join(dat,ci.m,by=c(".record.id","feature"),type="left",match="first")
    stopifnot(nrow(dat)==n.rec.dat)
  }
  
  if (is.null(features.order)) {
    dat$feature = sort.factor.by.total(dat$feature,dat[[value.name]])
  }
  else {
    dat$feature = factor(dat$feature,levels=features.order,ordered=F)
  }
  
  rownames.sorted = rownames(m_a$count)[
    do.call(order, -as.data.frame(m_a$count)[,levels(dat$feature),drop=F])
    ]
  dat$.record.id = factor(dat$.record.id,levels=rownames.sorted)
  
  if(make.summary.table) {
    dat.summary = eval(parse(text=sprintf('plyr::ddply(dat, c("feature",id.vars),
        plyr::summarise,
        mean = mean(%s,na.rm=T),
        sd = sd(%s,na.rm=T),
        median = median(%s,na.rm=T),
        incidence = mean(%s>0,na.rm=T)
        )',value.name,value.name,value.name,value.name)
    ))
    dat.melt = dat
  }
  else {
    dat.summary = NULL
    dat.melt = NULL
  }
  
  features = levels(dat$feature)
  if(!is.null(feature.names.conv)) {
    dat$feature = do.call(feature.names.conv,list(dat$feature))
  }
  if(is.factor(dat$feature)) {
    ##show only n.top
    features = levels(dat$feature)
    features = features[1:min(length(features),n.top)]
    dat = dat[dat$feature %in% features,,drop=F]
  }
  if(geom == "violin") {
    ##violin will fail the entire facet if one feature has zero variance,
    ##so we perturb data a tiny bit
    val = dat[,value.name]
    sd_dodge = max(abs(val),na.rm = T)*1e-6
    dat[,value.name] = val + rnorm(length(val),0,sd_dodge)
  }
  
  if(is.null(id.var.dodge)) {
    fill="feature"
    color="feature"
    palette.levels = factor(features,levels=features)
  }
  else {
    fill=id.var.dodge
    color=id.var.dodge
    palette.levels = factor(unique(dat[,id.var.dodge]))
  }
  
  make.facet.formula <- function(df,vars) {
    if(length(vars) == 0) {
      wr = NULL
      n.facet = 1
    }
    else if (length(vars) == 1) {
      wr = as.formula(paste("~",vars[1],sep=""))
      n.facet = num.levels(df[[vars[1]]])
    }
    else {
      splits = split.by.total.levels.data.frame(df[,vars,drop=F])
      wr = as.formula(paste(paste(splits[[2]],collapse = "+"),paste(splits[[1]],collapse = "+"),sep="~"))
      n.facet = c(num.levels.data.frame(df[,splits[[2]],drop=F]),
                  num.levels.data.frame(df[,splits[[1]],drop=F]))
    }
    return(list(facet.form=wr,n.facet=n.facet))
  }
  
  hide.panel.grid.major.x = F
  hide.panel.grid.major.y = F
  
  fc.form.res = make.facet.formula(dat,id.vars.facet)
  facet.form = fc.form.res$facet.form
  n.facet = fc.form.res$n.facet
  
  if(geom == "bar_stacked") {
    sqrt.scale = F
    dat$feature = factor(dat$feature,levels=rev(levels(dat$feature)))
    if(!is.null(record.label)) {
      x = record.label
    }
    else {
      x = ".record.id"
    }
    aes_s = aes_string(x=x,y=value.name,
                       fill = fill,color = color)
    
    gp = ggplot(dat, aes_s)
    
    gp = gp + #geom_col(position="stack",stat="summary",fun.y=stat_summary.fun.y) 
      stat_summary(fun.y=stat_summary.fun.y, geom="bar", 
                   position="stack")
    if(length(id.vars.facet) == 0) {
      wr = facet_null()
    }
    else if (length(id.vars.facet) == 1) {
      wr = facet_grid(facet.form,
                      drop=T,
                      scales="free_x", space = "free_x")
    }
    else {
      wr = facet_grid(facet.form,
                      drop=T,margins=facet_grid.margins,
                      scales="free_x", space = "free_x")
    }
    gp = gp + wr
    legend.position = "right"
  }
  else {
    ## here is what is going on with scale transformations in ggplot2 (v.1.0.0):
    ## scale_y_sqrt() - transforms data points before anything else is done like
    ## stats calculation or range detection; scale ticks are distributed quadratically
    ## (unevenly); this one can be combined with coord_flip().
    ## A (serious) downside of this method is that it can change relative magnitude
    ## of feature mean values as computed and shown by stat_summary(fun.y=mean) because large
    ## outliers will be compressed stronger.
    ## coord_trans(y = "sqrt") transforms only final geometric representation (including
    ## all glyphs) and creates a non-linear tick marks as well; it cannot be combined with
    ## coord_flip().
    ## Because of the above, we switch between two representations: vertical with
    ## coord_trans(y = "sqrt") and horizontal (flipped) with original linear coords
    ## 2016-10-05 - it looks like using 'group' breaks things now
    ## ,group = if(length(id.var.dodge) > 0) id.var.dodge else NULL
    aes_s = aes_string(x="feature",y=value.name,
                       fill = fill,color = color)
    gp = ggplot(dat, aes_s)
    
    if(geom == "bar") {
      # Create data.frame with shading info
      pos_dod = position_dodge(width=dodge.width) #0.9
      gp = gp + stat_summary(fun.y=stat_summary.fun.y, geom="bar", aes(width=0.5), 
                             position=pos_dod)
      #geom_obj = stat_summary(aes(label=round(..y..,2)), fun.y=mean, geom="text")
      #geom_obj = geom_bar(stat=stat_summary(fun.y="mean"),width=0.4)
      if(!is.null(ci) && stat_summary.fun.y=="identity") {
        gp = gp + geom_errorbar(position = pos_dod,
                                aes(ymin=lwr.ci, ymax=upr.ci,width=0.5),
                                color="black")
      }
    }
    else if(geom == "violin") {
      gp = gp + geom_violin(scale= "width", trim=TRUE, adjust=1)
    }
    else if(geom == "boxplot") {
      gp = gp + geom_boxplot(fill=NA,na.rm=T,notch=F)
    }
    else if(geom == "dotplot") {
      gp = gp + geom_dotplot(binaxis = "y", stackdir = "center", binpositions = "all",	
                             method="histodot")
    }
    else if(geom == "jitter") {
      gp = gp + geom_jitter()
    }
    else if(geom == "line") {
      #gp = gp + geom_point() + geom_path(aes_string(x="feature",y=value.name,
      #                                               group=".record.id"))
      ## We need aes(group=) here in stat_summary, but somehow not in other places
      gp = gp + stat_summary(aes_string(group=color),fun.y=stat_summary.fun.y, geom="line")
      
      #show.samp.n = F
    }
    else if(geom == "line_obs") {
      gp = gp + geom_path(aes_string(x="feature",y=value.name,
                                     group=".record.id"))
      if(line.show.points) {
        gp = gp + geom_point()
      }
      show.samp.n = F
    }
    else {
      stop(paste("Unexpected parameter value: geom = ",geom))
    }
    if(geom %in% c("line","line_obs")) {
      if(!is.na(line.legend.thickness)) {
        gp = gp + guides(colour = guide_legend(override.aes = list(size=line.legend.thickness)))
      }
    }
    else {
      if(!is.null(id.var.dodge)) {
        if(flip.coords) {
        ## create alternated shading and grid lines between the categories in order to visually separate
        ## groups of dodged geoms
        shading = data.frame(min = seq(from = 0.5, to = max(as.numeric(as.factor(dat$feature))), by = 2),
                             max = seq(from = 1.5, to = max(as.numeric(as.factor(dat$feature))) + 0.5, by = 2))
        gp = gp + geom_rect(data = shading,
                            aes(xmin = min, xmax = max, ymin = -Inf, ymax = Inf),
                            color="grey",alpha = 0.02,
                            inherit.aes=F)
        }
        if(flip.coords) hide.panel.grid.major.y = T
        else hide.panel.grid.major.x = T
      }
    }
    
    #stat_summary(fun.data = mean_cl_boot, geom = "pointrange",color="black")+
    #coord_flip()+
    #geom_boxplot(color="black")+
    #geom_point(position = "jitter")+
    #stat_identity(geom="bar")+
    #geom_bar(stat="identity")
    #scale_fill_brewer(type = "seq", palette = 1)
    #labels facet with number of cases
    if(flip.coords) {
      gp = gp + coord_flip()
    }

    ## flip the legend order so that both the legend and
    ## dodged geomes go from bottom to top if flipped coords,
    ## and the legend order would be the same in non-flipped
    gp = gp + guides(fill = guide_legend(reverse=TRUE),
                     color = guide_legend(reverse=TRUE))
     
    if(sqrt.scale) {
      gp = gp + coord_trans(y = "signed_sqrt")
    }
    if(length(id.vars.facet) == 0) {
      wr = facet_null()
    }
    else if (length(id.vars.facet) == 1) {
      wr = facet_wrap(facet.form,ncol = facet_wrap_ncol)
    }
    else {
      wr = facet_grid(facet.form,
                      drop=T,margins=facet_grid.margins)
    }
    gp = gp + wr
    
    if(!is.null(id.var.dodge)) {
      legend.position = "right"
    }
    else {
      legend.position = "none"
    }
  }  
  
  theme_font_size_abs = ggplot2::theme_get()$text$size
  #theme_font_size = 0.5 #if output is png
  n_feat_mult = (20/length(features))/(n.facet[[1]]/(if(length(id.vars.facet) == 1) facet_wrap_ncol else 1.))
  fontsize = theme_font_size*(n_feat_mult**0.33)
  
  if(length(id.vars.facet) > 0) {
    feature.names = as.character(features)
    #this will be used to label each facet with number of cases in it
    facet.cnt <- plyr::ddply(.data=data, id.vars.facet, function(x,feature.names) 
    { c(.n=nrow(x),
        .y=mean(apply(as.matrix(x[,feature.names]),2, function(z) max(z,na.rm=T)),
                names=F,na.rm=T)*(if(!flip.coords) 2./3 else 1./3),
        .x=max(length(feature.names)/3,1)) },
    feature.names)
    if(geom=="bar_stacked") {
      facet.cnt$.x = max(facet.cnt$.n/3,1)
      max.x.val = max(facet.cnt$.n)
    }
    facet.cnt$.n = paste("n =", facet.cnt$.n)
    #facet.cnt$y = facet.cnt$V2
    
    if(identical(stat_summary.fun.y,"identity")) {
      show.samp.n = F
    }
    if(show.samp.n) {
      # there is some bug in how ggplot2 treats size=theme_font_size_abs here if it is
      # not in the aes - it either gets tiny with rel() or huge with the absoulte
      # value. In the aes, size comes out as expected, but then it gets into
      # the legend.
      gp = gp +
        geom_text(aes(x=.x, y=.y, label=.n), 
                  data=facet.cnt, 
                  colour="black", 
                  inherit.aes=F, 
                  size=theme_font_size_abs/3,
                  parse=FALSE)
    }
  }
  else {
    if(geom=="bar_stacked") {
      max.x.val = nrow(data)
    }
  }
  
  if(geom=="bar_stacked") {
    if(max.x.val > 40) {
      hide.ticks.x = T
    }
  }
  
  #gp = gp + guides(colour = NULL, fill = guide_legend("XXX"))
  if(!is.null(legend.title)) {
    gp = gp + labs(fill=legend.title,colour=legend.title)
  }
  
  color.palette="brew"
  
  if(color.palette=="brew") {
    #n.color.orig = 8
    #palette = brewer.pal(n.color.orig, "Accent")
    #get.palette = colorRampPalette(palette)
    #palette = get.palette(max(length(features),n.color.orig))
    #palette = rep_len(palette,max(length(features),1000))
    palette = generate.colors.mgsat(palette.levels,value="palette")    
    gp = gp + 
      scale_fill_manual(values = palette) +
      scale_color_manual(values = palette)
  }
  else if(color.palette=="hue") {
    gp = gp +     
      scale_fill_hue(c = 50, l = 70, h=c(0, 360)) +
      scale_color_hue(c = 50, l = 70, h=c(0, 360))
  }
  #+ theme_grey(base_size = fontsize*theme_font_size_abs) + 
  theme_args = list(text=element_text(color=c("black","black"),size = fontsize*theme_font_size_abs),
                    legend.position = legend.position,
                    axis.title=element_blank(),
                    axis.text.y=element_text(color=c("black","black"),size = rel(if(sqrt.scale) axis.text.rel.y else axis.text.rel.y)),
                    plot.title = element_text(size = rel(1)),
                    axis.text.x = element_text(size = rel(if(!flip.coords || geom == "bar_stacked") axis.text.rel.x else axis.text.rel.x*1.25),
                                               angle=if(flip.coords) 0 else 90, hjust = 1))
  if(hide.ticks.x) {
    hide.axis.text.x = TRUE
    gp = gp + 
      scale_x_discrete(breaks=NULL) +
      scale_y_continuous(expand = c(0,0))
    theme_args[["axis.text.x"]]=element_blank()
    theme_args[["axis.ticks.x"]]=element_blank()
  }
  if(hide.axis.text.x) {
    theme_args[["axis.text.x"]]=element_blank()
  }
  if(hide.panel.grid.major.y) theme_args[["panel.grid.major.y"]]=element_blank()
  if(hide.panel.grid.major.x) theme_args[["panel.grid.major.x"]]=element_blank()
  
  gp = gp + do.call(theme,theme_args)
  if (!is.null(ggp.comp)) {
    for (g.c in ggp.comp) {
      gp = gp + g.c
    }
  }
  if (!is.null(file_name)) {
    ggsave(file_name)
  }
  return (new_mgsatres(plot=gp,dat.summary=dat.summary,dat.melt=dat.melt))
}

read.table.m_a <- function(file.base) {
  fn.count = paste(file.base,"count.tsv",sep=".")
  count = as.matrix(read.table(fn.count,
                               sep="\t",
                               header=T,
                               stringsAsFactors=T))
  fn.attr = paste(file.base,"attr.tsv",sep=".")
  attr = read.table(fn.attr,
                    sep="\t",
                    header=T,
                    stringsAsFactors=T)
  rownames(count) = attr$SampleID
  return (list(count=count,attr=attr))
}

write.table.m_a <- function(m_a,file.base,row.names=F) {
  fn.count = paste(file.base,"count.tsv",sep=".")
  write.table(m_a$count,
              fn.count,
              sep="\t",
              row.names = row.names)
  fn.attr = paste(file.base,"attr.tsv",sep=".")
  write.table(m_a$attr,
              fn.attr,
              sep="\t",
              row.names = row.names)
  return (list(fn.count=fn.count,fn.attr=fn.attr))
}

write.table.file.report.m_a = function(m_a,name.base,descr=NULL,row.names=F) {
  ## if we write row.names, Excel shifts header row to the left when loading
  file.base = report$make.file.name(name.base)
  files = write.table.m_a(m_a=m_a,
                          file.base=file.base,
                          row.names = row.names)
  if (!is.null(descr)) {
    links = sapply(files,
                   pandoc.link.verbatim.return
    )
    report$add.descr(paste("Wrote counts and metadata for",
                           descr,
                           "to files",
                           paste(links,collapse=",")))
  }
  return (files) 
}


export.taxa.meta <- function(m_a,
                             label,
                             descr="",
                             row.proportions=T,
                             row.names=F) {
  write.table.file.report.m_a(m_a=m_a,
                              name.base=paste("samples.raw",label,sep="."),
                              descr=paste("raw counts",descr),
                              row.names=row.names)  
  
  if(row.proportions) {
    
    write.table.file.report.m_a(m_a=norm.prop.m_a(m_a),
                                name.base=paste("samples.proportions",label,sep="."),
                                descr=paste("proportions counts",descr),
                                row.names=row.names)
  }
}

## Generate index into x that selects equal number
## of elements for each level of x. If target==0,
## target will be set to the smallest level count.
## If target is a character vector of levels,
## target will be set to the smallest count among
## only those levels
## If drop.smaller, indexes for levels below target
## will be removed. Otherwise, they will be returned
## with the original number of elements, but in random
## order.
## If there is a need to also balance within a specific
## strata (e.g. by phenotype but preserving as many
## family connections as possible), then we will have
## to use methods from package `sampling`.

balanced.sample <- function(x, target=NULL, drop.smaller=T) {
  library(data.table)
  x = data.table(lev=factor(x))[,ind:=.I]
  if(is.null(target)) target = levels(x$lev)
  if(is.character(target)) {
    n.target = min(x[lev %in% target,.N,by=lev]$N)
  }
  else {
    n.target = target
  }
  x = x[,.(ind=ind[sample(.N,min(.N,n.target))],
           n.target=n.target,
           n.lev=.N),by = lev]
  if(drop.smaller) {
    x = x[n.lev>=n.target]
  }
  x$ind
}

mgsat.lof <- function (data, k, cores = NULL, ...) 
{
  library(Rlof) #we use internal function here
  library(foreach)
  library(doParallel)
  
  if (is.null(k)) 
    stop("k is missing")
  if (!is.numeric(k)) 
    stop("k is not numeric")
  if (!is.numeric(cores) && !is.null(cores)) 
    stop("cores is not numeric")
  is.dist = inherits(data,"dist")
  data <- as.matrix(data)
  if (!is.numeric(data)) 
    stop("the data contains non-numeric data type")
  v.k <- as.integer(k)
  if (max(v.k) >= dim(data)[1]) 
    stop("the maximum k value has to be less than the length of the data")
  distdata <- mgsat.f.dist.to.knn(data, max(v.k), cores, is.dist, ...)
  p <- dim(distdata)[2L]
  dist.start <- as.integer((dim(distdata)[1])/2)
  dist.end <- dim(distdata)[1]
  ik <- numeric()
  registerDoParallel(cores = cores)
  m.lof <- foreach(ik = v.k, .combine = cbind) %dopar% {
    lrddata <- Rlof:::f.reachability(distdata, ik)
    v.lof <- rep(0, p)
    for (i in 1:p) {
      nneigh <- sum(!is.na(distdata[c((dist.start + 1):dist.end), 
                                    i]) & (distdata[c((dist.start + 1):dist.end), 
                                                    i] <= distdata[(dist.start + ik), i]))
      v.lof[i] <- sum(lrddata[distdata[(1:nneigh), i]]/lrddata[i])/nneigh
    }
    v.lof
  }
  if (length(v.k) > 1) 
    colnames(m.lof) <- v.k
  return(m.lof)
}

mgsat.f.dist.to.knn <- function (dataset, neighbors, cores, is.dist=F, ...) 
{
  library(foreach)
  library(doParallel)
  if(is.dist) {
    m.dist <- dataset
  }
  else {
    m.dist <- as.matrix(distmc(dataset, ...))
  }
  num.col <- dim(m.dist)[2]
  l.knndist <- lapply(c(1:num.col), function(i) {
    order.x <- order(m.dist[, i])
    kdist <- m.dist[, i][order.x[neighbors + 1]]
    numnei <- sum(m.dist[, i] <= kdist)
    data.frame(v.order = order.x[2:numnei], v.dist = m.dist[, 
                                                            i][order.x[2:numnei]])
  })
  rm(m.dist)
  maxnum <- max(unlist(lapply(l.knndist, function(x) {
    dim(x)[1]
  })))
  registerDoParallel(cores = cores)
  i <- numeric()
  knndist <- foreach(i = 1:num.col, .combine = cbind) %dopar% 
  {
    len <- dim(l.knndist[[i]])[1]
    c(l.knndist[[i]]$v.order, rep(NA, (maxnum - len)), 
      l.knndist[[i]]$v.dist, rep(NA, (maxnum - len)))
  }
  knndist
}

mgsat.find.outliers <- function(x,k=0.5,pval.adjust="BH",lof.dist.args=list(),alpha=alpha) {
  library(Rlof)
  library(fitdistrplus)
  if(inherits(x,"dist")) {
    n.obs = attr(x,"Size")
    rn = labels(x)
  } else {
    n.obs = nrow(x)
    rn = rownames(x)
  }
  stopifnot(length(k)==1)
  if(k<1) {
    k = round(n.obs*k)
  }
  lof.val = do.call(mgsat.lof,c(
    list(x,k,cores=NULL),
    lof.dist.args
  ))
  dist.par = fitdist(lof.val,"gamma")
  p.val = p.adjust(pgamma(lof.val,dist.par$estimate["shape"],dist.par$estimate["rate"],lower.tail = F),pval.adjust)
  if(!is.null(rn)) {
    names(p.val) = rn
  }
  return (list(p.val.adj=p.val,lof.val=lof.val,dist.par=dist.par,alpha=alpha,
               ind.outlier=which(p.val<=alpha)))
}

##TODO: comparing MDS and NMDS plots with outliers selected from the original
##distance matrix for sequence dissimilarity of HBV HBe protein, it appears
##that NMDS matches the original distances much better than MDS - the detected outliers
##are on the periphery of the plot, while for MDS they are often inside the large clusters.
##On the other hand, MDS plots seem to visually separate HBV genotypes better.
##We need to test on synthetic datasets that MDS ordination as implemented does not
##mismatch IDs and order of the observations.
mgsat.find.outliers.ordinate <- function(x,
                                         k=0.5,
                                         lof.on.ordinate=T,
                                         pval.adjust="BH",
                                         lof.dist.args=list(),
                                         ordinate.args=list(method="NMDS",k=3),do.plot=T,
                                         alpha=0.01,
                                         do.report=T) {
  library(phyloseq)
  if(inherits(x,"dist")) {
    m_a = list(attr=data.frame(name=labels(x),row.names = labels(x)))
    ordinate.args$distance = x
  }
  else {
    m_a = list(count=as.matrix(x,rownames.force=T))
  }
  ph = m_a.to.phyloseq(m_a)
  if(is.null(ordinate.args$method)) {
    ordinate.args$method = "NMDS"
    ordinate.args$k = 3
  }
  ord = do.call(phyloseq.ordinate,c(
    list(ph),
    ordinate.args
  ))
  
  ndim = ord$ndim
  if(is.null(ndim)) {
    ndim = 3
  }
  
  if(lof.on.ordinate) {
    ord.df = plot_ordination(ph, ord, type="samples",axes=1:ndim,justDF = T)
    out.data = ord.df[,1:ndim,drop=F]
  }
  else {
    out.data = x
  }
  out.res = mgsat.find.outliers(out.data,
                                k=k,
                                pval.adjust = pval.adjust,
                                lof.dist.args=lof.dist.args,
                                alpha=alpha)
  out.res$ord = ord
  if(do.plot) {
    if(ndim>=3) {
      out.res$pl3d = plot_ordination.3d(ph,ord,type="samples",
                                        color=log(out.res$p.val.adj),
                                        axes=1:ndim,
                                        labels=paste("ID:",names(out.res$p.val.adj),
                                                     "p.val.adj: ",format(out.res$p.val.adj,digits=3))
      )
      if(do.report) {
        report$add.widget(out.res$pl3d,
                          caption = "Ordination plot with observations colored according to p-value from
                   outlier detection procedure")
      }
      out.res$pl3d.outlier = plot_ordination.3d(ph,ord,type="samples",
                                                color=(seq_along(out.res$p.val.adj) %in% out.res$ind.outlier),
                                                axes=1:ndim,
                                                labels=paste("ID:",names(out.res$p.val.adj),
                                                             "p.val.adj: ",format(out.res$p.val.adj,digits=3))
      )
      if(do.report) {
        report$add.widget(out.res$pl3d.outlier,
                          caption = "Ordination plot with observations colored according to assigned
                   outlier status")
      }
    }
  }
  if(do.report) {
    out.df = data.frame(id=names(out.res$p.val.adj),
                        p.val.adj=out.res$p.val.adj,
                        is.outlier=(seq_along(out.res$p.val.adj) %in% out.res$ind.outlier))
    if(is.null(out.df$id)) {
      out.df$id = seq_along(out.df$p.val.adj)
    }
    out.df = out.df[order(out.df$p.val.adj),,drop=F]
    report$add.table(out.df,
                     caption = "Results of outlier detection")
  }
  return (out.res)
}

## If hill=T (default) return Hill numbers, otherwise - 
## Renyi entropies. Column names are prepended with N for Hill and
## with H for Renyi (to make them work nicely in formulas etc where naked numbers
## as column names cause problems). For evenness, E is also appended (e.g. NE_0.25).
mgsat.hill <- function(count,
                       scales = c(0,0.25, 0.5, 1, 2, 4, 8, Inf), 
                       evenness = F,
                       hill = T,
                       ...) {
  require(vegan)
  result <- renyi(count, scales = scales, ...)
  if (attributes(result)$class[2] == "numeric") {
    result <- as.matrix(result)
  }
  sep="_"
  if (evenness == T) {
    ## Jost 2010 http://dx.doi.org/10.3390/d2020207
    ## taking into an account that Hill = exp(Renyi)
    result = result[,colnames(result) != "0"] - renyi(count, scales = c(0))
    colnames(result) = paste("E",colnames(result),sep=sep)
    sep=""
  }
  if(hill == T) {
    result = exp(result)
    ##note that for q=Inf, this is a reciprocal of the frequency of most common
    ##species, and for q=-Inf - most rare species (PMC3470749)
    colnames(result) = paste("N",colnames(result),sep=sep)
  }
  else {
    colnames(result) = paste("H",colnames(result),sep=sep)    
  }
  return(result)
}

mgsat.diversity.hill.counts <- function(m_a,n.rar.rep=400,is.raw.count.data=T,do.rarefy=T,
                                        hill.args=list()) { 
  require(vegan)
  n.rar = min(rowSums(m_a$count))
  plus_here = function(x,y) (x+y)
  if(is.raw.count.data && do.rarefy) {
    x = foreach(seq(n.rar.rep),.packages=c("vegan"), .export=c("mgsat.hill"), .combine=plus_here,
                .final=function(x) (x/n.rar.rep)) %dopar% 
                {
                  do.call(
                    mgsat.hill,c(
                      list(rrarefy(m_a$count,n.rar)),
                      hill.args)
                  )
                }
  }
  else {
    x = do.call(mgsat.hill,c(
      list(m_a$count),
      hill.args)
    )
  }
  return(list(e=x))
}


mgsat.richness.counts <- function(m_a,
                                  n.rar.rep=400,
                                  do.rarefy=T,
                                  filtered.singletons=F) {
  
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  #S.ACE & all se.* give NaN often
  if(!filtered.singletons) {
    ind.names = c("S.obs","S.chao1")
  }
  else {
    ind.names = c("S.obs")
  }
  if(do.rarefy) {
    x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
                .final=function(x) (x/n.rar.rep)) %dopar% 
                {estimateR(rrarefy(m_a$count,n.rar))[ind.names,,drop=F]}
  }
  else {
    x = estimateR(m_a$count)[ind.names,,drop=F]
  }
  return(list(e=t(x)))
}

## This uses incidence data and therefore should be applicable to both raw count data as 
## well as to proportions
## or other type of measurement where non-zero value means "present"
mgsat.richness.samples <- function(m_a,group.attr=NULL,n.rar.rep=400,do.rarefy=T) {
  
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  
  if(is.null(group.attr)) {
    pool = factor(rep("All",nrow(m_a$count)))
    do.stratify = F
  }
  else {
    pool = m_a$attr[,group.attr]
    do.stratify = T
  }
  count = m_a$count
  if(!(do.stratify || do.rarefy)) {
    n.rar.rep = 1
  }
  #make.global(name="dbg1")
  ##somehow just supplying .combine="+" generates an error,
  ##but both the function below or skipping .combine and
  ##applying Reduce("+",...) on the returned list work fine
  plus <-function(x,y) (x+y)
  x = foreach(seq(n.rar.rep),.packages=c("vegan"),
              .combine=plus,
              .export=c("balanced.sample"),
              #DEBUG
              .final=function(x) (x/n.rar.rep)) %do% #%dopar% 
              {
                if(do.stratify) {
                  strat.ind = balanced.sample(pool)
                  count = count[strat.ind,,drop=F]
                  pool=pool[strat.ind]
                }
                if(do.rarefy) {
                  count.run = rrarefy(count,n.rar)
                }
                else {
                  count.run = count
                }
                specpool(count.run,pool=pool)
              }
  
  se.ind = grep(".*[.]se",names(x))
  e = x[,-se.ind,drop=F]
  se = x[,se.ind,drop=F]
  names(se) = sub("[.]se","",names(se))
  return(list(e=e,se=se,e.se=x))
}

## This returns Hill numbers
mgsat.diversity.alpha.counts <- function(m_a,n.rar.rep=400,is.raw.count.data=T,do.rarefy=T) {
  
  require(vegan)
  n.rar = min(rowSums(m_a$count))
  
  f.div = function(m) { cbind(N1=exp(diversity(m,index="shan")),
                              N2 = diversity(m,index="invsimpson")) }
  
  if(is.raw.count.data && do.rarefy) {
    x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
                .final=function(x) (x/n.rar.rep)) %dopar% 
                {f.div(rrarefy(m_a$count,n.rar))}
  }
  else {
    x = f.div(m_a$count)
  }
  
  if(is.raw.count.data) {
    #this is already unbiased (determenistic wrt rarefication)
    div.unb.simpson = rarefy(m_a$count,2)-1
    #div.fisher.alpha = fisher.alpha(m_a$count)
    x = cbind(x,div.unb.simpson=div.unb.simpson)
  }
  
  return(list(e=x))
}

mgsat.diversity.beta.dist.complexity <- function(n.r,n.c) {
  n.c*(n.r**2)
}

mgsat.diversity.beta.dist <- function(m_a,n.rar.rep=400,method="-1",do.rarefy=T) {
  
  require(vegan)
  n.rar = min(rowSums(m_a$count))
  
  if(do.rarefy) {
    
    one.iter.compl = mgsat.diversity.beta.dist.complexity(nrow(m_a$count),ncol(m_a$count))
    one.iter.compl.medium = mgsat.diversity.beta.dist.complexity(200,200)
    
    n.rar.rep = round(min(n.rar.rep,max(1,n.rar.rep * (one.iter.compl.medium / one.iter.compl))))
    
    x = foreach(seq(n.rar.rep),.packages=c("vegan"),.combine="+",
                .final=function(x) (x/n.rar.rep)) %dopar% 
                {
                  betadiver(rrarefy(m_a$count,n.rar),method=method)
                }
  }
  else {
    x = betadiver(m_a$count,method=method)
  }
  return(list(e=x))
}

mgsat.diversity.beta <- function(m_a,n.rar.rep=400,method="-1",
                                 group.attr=NULL,
                                 betadisper.task=list(),
                                 adonis.task=NULL,
                                 do.rarefy=T) {
  
  require(vegan)
  
  res = new_mgsatres()
  
  beta.dist = mgsat.diversity.beta.dist(m_a,n.rar.rep=n.rar.rep,method=method,do.rarefy=do.rarefy)$e
  
  method.help = paste(grep(sprintf('\"%s\"',method),capture.output(betadiver(help=T)),value=T),
                      ", where number of shared species in two sites is a, 
                      and the numbers of species unique to each site are b and c.",sep="")
  
  report$add.descr(sprintf("Computed beta-diversity matrix using function betadiver {vegan}
                   with method %s",
                           method.help))
  
  if(!is.null(group.attr)) {
    betadisp = do.call(betadisper,
                       c(
                         list(beta.dist,group=m_a$attr[,group.attr]),
                         betadisper.task
                       )
    )
    report$add.descr(sprintf("Results of function betadisper {vegan}
                       for the analysis of multivariate homogeneity of group dispersions.
                       This is applied to sample beta diversity matrix to analyze it with
                       respect to a grouping variable %s. Arguments for the call are: %s",
                             group.attr,
                             arg.list.as.str(betadisper.task)))
    anova.betadisp = anova(betadisp)
    report$add(anova.betadisp)
    res$anova.betadisp = anova.betadisp
    report$add(plot(betadisp),caption=sprintf("Results of betadisper {vegan}. Distances from samples 
               to the group
               centroids are shown in the first two principal coordinates.
               Groups are defined by the variable %s.
               Sample beta-diversity matrix was generated with method %s",
                                              group.attr,method.help))
  }
  
  if(!is.null(adonis.task)) {
    adonis.task$data.descr = sprintf("Beta-diversity dissimilarity matrix created with method %s",
                                     method.help)
    tryCatchAndWarn({
      m_a.bd = m_a
      m_a.bd$count = beta.dist
      res$adonis = do.call(test.counts.adonis.report,
                           c(
                             list(m_a=m_a.bd),
                             adonis.task
                           )
      )
    })
  }
  return(res)
}

## Modified vegan::rarefy code.
## Parallelized with nested foreach loops (requires active parallel backend).
## Sets to NA all output matrix elements where requested rarefication sample count
## is less than total available.
## include.max - if T, add maximum sample counts to the grid of sample sizes (adopted from 
## www.joshuajacobs.org/R/rarefaction.html)
## Value: named list(e=expected species count,se=standard error of e (if se=T in arguments))
rarefy.mgsat <- function (x, sample, se = FALSE, MARGIN = 1, include.max = T) 
{
  library(foreach)
  x <- as.matrix(x)
  if(MARGIN==2) {
    x = t(x)
  }
  sample.sums = rowSums(x)
  minsample <- min(sample.sums)
  if (any(sample > minsample)) 
    warning(gettextf("Requested 'sample' was larger than smallest site maximum (%d)", 
                     minsample))
  if (ncol(x) == 1 && MARGIN == 1) 
    x <- t(x)
  if (!identical(all.equal(x, round(x)), TRUE)) 
    stop("function accepts only integers (counts)")
  if(missing(sample)) {
    if(include.max) {
      sample = c()
    }
  }
  if(include.max) {
    sample = unique(sort(c(sample.sums,sample)))
  }
  if (missing(sample)) {
    sample <- minsample
    info <- paste("The size of 'sample' must be given --\nHint: Smallest site maximum", 
                  sample)
    stop(info)
  }
  rarefun <- function(x, sample) {
    x <- x[x > 0]
    J <- sum(x)
    if(J<sample) {
      out = NA
    }
    else {
      ldiv <- lchoose(J, sample)
      p1 <- ifelse(J - x < sample, 0, exp(lchoose(J - x, sample) - 
                                            ldiv))
      out <- sum(1 - p1)
    }
    if (se) {
      if(is.na(out)) {
        out = cbind(NA,NA)
      }
      else {
        V <- sum(p1 * (1 - p1))
        Jxx <- J - outer(x, x, "+")
        ind <- lower.tri(Jxx)
        Jxx <- Jxx[ind]
        V <- V + 2 * sum(ifelse(Jxx < sample, 0, exp(lchoose(Jxx, 
                                                             sample) - ldiv)) - outer(p1, p1)[ind])
        out <- cbind(out, sqrt(V))
      }
    }
    out
  }
  S.rare <- foreach(n=sample,.combine=cbind) %:% 
    foreach(i=1:nrow(x),.combine=rbind) %dopar% {
      rarefun(x[i,],sample=n)
    }
  S.rare <- matrix(S.rare, ncol = length(sample))
  colnames(S.rare) <- paste("", sample, sep = "")
  if (se) {
    S.rare.se = S.rare[seq(2,nrow(x),2),]
    S.rare.e = S.rare[seq(1,nrow(x),2),]
    rownames(S.rare.se) <- rownames(x)
  }
  else {
    S.rare.e = S.rare
  }
  rownames(S.rare.e) <- rownames(x)
  attr(S.rare.e, "Subsample") <- sample
  return (list(e=S.rare.e,se=if(se) S.rare.se else NULL))
}

mgsat.divrich.accum.plots <- function(m_a,is.raw.count.data=T,do.plot.profiles=T,plot.profiles.task=list(),do.rarefy=T) {
  
  res = new_mgsatres()
  
  if(is.raw.count.data) {
    max.n = max(rowSums(m_a$count))
    res.rar <- rarefy.mgsat(m_a$count, sample = seq(1,max.n,length.out = 30), se = FALSE, include.max = T) 
    if(do.plot.profiles) {
      plot.profiles.task$show.profile.task = within(plot.profiles.task$show.profile.task,{
        geoms = c("line_obs")
        dodged=T
        faceted=F
        line.show.points=F
        sqrt.scale=F
      })
      plot.profiles.task$feature.names.conv = function(x) as.numeric(as.character(x))
      do.call(plot.profiles,
              c(list(m_a=list(count=res.rar$e,attr=m_a$attr),
                     feature.descr="Abundance-based rarefaction curves",
                     value.name="Richness"),
                plot.profiles.task
              )
      )
    }
    res$rarecurve.abund = res.rar
  }
  return(res)
}


mgsat.divrich.accum.extra.plots <- function(m_a,is.raw.count.data=T,do.rarefy=T) {
  require(vegan)
  
  n.rar = min(rowSums(m_a$count))
  count = m_a$count
  rar.descr = ""
  if(do.rarefy) {
    count = rrarefy(m_a$count,n.rar)
    rar.descr = sprintf(" Samples were rarefied
             to the the minimum sample size (%s).",n.rar)
  }
  y = poolaccum(count)
  report$add(plot(y),caption=sprintf("Accumulation curves for extrapolated richness indices 
        for random ordering of samples (function poolaccum of package vegan;
             estimation is based on incidence data).%s",rar.descr))
  
  if(is.raw.count.data) {
    y = estaccumR(count)
    report$add(plot(y),caption=sprintf("Accumulation curves for extrapolated richness indices 
        for random ordering of samples (function estaccumR of package vegan;
             estimation is based on abundance data).%s",rar.descr))
  }
  
  y = specaccum(count,method="exact")
  report$add(plot(y, ci.type="polygon", ci.col="yellow",xlab="Size",ylab="Species"),
             caption=sprintf("Accumulation curve for expected number of species (features)
             for a given number of samples (function specaccum of package vegan,
             using method 'exact'.%s",rar.descr))
}

## Comparison of multiple correlated variables between two or more groups
## with Westfall & Young correction for multiple testing.
## Taken from simboot:::mcpHill. The original simboot function computes
## Hill numbers describing ecological diversity. In the present modified function, the variables (e.g.
## Hill numbers) must be supplied in a matrix `data` (observations are rows).
## The meaning of the remaining parameters and assumptions are the same as
## in the mcpHill - see its help page as well as the original Pallmann, P. et al. (2012)
## paper. In particular, the permutation based test assumes homoscedasticity.
mcp.wy <- function (data, fact, align = FALSE, block, boots = 5000, udmat = FALSE, 
                    usermat, mattype = "Dunnett", dunbase = 1, opt = "two.sided") 
{
  require(simboot) #or multicomp
  data = as.matrix(data)
  if (!is.matrix(data)) {
    stop("data must be an object of class 'matrix'.")
  }
  if (length(fact) != dim(data)[1]) {
    stop("The length of fact must equal the number of rows in dataf.")
  }
  if (length(levels(fact)) <= 1) {
    stop("The factor variable fact should have at least 2 levels to be compared.")
  }
  ni <- as.vector(summary(fact))
  if (udmat == FALSE) {
    if (mattype == "Dunnett") {
      cmat <- contrMat(ni, type = "Dunnett", base = dunbase)
    }
    else {
      cmat <- contrMat(ni, type = mattype)
    }
  }
  else {
    cmat <- usermat
  }
  qval = colnames(data)
  tabtab <- data
  group <- factor(fact)
  tabelle <- cbind(tabtab, group)
  if (align == TRUE) {
    alignfunc <- function(aliblock) {
      tfit <- lm(tabelle[, c(1:length(qval))] ~ aliblock)
      tabelleneu <- tabelle[, c(1:length(qval))] - predict(tfit)
    }
    tabelle2 <- alignfunc(block)
  }
  else {
    tabelle2 <- tabelle
  }
  funcA <- function(f) {
    fit <- lm(tabelle2[, c(1:length(qval))] ~ f - 1)
    epsilon <- residuals(fit)
  }
  epstabelle <- funcA(f = group)
  tstatshort <- function(mytab, i, f, cmat, ni) {
    mytab <- mytab[i]
    ni <- ni
    FIT <- lm(mytab ~ f - 1)
    mi <- coefficients(FIT)
    res <- residuals(FIT)
    varpool <- sum(res^2)/(sum(ni) - length(mi))
    estC <- (cmat %*% mi)
    varC <- (cmat^2) %*% (varpool/ni)
    ti <- estC/sqrt(varC)
    return(ti)
  }
  if (length(qval) == 1) {
    funcfunc <- function(abc, i, f, cmat, ni) {
      tstatshort(mytab = abc, f = f, i = i, cmat = cmat, 
                 ni = ni)
    }
  }
  else {
    funcfunc <- function(abc, i, f, cmat, ni) {
      apply(abc, 2, FUN = function(xx) {
        tstatshort(mytab = xx, f = f, i = i, cmat = cmat, 
                   ni = ni)
      })
    }
  }
  wyboot <- boot(epstabelle, funcfunc, R = boots, stype = "i", 
                 f = group, cmat = cmat, ni = ni)
  laenge <- 1:(length(qval) * dim(cmat)[1])
  tact <- as.vector(funcfunc(abc = tabelle2[, c(1:length(qval))], 
                             i = 1:sum(ni), f = group, cmat = cmat, ni = ni))
  bothfunc <- function(bobo) {
    apply(bobo, 1, function(huhu) max(abs(huhu)))
  }
  tboth <- bothfunc(wyboot$t)
  bothpval <- function(k) {
    tact <- tact[k]
    sapply(k, function(ka) sum(tboth > abs(tact[ka]))/boots)
  }
  maxfunc <- function(mama) {
    apply(mama, 1, max)
  }
  tmax <- maxfunc(wyboot$t)
  grpval <- function(k) {
    tact <- tact[k]
    sapply(k, function(ka) sum(tmax > tact[ka])/boots)
  }
  minfunc <- function(mimi) {
    apply(mimi, 1, min)
  }
  tmin <- minfunc(wyboot$t)
  lepval <- function(k) {
    tact <- tact[k]
    sapply(k, function(ka) sum(tmin < tact[ka])/boots)
  }
  allqs <- rep(qval, each = dim(cmat)[1])
  names(allqs) <- rep(rownames(cmat), times = length(qval))
  switch(opt, two.sided = cbind(q = allqs, `p-value` = bothpval(laenge)), 
         greater = cbind(q = allqs, `p-value` = grpval(laenge)), 
         less = cbind(q = allqs, `p-value` = lepval(laenge)))
}

mgsat.divrich.counts.glm.test <- function(m_a.divrich,
                                          divrich.names=NULL,                                          
                                          formula.rhs,
                                          glm.args=list(),
                                          do.plot.profiles=F,
                                          plot.profiles.task=NULL) {
  if(is.null(divrich.names)) {
    divrich.names = colnames(m_a.divrich$count)
  }
  
  res = new_mgsatres()
  
  for(divrich.name in divrich.names) {
    
    family = "gaussian"
    form.str = sprintf("%s~%s",divrich.name,formula.rhs)
    do.model = F
    
    if(divrich.name %in% c("S.obs","S.chao1","S.ACE")) {
      descr = paste("Richness estimate",divrich.name)
      do.model = T
    }
    else {
      descr = paste("Hill number of order",divrich.name)
      do.model = T      
    }
    
    if(do.model) {
      ##To see how well normal fits:
      ##library(fitdistrplus)
      ##descdist(x,boot=1000)
      ##gof = gofstat(fitdist(x,"norm"))
      ##gof$kstest #conservative, will properly reject only at high sample count
      mod = do.call(glm,
                    c(
                      list(
                        formula(form.str),
                        family=family,
                        data=cbind(as.data.frame(m_a.divrich$count),m_a.divrich$attr)
                      ),
                      glm.args
                    )
      )
      
      report$add(summary(mod),
                 caption=sprintf("Association of abundance based %s with sample metadata.
                                 GLM with family %s and formula %s",
                                 descr,family,form.str)
      )

      if(do.plot.profiles) {
        plot.profiles.task$show.profile.task = within(plot.profiles.task$show.profile.task,{
          geoms = c("boxplot","dotplot")
          dodged=T
          faceted=F
          theme_font_size=0.5
          hide.ticks.x=T
          flip.coords=F
        })
        do.call(plot.profiles,
                c(list(m_a=list(count=m_a.divrich$count[,divrich.name,drop=F],attr=m_a.divrich$attr),
                       feature.descr=sprintf("Abundance-based %s",
                                             descr),
                       value.name="index"),
                  plot.profiles.task
                )
        )
      }
            
      res[[divrich.name]] = mod
    }  
  }
  
  return(res)
  
}

mgsat.plot.richness.samples <- function(rich,var.names=c("chao","jack1","boot")) {
  
  var.names.se = var.names
  e = rich$e[,var.names]
  e$Group = rownames.as.factor(e)
  
  se = rich$se[,var.names.se]
  se$Group = rownames.as.factor(se)
  
  e.m = melt(e,"Group",var.names,variable.name="Index",value.name="e")
  se.m = melt(se,"Group",var.names.se,variable.name="Index",value.name="se")
  
  data = plyr::join(e.m,se.m,by=c("Group","Index"))
  dodge <- position_dodge(width=0.9)  
  pl = ggplot(data, aes(x = Index, y = e, fill = Group)) +  
    geom_bar(position = dodge,stat="identity",width=0.8) + 
    geom_errorbar(position = dodge,aes(ymin=e-se, ymax=e+se,width=0.5)) +
    xlab("Index name") +
    ylab("Index estimate and standard error")
  
  return(pl)
}

new_divrich <- function(...) {
  x = new_mgsatres(...)
  class(x) <- append(class(x),"divrich",0)
  return(x)
}

mgsat.divrich.report <- function(m_a,
                                 n.rar.rep=400,
                                 is.raw.count.data=T,
                                 filtered.singletons=F,
                                 group.attr=NULL,
                                 counts.glm.task=NULL,
                                 counts.genesel.task=NULL,
                                 beta.task=NULL,
                                 plot.profiles.task=list(),
                                 do.plot.profiles=T,
                                 do.plot.profiles.glm=T,
                                 do.incidence=T,
                                 do.abundance=T,
                                 do.abundance.richness=F,
                                 do.beta=T,
                                 do.rarefy=T,
                                 do.accum=T,
                                 extra.header="") {
  
  report.section = report$add.header(sprintf("Richness and diversity estimates %s",extra.header),
                                     section.action="push", sub=T)
  group.descr = " Samples are not grouped."
  group.descr.short = " for all samples"
  if(!is.null(group.attr)) {
    if(do.incidence) {
      group.descr = sprintf(" Incidence-based estimates are computed on sample pools split by
                          metadata attribute %s, and in each repetition, samples are also
                          stratified to balance the number of samples at each level
                          of the grouping variable.", group.attr)
    }
    else {
      group.descr = sprintf(" Samples are grouped by %s.", group.attr)
    }
    group.descr.short = sprintf(" for samples grouped by %s",group.attr)
  }
  if(!is.raw.count.data) {
    do.rarefy = F
  }
  
  rar.descr = " Counts are not rarefied."
  rar.descr.short = "No rarefication."
  if(is.raw.count.data && do.rarefy) {
    
    n.rar = min(rowSums(m_a$count))
    
    rar.descr = sprintf(" Counts are rarefied to the lowest library size (%s), abundance-based and
                   incidence-based alpha diversity indices and richness estimates are computed
                   (if requested).
                   This is repeated multiple times (n=%s), and the results are averaged.
                   Beta diversity matrix is also computed by averaging over multiple 
                   rarefications.",
                        n.rar, n.rar.rep)
    rar.descr.short = "With rarefication."
  }
  
  descr = sprintf("%s%s Dimensions of input count matrix (%s x %s).",rar.descr,group.descr,nrow(m_a$count),ncol(m_a$count))
  
  if(!str_blank(descr)) {
    report$add.descr(descr)
  }
  
  report$add.package.citation("vegan")
  
  res = new_divrich()
  
  if( do.incidence) {
    ##even if we filtered singletons, incidence based estimators are probably still fine - singletons
    ##for them means a species that was observed only in a single sample, and filtering out singletons
    ##merely changes our definition of "observed"
    res$rich.samples = mgsat.richness.samples(m_a,group.attr=group.attr,n.rar.rep=n.rar.rep,do.rarefy=do.rarefy)
    caption.inc.rich=sprintf("Incidence based rihcness estimates and corresponding standard errors%s",
                             group.descr.short)
    report$add.table(res$rich.samples$e.se,
                     caption=caption.inc.rich
    )
    report$add(mgsat.plot.richness.samples(res$rich.samples),
               caption=caption.inc.rich
    )
  }
  
  if(do.abundance) {
    if(do.abundance.richness && is.raw.count.data) {
      res$rich.counts = mgsat.richness.counts(m_a,n.rar.rep=n.rar.rep,
                                              do.rarefy=do.rarefy,
                                              filtered.singletons=filtered.singletons)
      if(do.plot.profiles) {
        do.call(plot.profiles,
                c(list(m_a=list(count=as.matrix(res$rich.counts$e),attr=m_a$attr),
                       feature.descr=sprintf("Abundance-based richness estimates %s",rar.descr.short),
                       value.name="Richness.Estimate"),
                  plot.profiles.task
                )
        )
      }
    }
    
    res$div.counts = list()
    for(div.task in list(list(id.task="diversity",evenness=F,descr="diversity indices (Hill numbers)"),
                         list(id.task="evenness",evenness=T,descr="evenness indices (Hill numbers / Observed 'species')")
    )
    ) {
      div.counts = mgsat.diversity.hill.counts(m_a,n.rar.rep=n.rar.rep,
                                               is.raw.count.data=is.raw.count.data,
                                               do.rarefy=do.rarefy,
                                               hill.args=list(evenness=div.task$evenness))
      if(do.plot.profiles) {
        plot.profiles.task$show.profile.task = within(plot.profiles.task$show.profile.task,{
          geoms = c("line","line_obs")
          dodged=T
          faceted=F
        })
        do.call(plot.profiles,
                c(list(m_a=list(count=div.counts$e,attr=m_a$attr),
                       feature.descr=sprintf("Abundance-based %s %s",
                                             div.task$descr, rar.descr.short),
                       value.name="index"),
                  plot.profiles.task
                )
        )
      }
      res$div.counts[[div.task$id.task]] = div.counts
      if(!is.null(group.attr)) {
        attr = factor(m_a$attr[,group.attr])
        lev.descr = data.frame(level=seq_along(levels(attr)),label=levels(attr))
        mcp.res = mcp.wy(div.counts$e,attr)
        report$add.package.citation("simboot")
        report$add.table(mcp.res,caption=sprintf("Comparison of %s
                         with Westfall and Young correction for multiple testing
                          across levels of attribute %s",div.task$descr,group.attr),
                         show.row.names=T)
        report$add.table(lev.descr,caption="Levels that defined contrasts 
                         in the previous table")
      }
    }
    divrich.counts = cbind(res$div.counts[["diversity"]]$e,res$div.counts[["evenness"]]$e)
    if(do.abundance.richness && is.raw.count.data) {
      divrich.counts = cbind(divrich.counts,res$rich.counts$e)
    }
    
    m_a.dr=list(count=as.matrix(divrich.counts),attr=m_a$attr)
    
    write.table.file.report.m_a(m_a=m_a.dr,
                                name.base="divrich.counts",
                                descr="Abundance based richness and diversity")
    
    if(!(is.null(counts.glm.task) || is.null(counts.glm.task$formula.rhs))) {
      
      res$glm.res = do.call(mgsat.divrich.counts.glm.test,
                            c(list(m_a.dr),
                              counts.glm.task,
                              list(
                                do.plot.profiles=do.plot.profiles.glm,
                                plot.profiles.task=plot.profiles.task
                              )
                            )
      )
      
    }
    
    if(!is.null(counts.genesel.task)) {
      
      ## note that while group.attr above may have any number of
      ## levels, the counts.genesel.task$group.attr must be two-level
      tryCatchAndWarn({ 
        res$genesel.res = do.call(genesel.stability.report,
                                  c(list(m_a.dr),
                                    counts.genesel.task))
      })
      
    }
  }
  
  if(do.accum) {
    report$push.section(report.section)
    mgsat.divrich.accum.plots(m_a,
                              is.raw.count.data=is.raw.count.data,
                              do.rarefy=do.rarefy,
                              plot.profiles.task=plot.profiles.task,
                              do.plot.profiles=do.plot.profiles)
    report$pop.section()
  }
  
  if(do.beta && !is.null(beta.task)) {
    res$beta = do.call(mgsat.diversity.beta,
                       c(list(m_a,
                              n.rar.rep=n.rar.rep,
                              group.attr=group.attr),
                         beta.task
                       )
    )
  }
  
  report$pop.section()
  
  return(res)
}



plot.profiles <- function(m_a,
                          ci=NULL,
                          feature.order=NULL,
                          feature.names.conv = NULL,
                          id.vars.list=list(c()),
                          feature.meta.x.vars=c(),
                          do.profile=T,
                          do.feature.meta=T,
                          value.name="Abundance",
                          show.profile.task=list(
                            geoms=c("bar","violin","boxplot","bar_stacked","line","line_obs"),
                            dodged=T,
                            faceted=T,
                            stat_summary.fun.y="mean",
                            sqrt.scale=F,
                            line.show.points=T,
                            line.legend.thickness=rel(3),
                            legend.title=NULL,
                            record.label=NULL,
                            theme_font_size=0.8,
                            width=NULL,
                            height=NULL,
                            hi.res.width=NULL,
                            show.samp.n=T,
                            axis.text.rel=1.,
                            axis.text.rel.x=1.,
                            axis.text.rel.y=1.,
                            flip.coords=c(T,F),
                            n.top=20
                          ),
                          show.feature.meta.task=list(),
                          feature.descr="Abundance.") {
  
  if(!grepl("\\.$",feature.descr)) {
    feature.descr = paste(feature.descr,".",sep="")
  }
  report.section = report$add.header(sprintf("Plots of %s",feature.descr),
                                     section.action="push", sub=T)
  report$add.descr("Plots are shown with relation to various combinations of meta 
                   data variables and in different graphical representations. Lots of plots here.")
  
  report$add.header("Iterating over all combinations of grouping variables")
  report$push.section(report.section)
  
  if(is.null(feature.order)) {
    feature.order = list(list(ord=NULL,ord_descr="original"))
  }
  
  if(length(feature.meta.x.vars)==0) {
    do.feature.meta = F
  }
  
  for (id.vars in id.vars.list) {
    
    if(length(id.vars)>0) {
      msg = paste("Grouping variables",paste0(id.vars,collapse=","))
    }
    else {
      msg = "Entire pool of samples"
    }
    report$add.header(msg)
    
    
    report$add.header(sprintf("Iterating over %s profile sorting order",feature.descr))
    report$push.section(report.section)
    
    for(pl.par in feature.order) {
      report$add.header(sprintf("%s profile sorting order: %s",feature.descr,pl.par$ord_descr))
      if(do.feature.meta) {
        feature.names.meta=if(is.null(pl.par$ord)) colnames(m_a$count) else pl.par$ord
        feature.names.meta = feature.names.meta[1:min(length(feature.names.meta),10)]
        
        report$add.header("Iterating over meta data variables")
        report$push.section(report.section)
        
        for(x.var in feature.meta.x.vars) {
          
          group.var = id.vars[id.vars != x.var]
          if(length(group.var)>0) {
            group.var = group.var[1]
          }
          else {
            group.var = NULL
          }
          tryCatchAndWarn({
            do.call(show.feature.meta,
                    c(
                      list(m_a=m_a,
                           feature.names=feature.names.meta,
                           x.var=x.var,
                           group.var=group.var,
                           value.name=value.name,
                           vars.descr=feature.descr),
                      show.feature.meta.task
                    )
            )
          })
          
        }
        
        report$pop.section()
      }
      
      if(do.profile) {
        
        dat.summary.done = list()
        
        report$add.header("Iterating over dodged vs faceted bars")
        report$add.descr("The same data are shown in multiple combinations of graphical representations. 
                         This is the same data, but each plot highlights slightly different aspects of it.
                         It is not likely that you will need every plot - pick only what you need.")
        report$push.section(report.section)
        
        within(show.profile.task, {
          if(!(dodged || faceted)) {
            faceted = T
          }
        })
        id.vars.dodge = list()
        if(show.profile.task$faceted) {
          id.vars.dodge[["faceted"]] = list(dodge=NULL,descr="faceted")
        }
        if(show.profile.task$dodged && length(id.vars)>0) {
          id.vars.dodge[["dodged"]] = list(dodge=id.vars[1],descr="dodged")
        }
        
        for(id.var.dodge in id.vars.dodge) {
          
          report$add.header(paste(id.var.dodge$descr,"plots. Iterating over orientation and, optionally, scaling"))
          report$push.section(report.section)
          sqrt.scale = show.profile.task$sqrt.scale
          for(flip.c in show.profile.task$flip.coords) {
            other.params = list(
              flip.coords=flip.c,
              sqrt.scale=if(flip.c) F else sqrt.scale,
              descr=if(flip.c) "Plot is in flipped orientation, Y axis not scaled." else 
                paste("Plot is in original orientation", 
                                        if(sqrt.scale) ", Y axis SQRT scaled." else ".", sep=""))

            report$add.header(paste(feature.descr, 
                                    other.params$descr, "Iterating over plot geometry"))
            report$push.section(report.section)
            
            for(geom in show.profile.task$geoms) {
              ## "bar_stacked" is only compatible with some combinations of other
              ## parameters, skip otherwise; same for "dotplot"
              skip.bar_stacked = F
              if(other.params$flip.coords || 
                 !is.null(id.var.dodge$dodge) || 
                 #length(id.vars) > 1 ||
                 ncol(m_a$count) < 2) {
                skip.bar_stacked = T
              }
              skip.dotplot = F
              if(!is.null(id.var.dodge$dodge)) {
                skip.dotplot = T
              }
              
              if(!((geom == "bar_stacked" && skip.bar_stacked) ||
                   (geom == "dotplot" && skip.dotplot))) {
                
                tryCatchAndWarn({
                  id.vars.key = paste(id.vars,collapse="#")
                  if(!id.vars.key %in% dat.summary.done) {
                    make.summary.table = T
                    dat.summary.done[[length(dat.summary.done)+1]] = id.vars.key
                  }
                  else {
                    make.summary.table = F
                  }
                  pl.abu = plot.abund.meta(m_a=m_a,
                                           ci=ci,
                                           id.vars=id.vars,
                                           features.order=pl.par$ord,
                                           feature.names.conv = feature.names.conv,
                                           geom=geom,
                                           line.show.points = show.profile.task$line.show.points,
                                           line.legend.thickness = show.profile.task$line.legend.thickness,
                                           file_name=NULL,
                                           id.var.dodge=id.var.dodge$dodge,
                                           flip.coords=other.params$flip.coords,
                                           sqrt.scale=other.params$sqrt.scale,
                                           value.name=value.name,
                                           stat_summary.fun.y=show.profile.task$stat_summary.fun.y,
                                           make.summary.table = make.summary.table,
                                           legend.title = show.profile.task$legend.title,
                                           record.label = show.profile.task$record.label,
                                           theme_font_size = show.profile.task$theme_font_size,
                                           hide.ticks.x = show.profile.task$hide.ticks.x,
                                           hide.axis.text.x = show.profile.task$hide.axis.text.x,
                                           show.samp.n = show.profile.task$show.samp.n,
                                           axis.text.rel = show.profile.task$axis.text.rel,
                                           axis.text.rel.x = show.profile.task$axis.text.rel.x,
                                           axis.text.rel.y = show.profile.task$axis.text.rel.y,
                                           n.top = show.profile.task$n.top,
                                           dodge.width = show.profile.task$dodge.width
                  )
                  
                  pl.hist = pl.abu$plot
                  #env=as.environment(as.list(environment(), all.names=TRUE))
                  #print(names(as.list(env)))
                  #print(evals("pl.hist",env=env))
                  if(length(id.vars)>0) {
                    gr.by.msg = sprintf("Data grouped by %s.", paste(id.vars,collapse=","))
                  }
                  else {
                    gr.by.msg = "Data for all pooled samples."
                  }
                  geom.descr = geom
                  if(geom == "bar") {
                    geom.descr = sprintf("%s (sample %s)",
                                         geom.descr,
                                         show.profile.task$stat_summary.fun.y
                    )
                  }
                  
                  dat.melt = pl.abu$dat.melt
                  if(!is.null(dat.melt)) {
                    report$add.table(dat.melt,caption=paste("Data table used for plots.",gr.by.msg))
                  }
                  
                  dat.summary = pl.abu$dat.summary
                  if(!is.null(dat.summary)) {
                    report$add.table(dat.summary,caption=paste("Summary table.",gr.by.msg))
                  }
                  caption=paste(sprintf("%s %s",feature.descr,gr.by.msg),
                                if(!is.null(pl.par$ord)) 
                                {sprintf("Sorting order of features is %s.",pl.par$ord_descr)} 
                                else {""},
                                geom.descr,"plot.")
                  report$add(pl.hist,caption = caption,
                             width = first_non_null(show.profile.task$width,evalsOptions("width")),
                             height = first_non_null(show.profile.task$height,evalsOptions("height")),
                             hi.res.width = first_non_null(show.profile.task$height,evalsOptions("hi.res.width"))
                  )
                  
                })
              }              
            }
            report$pop.section()
          }
          report$pop.section()        
        }
        report$pop.section()
      }
    }
    report$pop.section()
  }
  report$pop.section()
  report$pop.section()
}


dirmult.from.df <- function(data,col_ignore=c()) {
  cnt_m = count_matr_from_df(data,col_ignore=col_ignore)
  #We do not use dirmult() from dirmult package directly
  #because it drops pi elements that are zero, which makes
  #further work with the result very clumsy.
  #However, keep in mind that various methods from HMP package
  #will return NaNs if you pass dirmult parameters that contain
  #zero elements. You will have to filter the parameter object
  #before passing it to such methods.
  dm.par = DM.MoM(cnt_m)
  return (dm.par)
}


dirmult.kelvin <- function(file_name,group_sel) {
  
  attr_names = c("id_repl","group")
  
  data_all = read.kelvin.summary.matr(file_name)
  data_all = count.filter(data_all,col_ignore=attr_names,min_max_frac=0.25,min_row_sum=5000)  
  dm.par = list()
  for (gr in group_sel) {
    data = data_all[data_all$group == gr,]
    dm.par[[gr]] = dirmult.from.df(data,col_ignore=attr_names)
  }
  return (dm.par)
}

get.feature.names <- function(data,attr.names) {
  colnames(data)[!(colnames(data) %in% attr.names)]
}

power.choc<-function(taxa.meta.data,taxa.meta.attr.names) {
  
  
  dm.par.kelv = dirmult.kelvin("v35.16sTaxa.TotFilt_1000.allSamples.summary_table.xls",c("Stool"))[[1]]
  
  
  abund_matr_norm = norm.prop(count_matr_from_df(taxa.meta.data,col_ignore=taxa.meta.attr.names))
  
  #assign("abund_matr_norm",abund_matr_norm,envir=globalenv())
  
  #We have only one paired sample. 
  #set pi as proportions from a patient before therapy. Keep theta form HMP dataset.
  dm.par.1 <- dm.par.kelv
  dm.par.1$pi = abund_matr_norm["LK3_91",]
  dm.par.1$gamma = dirmult.comp.gamma(dm.par.1$pi,dm.par.1$theta)
  
  #set pi from the same patient after therapy
  dm.par.2 <- dm.par.kelv
  dm.par.2$pi = abund_matr_norm["LK4_92",]
  dm.par.2$gamma = dirmult.comp.gamma(dm.par.2$pi,dm.par.2$theta)
  
  dm.par.orig = list(dm.par.1,dm.par.2)
  
  groups = as.factor(c("Before","After"))
  
  all_features = get.feature.names(taxa.meta.data,taxa.meta.attr.names)
  
  #test_features = as.factor(c("Actinobacteria_0.1.1.1","Clostridia_0.1.6.2"))
  test_features = as.factor(c("Lachnospiracea_incertae_sedis_0.1.6.2.1.3.7","Faecalibacterium_0.1.6.2.1.5.4","unclassified_0.1.6.2.1.5.9"))
  
  
  n.samp.grp = c(50,60)
  n.seq.range = seq(1500,2000,by=1)
  #n.seq.range <- rep(50,100)
  #n.samp.grp = c(8,5)
  
  #set number of sequences here so that we can estimate
  #effect size and DM test statistics outside of sample generation loop
  n.seq = list()
  
  for (i.group in seq(length(groups))) {
    group = groups[i.group]
    n.seq[[i.group]] <- sample(n.seq.range, size=n.samp.grp[i.group]) 
  }
  
  res.power = power.dirmult.range(
    dm.par.orig=dm.par.orig,
    n.seq=n.seq,
    groups=groups,
    all_features=all_features,
    test_features=test_features,
    n.samp.grp=n.samp.grp,
    effect.size.range = c(0, 0.1, 0.2, 0.3, 0.4, 1.0),
    n.rep = 100
  )    
  write.csv(res.power,"res.power.csv")
}

## Default method for loading metadata file
## This will work OK if your file follows certain conventions:
## 1. Tab delimited with a header row
## 2. All strings can be treated as factors
## 3. All fields that look like numbers are OK to treat like numbers.
##    Specifically, make sure that all ID fields (SampleID, SubjectID etc)
##    start with a letter, not a number.
## 4. There has to be a unique key field SampleID, that matches SampleID in abundance tables
##
## If your file cannot follow the conventions above, then write your own method and pass it to 
## read.data.project(). The method should create SampleID factor key field and also set row.names to
## the values of that field.

load.meta.default <- function(file.name) {
  meta = read.delim(file.name,header=T,stringsAsFactors=T)
  meta$SampleID = as.factor(meta$SampleID)
  row.names(meta) = meta$SampleID
  return (meta)
}

load.meta.american.gut <- function(file.name) {
  meta = read.delim(file.name,header=T,stringsAsFactors=T)
  ##header starts with #, which gets replaced with `X.`
  meta$SampleID = as.factor(meta$X.SampleID)
  row.names(meta) = meta$SampleID
  return (meta)
}


read.data.project.yap <- function(taxa.summary.file,
                                  otu.shared.file,
                                  cons.taxonomy.file,
                                  taxa.summary.file.otu,
                                  meta.file,
                                  load.meta.method,
                                  load.meta.options=list(),
                                  count.filter.options=NULL,
                                  count.basis=c("seq","otu"),
                                  sanitize=T,
                                  taxa.count.source=c("shared","summary"),
                                  otu.count.filter.options=NULL,
                                  taxa.level=3,
                                  taxa.levels.mix=0) {
  taxa.count.source = taxa.count.source[[1]]
  count.count.source.descr = sprintf(" with taxa count source %s",taxa.count.source)
  count.basis = count.basis[[1]]
  count.basis.descr = sprintf(" with count basis %s",count.basis)
  if (taxa.count.source == "shared" || taxa.level == "otu") {
    ## glob() to make it easy just providing directory name with extension-based globs
    otu.shared.file = Sys.glob(otu.shared.file)
    stopifnot(length(otu.shared.file)==1)
    cons.taxonomy.file = Sys.glob(cons.taxonomy.file)
    stopifnot(length(cons.taxonomy.file)==1)
    if(taxa.level=="otu") {
      ## Since we label taxa with OTU name anyway, just use the first deined rank name
      ## for the prefix
      taxa.levels.mix = Inf
    }
    taxa.lev.all = read.mothur.otu.with.taxa.m_a(otu.shared.file=otu.shared.file,
                                                 cons.taxonomy.file=cons.taxonomy.file,
                                                 sanitize=sanitize,
                                                 taxa.level=taxa.level,
                                                 count.basis=count.basis,
                                                 otu.count.filter.options=otu.count.filter.options,
                                                 taxa.levels.mix=taxa.levels.mix)
    report$add.p(sprintf("Loaded OTU taxonomy file %s.",
                         pandoc.link.verbatim.return(cons.taxonomy.file)
    ))
    count.file = otu.shared.file
  }
  else {
    count.file = switch(count.basis,
                        seq=taxa.summary.file,
                        otu=taxa.summary.file.otu
    )
    count.file = Sys.glob(count.file)
    stopifnot(length(count.file)==1)
    moth.taxa <- read.mothur.taxa.summary(count.file,sanitize=sanitize)
    taxa.lev.all = multi.mothur.to.abund.m_a(moth.taxa,taxa.level)
    
  }  
  report$add.p(sprintf("Loaded %i records for %i features from count file %s for taxonomic level %s 
                       with taxa name sanitize setting %s%s%s",
                       nrow(taxa.lev.all$count),ncol(taxa.lev.all$count),
                       pandoc.link.verbatim.return(count.file),
                       taxa.level,
                       sanitize,
                       count.basis.descr,
                       count.count.source.descr))
  taxa.lev = report.count.filter.m_a(taxa.lev.all,count.filter.options=count.filter.options)
  meta = do.call(load.meta.method,c(file.name=meta.file,load.meta.options))  
  m_a = merge.counts.with.meta(taxa.lev$count,meta)
  report$add.p(sprintf("After merging with metadata, %i records left",
                       nrow(m_a$count)))
  return (m_a)
}

## Note that physeq constructor will convert all fields in attr.taxa to strings
m_a.to.phyloseq <- function(m_a,attr.taxa=NULL) {
  require(phyloseq)
  if(is.null(m_a$count) && is.null(m_a$attr)) {
    stop("Need at least some data for the samples")
  }
  attr = m_a$attr
  if(inherits(attr,"data.table")) {
    attr = as.data.frame(attr)
    ##input data table must have row names in 'rn' field
    rownames(attr) = attr$rn
  }
  ##need at least two columns or some phyloseq methods lose the dimension (not using drop=F)
  if(is.null(m_a$count)) {
    count = matrix(0,nrow(attr),2)
    rownames(count) = rownames(attr)
    colnames(count) = paste("Dummy",seq(ncol(count)),sep=".")
  }
  else {
    count = m_a$count
  }
  otu = otu_table(count, taxa_are_rows = F)
  tax = data.frame(Feature=colnames(count),Dummy=colnames(count))
  rownames(tax) = tax[,"Feature"]
  if(!is.null(attr.taxa)) {
    attr.taxa = as.data.frame(attr.taxa)
    tax.m = merge(tax,attr.taxa,by="row.names")
    stopifnot(nrow(tax.m)==nrow(tax) && nrow(tax.m)==nrow(attr.taxa))
    rownames(tax.m) = tax.m$Row.names
    tax.m$Row.names=NULL
    tax = tax.m[rownames(tax),]
    tax.d = tax
  }
  tax = tax_table(as.matrix(tax))
  attr = sample_data(attr)
  phyloseq(otu,tax,attr)
}

read.american.gut.m_a <- function(biom.file,meta.file,taxa.level) {
  ph = import_biom(biom.file)
  meta = sample_data(load.meta.american.gut(meta.file))
  ph = merge_phyloseq(ph,meta)
  stopifnot(as.character(taxa.level) != "otu")
  taxa.level.ind = as.numeric(taxa.level)
  ph = tax_glom(ph,taxrank=rank_names(ph)[taxa.level.ind])
  count = otu_table(ph)
  taxa = substring(tax_table(ph)[,taxa.level.ind],4)
  taxa[str_blank(taxa)] = "Unclassified"
  rownames(count) = taxa
  count = t(count)
  attr = sample_data(ph)
  stopifnot(nrow(count)==nrow(attr))
  stopifnot(all(rownames(count)==rownames(attr)))
  list(count=count,attr=attr)
}

summary.meta.method.default <- function(taxa.meta) {
  
  report$add.header("Summary of metadata variables")
  
  m_a = split_count_df(taxa.meta$data,col_ignore=taxa.meta$attr.names)
  
  report$add.printed(summary(m_a$attr),caption="Summary of metadata variables")
  
}


mgsat.16s.task.template = within(list(), {
  
  label.base = "16s"
  
  main.meta.var = "Group"
  
  read.data.method=read.data.project.yap  
  
  read.data.task = list(
    taxa.summary.file=NULL,
    otu.shared.file=NULL,
    cons.taxonomy.file=NULL,
    count.basis="seq",
    sanitize=T,
    count.filter.options=list(),
    taxa.count.source=c("shared"),
    ##drop all low abundant OTUs before we aggregate them into taxonomic rank counts
    otu.count.filter.options=list(), #list(min_max_frac = 0.0001),
    meta.file=NULL,
    load.meta.method=load.meta.default,
    load.meta.options=list()
  )
  
  taxa.levels = c("otu",2,3,4,5,6)
  
  get.taxa.meta.aggr<-function(m_a) { return (m_a) }
  
  count.filter.sample.options=list(min_row_sum=2000,max_row_sum=400000)
  
  summary.meta.method=summary.meta.method.default
  
  do.summary.meta = F
  
  do.plots = T
  
  do.tests = T
  
  summary.meta.task = list(
    meta.x.vars = c(),
    group.vars = c(main.meta.var)
  )
  
  test.counts.task = within(list(), {
    
    do.deseq2 = T
    do.genesel=T
    do.stabsel=T
    do.glmer=T
    do.adonis=T
    do.divrich=c("otu",6)
    do.divrich.pre.filter=T
    do.divrich.post.filter=F    
    do.plot.profiles.abund=T
    do.heatmap.abund=T
    do.ordination=T
    do.network.features.combined=T
    do.select.samples=c()
    do.extra.method=c()
    do.aggr.after.norm=c()
    
    feature.ranking = "stabsel"
    
    alpha = 0.05
    
    divrich.task = within(list(),{
      n.rar.rep=400
      is.raw.count.data=T
      filtered.singletons=F
      do.abundance.richness=F
      group.attr = main.meta.var
      counts.glm.task = within(list(),{
        formula.rhs = main.meta.var
      })
      beta.task = within(list(),{
        method="-1"
        betadisper.task=list()
        adonis.task=NULL #will be replaced by global adonis.task
      })
      ## Same structure as the task-wide genesel.task; if NULL,
      ## this will be taken from task-wide structure
      counts.genesel.task = NULL
      do.plot.profiles = T
      do.plot.profiles.glm = T
      ## Computing beta-diversity matrix on multiple rarefications can take a while
      do.beta = T
      do.accum = T
      do.incidence = T
    })
    
    count.filter.feature.options=list(drop.unclassified=T,min_mean_frac=0.0005)
    
    norm.count.task = within(list(), {
      method = "norm.ihs.prop"
      method.args = list(theta=1)
      drop.features=list("other")
      ##method="norm.rlog.dds"
      ##method.args=list(dds=NA) #signals to pull Deseq2 object
    })
    
    deseq2.task = within(list(), {
      formula.rhs=main.meta.var
      ##parameters of DESeq call
      test.task=list(
        fitType="local"
      )
      ##list of parameter sets for DESeq2::result calls
      ##can leave empty, in which case one call to result()
      ##with default contrast will be made
      result.tasks=list(list())
    })
    
    stabsel.task = list(
      resp.attr = main.meta.var,      
      args.fitfun = list(
        family="binomial",
        standardize=T                                     
      ),
      args.stabsel = list(
        PFER=0.05,
        sampling.type="SS",
        assumption="r-concave",
        B=400
      )
    )
    
    
    genesel.task = within(list(), {
      group.attr = main.meta.var
      do.nmds = F
      do.plot.profiles = T
      norm.count.task = NULL
      genesel.param = within(list(), {
        block.attr = NULL
        type="unpaired"
        replicates=400
        samp.fold.ratio=0.5
        maxrank=20
        comp.log.fold.change=T
      })
    })
    
    select.samples.task = list (
      n.species = 20,
      n.samples = 20
    )
    
    adonis.task = within(list(), {
      
      tasks=list(list(
        formula.rhs=main.meta.var,
        strata=NULL,
        descr=paste("Association with",main.meta.var)
      ))
      
      n.perm=4000
      dist.metr="bray"
      col.trans="ident"
      norm.count.task = within(norm.count.task, {
        method = "norm.prop"
        method.args = list()
      })      
    })
    
    glmer.task = within(list(), {
      
      tasks = list(list(
        descr.extra = "",
        formula.rhs = paste(main.meta.var,"(1|SampleID)",sep="+"),
        linfct=c("YearsSinceDiagnosis = 0")
      ))
      
    })
    
    plot.profiles.task = within(list(), {
      id.vars.list = list(c(main.meta.var))
      feature.meta.x.vars=c()
      do.profile=T
      do.feature.meta=F
      show.profile.task=list(
        ##geoms can be also c("dotplot","jitter","line","line_obs")
        geoms=c("bar_stacked","bar","violin","boxplot"),
        dodged=T,
        faceted=T,
        stat_summary.fun.y="mean",
        sqrt.scale=T,
        line.show.points=F,
        line.legend.thickness=rel(3),
        facet_wrap_ncol=3,
        legend.title=NULL,
        ## record.label will control the order of bars in stacked bar plot,
        ## without that argument they will be ordered by abundance of the most
        ## dominant taxa
        record.label=NULL,
        hide.ticks.x=F,
        hide.axis.text.x=F,
        theme_font_size = 0.8,
        show.samp.n = T,
        width = NULL,
        height = NULL,
        hi.res.width = NULL,
        axis.text.rel = 1.,
        axis.text.rel.x = 1.,
        axis.text.rel.y = 1.,
        flip.coords = c(F,T),
        n.top = 20,
        dodge.width = 0.7
      )
      show.feature.meta.task=list()
    })
    
    plot.profiles.abund.task = within(list(), {
      
      norm.count.task = within(norm.count.task, {
        method="norm.prop"
        method.args=list()
      })
      
    })
    
    heatmap.abund.task = within(list(), {
      attr.annot.names=c(main.meta.var)
      attr.row.labels=NULL
      trans.clust=NULL
      stand.clust="range"
      dist.metr="bray"
      caption="Heatmap of abundance profile"
      trans.show="norm.boxcox"
      stand.show="range"
    })
    
    heatmap.combined.task = within(list(), {
      norm.count.task=NULL
      hmap.width=1000
      hmap.height=hmap.width*0.8
      attr.annot.names=c()
      clustering_distance_rows="pearson"
      km.abund=0
      km.diversity=0
      show_row_names=F
      max.n.columns=NULL
    })
    
    ordination.task = within(list(), {
      distance="bray"
      ord.tasks = list(
        list(
          ordinate.task=list(
            method="RDA"
            ##other arguments to phyloseq:::ordinate
          ),
          plot.task=list(
            type="samples",
            legend.point.size = ggplot2::rel(4),
            legend.position="right",
            ggplot.extra=list()
            ##line.args, axis.scale=c(1,1,1)
            ##other arguments to phyloseq:::plot_ordination
          )
        )
      )
    })
    
    network.features.combined.task = within(list(), { 
      count.filter.options=NULL                                                
      drop.unclassified=T
      method="network.spiec.easi"
      ## see network.spiec.easi.options for default options that can be overriden here
      method.options=list()
      descr=""
    }) 
    
    extra.method.task = within(list(), {
      func = function(m_a,m_a.norm,res.tests,...) {}
      ##possibly other arguments to func()
    })
    
    aggr.after.norm.task = within(list(), {
      func = function(m_a,m_a.norm,m_a.abs,res.tests,...) 
      {return(list(m_a=m_a,m_a.norm=m_a.norm,m_a.abs=m_a.abs))}
      ##possibly other arguments to func()
    })
    
  })
  
})


start.cluster.project <- function(parallel.type="PSOCK",bioc.backend=c("snow","parallel"),...) {
  library(foreach)
  library(doParallel)
  library(parallel)
  ## number of cores to use on multicore machines
  node.cores = getOption("mc.cores", 4)
  options(boot.ncpus = node.cores)
  ## parallel backend
  options(boot.parallel = "multicore")
  cl<-parallel::makeCluster(node.cores,type=parallel.type,...)
  registerDoParallel(cl)
  library("BiocParallel")
  ##should use MulticoreParam to use "parallel",
  ##but currently gives send/receive errors on MacOS R 3.4.1 Conda
  bioc.backend = bioc.backend[1]
  if(bioc.backend=="snow") {
    register(SnowParam(node.cores))
  }
  else if(bioc.backend=="parallel") {
    register(MulticoreParam(node.cores))
  }
  else {
    stop(sprintf("Unknown value for bioc.backend: %s",bioc.backend))
  }
  return(cl)
}

stop.cluster.project <- function(cl) {
  if(!is.null(cl)) {
    stopCluster(cl)
  }
}

## create new MGSAT result object
new_mgsatres <- function(...) {
  x = list(...)
  class(x) <- append(class(x),"mgsatres",0)
  return(x)
}

## extract ranking of features created by a specific method
## Value: named list with at least name `ranked` that is an ordered 
## vector of feature names.
get.feature.ranking <- function(x, ...) {
  UseMethod('get.feature.ranking', x)
}

get.feature.ranking.default <- function(x.f) {
  stop("Not defined for arbitrary objects")
}

get.feature.ranking.genesel <- function(x.f,only.names=T) {
  ranked = x.f$stab_feat
  if(only.names) {
    ranked = ranked$name
  }
  return(list(ranked=ranked))
}

get.feature.ranking.stabsel <- function(x.f,only.names=T) {
  ranked = x.f$max[order(x.f$max,decreasing = T)]
  if(only.names) {
    ranked = names(ranked)
  }
  return(list(ranked=ranked))
}

get.feature.ranking.deseq2 <- function(x.f,only.names=T,id.result=1) {
  if(is.numeric(id.result) && id.result == 0) {
    return (list(ranked=x.f))
  }
  else {
    ranked = as.data.frame(x.f$results[[id.result]])
    ranked = cbind(name=rownames(ranked),ranked)
    if(only.names) {
      ranked = ranked$name
    }
    return (list(ranked=ranked))
  }
}

get.feature.ranking.mgsatres <- function(x.f,method="stabsel",...) {
  meth.res = x.f[[method]]
  if(is.null(meth.res)) {
    res = NULL
  }
  else {
    if(method %in% c("stabsel","genesel","deseq2")) {
      res = get.feature.ranking(meth.res,...)
    }
    else {
      stop(paste("I do not know what to do for method",method))
    }
  }
  return (res)
}

## extract results of diversity analysis
get.diversity <- function(x, ...) {
  UseMethod('get.diversity', x)
}

get.diversity.default <- function(x.f) {
  stop("Not defined for arbitrary objects")
}

get.diversity.divrich <- function(x.f,type) {
  if(type %in% c("diversity","evenness")) {
    x.f$div.counts[[type]]
  }
  else {
    stop(paste("Unknown diversity type",type))    
  }
}

get.diversity.mgsatres <- function(x.f,type=NULL,...) {
  dr = x.f[["divrich"]]
  if(is.null(dr)) return (NULL)
  if(is.null(type)) {
    return (dr)
  }
  else {
    return (get.diversity(dr,type=type,...))
  }
}

report.count.filter.m_a <- function(m_a,count.filter.options=NULL,descr="",warn.richness=T) {
  if(is.null(count.filter.options)) {
    count.filter.options = list()
  }
  report$add.p(sprintf("Filtering abundance matrix with arguments %s. %s",
                       arg.list.as.str(count.filter.options),
                       descr))
  if(warn.richness) {
    report$add.p("Note that some community richness estimators will not work correctly 
               if provided with abundance-filtered counts")
  }
  m_a = do.call(count.filter.m_a,c(list(m_a),count.filter.options))
  report$add.p(sprintf("After filtering, left %i records for %i features",
                       nrow(m_a$count),ncol(m_a$count)))
  return (m_a)
}


proc.project <- function(
  task.generator.method
) {
  
  report$add.descr("Set of analysis routines is applied
                   in nested loops over combinations of defined subsets of samples (if any) and 
                   taxonomic levels. For each output, look for the nearest
                   headers to figure out its place in the report hierarchy.
                   The entire analysis is typically split into sub-reports linked
                   from higher-level pages. Follow links called 'Subreport'.
                   If viewing HTML formatted report, you can click on the
                   images to view the hi-resolution picture. To make it easier to
                   extract pictures for downstream use, picture files are also
                   reported as direct links in the legends.
                   Various intermediate datasets are also saved as delimited files and reported as 
                   direct links.
                   The HTML report is viewed best with modern versions of Chrome or Firefox browsers.
                   Internet Explorer might fail to show left-pane contents menu.")
  
  report.section = report$get.section()
  
  report$add.header("Iterating over subsets of data")
  report$push.section(report.section) #1 {
  
  tasks = task.generator.method()
  
  cl = start.cluster.project()
  
  res.tasks = lapply(tasks,function(task) {
    
    report$add.header(paste("Subset:",task$descr),
                      report.section=report.section,sub=T) #2 {
    
    res.task = new_mgsatres(task=task)
    
    res.task$res.taxa.levels = with(task,{
      
      report$add.header("Iterating over taxonomic levels")
      report$push.section(report.section) #3 {
      
      res.taxa.levels = list()
      
      for (taxa.level in taxa.levels) {
        
        res.level = new_mgsatres(taxa.level=taxa.level)
        
        label = paste(label.base,"l",paste(taxa.level,collapse="_"),sep=".",collapse=".")
        report$add.header(sprintf("Taxonomic level: %s of Subset: %s",taxa.level, descr),
                          report.section=report.section,sub=T) #4 {
        
        report$add.header("Loading counts and metadata",
                          report.section=report.section,sub=T)
        
        m_a = do.call(read.data.method,
                      c(
                        list(taxa.level=taxa.level),
                        read.data.task
                      )
        )
        
        m_a = get.taxa.meta.aggr(m_a)
        
        report$add.p(paste("After aggregating/subsetting, sample count is:",nrow(m_a$count)))
        
        ## wrap in list to test variable itself rather than its elements
        if(!is.na(list(count.filter.sample.options))) {
          m_a = report.count.filter.m_a(m_a,
                                        count.filter.options=count.filter.sample.options,
                                        "sample")
        }
        
        export.taxa.meta(m_a,
                         label=label,
                         descr=paste(descr,"After initial filtering",sep="."),
                         row.proportions=T,
                         row.names=F)
        
        
        report$pop.section() # end loading counts
        
        if(do.summary.meta) {
          summary.meta.method(m_a)
          do.call(report.sample.count.summary,c(
            list(m_a),
            summary.meta.task
          )
          )
          ##only do summary once
          do.summary.meta = F
        }
        
        res.tests = NULL
        
        if (do.tests) {
          
          ## modify a copy
          test.counts.task.call = test.counts.task
          test.counts.task.call$do.select.samples = (taxa.level %in% test.counts.task.call$do.select.samples)
          test.counts.task.call$do.divrich = (taxa.level %in% test.counts.task.call$do.divrich)
          test.counts.task.call$do.extra.method = (taxa.level %in% test.counts.task.call$do.extra.method)
          test.counts.task.call$do.aggr.after.norm = (taxa.level %in% test.counts.task.call$do.aggr.after.norm)
          
          res.tests = tryCatchAndWarn(
            do.call(test.counts.project,
                    c(
                      list(
                        m_a=m_a,
                        label=label
                      ),
                      test.counts.task.call
                    )
            )
          )          
        }
        
        res.level$res.tests = res.tests
        
        report$pop.section() #4 }
        res.taxa.levels[[as.character(taxa.level)]] = res.level
      }
      report$pop.section() #3 }
      res.taxa.levels
    })
    
    report$pop.section() #2 }
    
    res.task
  })
  
  report$pop.section() #1 }
  
  stop.cluster.project(cl)
  
  res = new_mgsatres(res.tasks=res.tasks)
  
  return (res)
}


##Delta method to approximate mean and variance of
##a function of random variable
##(http://en.wikipedia.org/wiki/Taylor_expansions_for_the_moments_of_functions_of_random_variables)
##See ref for the caviats
mom.f <- function(var.x,mean.x,f,f.d1,f.d2) {
  mean.f = f(mean.x) + f.d2(mean.x)/2*var.x
  var.f = f.d1(mean.x)**2*var.x
  c(mean.f=mean.f,var.f=var.f,sd.f=var.f**0.5)
}

cohens.d.from.mom <- function(mean.gr,var.gr,n.gr) {
  lx <- n.gr[1] - 1
  ly <- n.gr[2] - 1
  md  <- abs(mean.gr[1] - mean.gr[2])        ## mean difference (numerator)
  csd <- lx * var.gr[1] + ly * var.gr[2]
  csd <- csd/(lx + ly)
  csd <- sqrt(csd)                     ## common sd computation
  
  md/csd                        ## cohen's d
}

show.distr <- function(x,binwidth=NULL,dens=T) {
  #because aes leaves its expressions unevaluated, we need
  #to bind the value of x as data frame parameter of ggplot
  pl = ggplot(data.frame(x=x),aes(x=x))
  if (dens) {
    pl = pl + geom_histogram(aes(y=..density..),
                             binwidth=binwidth,color="black",fill=NA) +
      geom_density()
  }
  else {
    pl = pl + geom_histogram(binwidth=binwidth,color="black",fill=NA)
  }
  pl
  #stat_density()
  #hist(x,
  #     freq = F,
  #     breaks = "FD",      # For more breaks than the default
  #     col = "darkslategray4", border = "seashell3")
  #lines(density(x),
  #      col = "firebrick2", lwd = 3)
}

show.distr.group <- function(x,group,binwidth=NULL,dens=T) {
  #because aes leaves its expressions unevaluated, we need
  #to bind the value of x as data frame parameter of ggplot
  pl = ggplot(data.frame(x=x,Group=group),aes(x=x,fill=Group))
  if(dens) {
    pl = pl + geom_density(alpha=0.2,aes(y = ..density..))
  }
  else {
    pl = pl + geom_histogram(alpha=0.5,binwidth=binwidth,position = "identity") + 
      geom_density(aes(y = ..count..,color=Group),fill=NA,alpha=1,size=rel(3))
  }
  pl
}

show.trend <- function(meta.data,x,y,group,title="Group trends") {
  df = meta.data[,c(x,y,group)]
  names(df) = c("x","y","group")
  #df$x = as.Date(timeDate(df$x))  
  df$y = boxcox.transform.vec(df$y)$x
  print(summary(lm(y~x,df)))
  smooth_method = "loess" #"lm"
  print(
    ggplot(df, aes(x=x, y=y,color=group)) +
      geom_point() +
      #geom_line(alpha=0.3, linetype=3) + 
      #geom_smooth(aes(group=group,color=group), method='lm', formula=y~x+I(x^2)+I(x^3)) + 
      stat_smooth(method=smooth_method, se = T,degree=1) + 
      #scale_x_date() +
      labs(title=title)
  )
}

make.sample.summaries <- function(m_a.abs) {
  x = data.frame(
    count.sum = rowSums(m_a.abs$count)
  )
  return (list(count=x,attr=m_a.abs$attr))
}


show.sample.summaries.meta <- function(m_a,
                                       x.var,
                                       group.var,
                                       summ.names=NULL,
                                       value.name="Y",
                                       trans="ident",
                                       vars.descr="sample summary properties") {
  if(is.null(summ.names)) {
    summ.names = names(m_a$count)
  }
  show.feature.meta(m_a=m_a,
                    feature.names=summ.names,
                    x.var=x.var,
                    group.var=group.var,
                    value.name=value.name,
                    trans=trans,
                    vars.descr=vars.descr) 
  
}


show.feature.meta <- function(m_a,
                              feature.names,
                              x.var,
                              group.var,
                              value.name="Abundance",
                              trans="boxcox",
                              vars.descr="Abundances") {
  
  
  count = m_a$count[,feature.names,drop=F]
  
  count = switch(trans,
                 boxcox=norm.boxcox(count),
                 ihs=ihs(count,1),
                 ident=count,
                 binary=(count > 0))  
  
  if(trans != "ident") {
    trans.msg = paste("(transformed with",trans,")")
  }
  else {
    trans.msg = ""
  }
  
  if(is.null(group.var)) {
    group.var.msg = ""
  }
  else {
    group.var.msg = paste("split by",group.var)
  }
  
  report.section = report$add.header(paste(vars.descr,
                                           trans.msg,
                                           "as a function of",
                                           x.var,
                                           group.var.msg),
                                     section.action="push")
  
  id.vars = unique(c(x.var,group.var))
  dat = cbind(m_a$attr[,id.vars,drop=F],count)
  dat = melt.abund.meta(dat,id.vars=id.vars,attr.names=id.vars,value.name=value.name)
  dat$.x.var = dat[,x.var]
  smooth_method = "loess" #"lm"
  for(feature.name in feature.names) {
    dat.feature = dat[dat$feature==feature.name,]
    x.var.range = max(dat.feature$.x.var) - min(dat.feature$.x.var)
    pl = ggplot(dat.feature, aes_string(x=x.var, y=value.name,color=group.var)) +
      geom_point() +
      #geom_line(alpha=0.3, linetype=3) + 
      #geom_smooth(aes(group=group,color=group), method='lm', formula=y~x+I(x^2)+I(x^3)) + 
      stat_smooth(method=smooth_method, size = 1, se = T,method.args=list(degree=1)) +
      ylab(paste(value.name,"of",feature.name))
    
    gr.range = plyr::ddply(dat.feature,group.var,plyr::summarise,
                           range = max(.x.var) - min(.x.var))
    for (group in gr.range[gr.range$range<.Machine$double.eps,group.var]) {
      pl = pl + stat_summary(aes_string(x=x.var, y=value.name,color=group.var),
                             fun.data = "mean_cl_boot", geom = "crossbar", width=x.var.range/20,
                             data = dat.feature[dat.feature[,group.var]==group,])
    }
    #scale_x_date() +
    #labs(title=title)+
    #facet_wrap(~feature,scales="free")
    report$add(pl,
               caption=paste("Value",
                             trans.msg,
                             "of",
                             feature.name,
                             "as a function of",
                             x.var,
                             group.var.msg)
    )
  }
}


##Code adapted from Jyoti Shankar. Find best alpha through cross-validation as
##per the recipe from cv.glmnet help page.
##Extra argument q=dfmax+1 to make this function compatible with 
##the parameter list for stabs::stabsel (see how this is used in
##stabsel.report to contrain the number of selected variables in the model).
cv.glmnet.alpha <- function(y, x, family, q=NULL, seed=NULL,  n.cvs = 400, standardize=T, foreach.errorhandling = "remove", ...) {
  ##on a single fold set picking alpha is extremely unstable - returned alpha is all
  ##over the place from one invocation to another. We run it on multiple sets
  ##of splits, average the error curve (see cv.glmnet help page), and pick the best
  ##parameters from the average.
  # This seed is taken from LASSO's example. 
  if (!is.null(seed) ) { set.seed(seed) }
  # Setting the number of folds to the number of samples (leave one out)
  #is not recommended by cv.glmnet help page
  numfolds <- min(30,max(3,round(nrow(x)/30)))
  # Grid for alpha crossvalidation
  alphas <- c(0.001, 0.005, 0.01, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.7, 0.9, 0.95, 0.99, 0.999,1.0)
  
  foldids <- foreach(seq(n.cvs),.combine=rbind) %do% sample(rep(seq(numfolds),length=nrow(x)))
  # Go through alpha grid
  # Run crossvalidation for lambda.
  # Each model for each alpha is run by a parallel core and put into a list called lassomodels
  if(!is.null(q)) {
    dfmax = q - 1
  }
  else {
    dfmax = ncol(x) + 1
  }
  
  .errorhandling = foreach.errorhandling
  if(getOption("try.debug",F)) {
    #.errorhandling = "stop"
  }
  
  cv.res <- foreach(alpha=alphas,
                    .packages=c("glmnet"),
                    .combine=rbind,
                    .errorhandling=.errorhandling) %do% {
                      #re-import required for windows parallelism with doSNOW
                      #library(glmnet)
                      set.seed(seed)
                      # the function finds the best lambda for a given alpha
                      # within each model there is cross-validation happening for lambda for each alpha.
                      # lambda1 = lambda*alpha 
                      # lambda2 = lambda*(1-alpha)
                      cv.glmnet.args = list(x=x, y=y, family=family,
                                            nfolds=numfolds, 
                                            type.measure="deviance", 
                                            standardize=standardize, 
                                            alpha=alpha,
                                            dfmax=dfmax,
                                            ...)
                      ## get a lambda vector here and pass it to all other cv.glmnet calls
                      ## inside the loop below
                      lambda = do.call(cv.glmnet,c(list(foldid=foldids[1,]),cv.glmnet.args))$lambda
                      plus = function(x,y) {
                        ## we do need to align lambdas here now that we pass a pre-computed lambda
                        ## vector, but this shoud not hurt in any case
                        lambda.common = intersect(x$lambda,y$lambda)
                        x = x[x$lambda %in% lambda.common,]
                        x$cvm = x$cvm + y$cvm[y$lambda %in% lambda.common]
                        x
                      }
                      cvmeans = foreach(i.cv = seq(n.cvs),
                                        .combine=plus,
                                        .final=function(x) (x/n.cvs),
                                        .packages=c("glmnet"),
                                        .errorhandling=.errorhandling) %dopar% {
                                          cv.res = do.call(cv.glmnet,c(list(foldid=foldids[i.cv,]),cv.glmnet.args))
                                          data.frame(lambda=cv.res$lambda,cvm=cv.res$cvm)
                                        }
                      data.frame(alpha=alpha,cvmin=min(cvmeans$cvm))
                    }
  alpha = cv.res$alpha[which.min(cv.res$cvmin)]
  list(param=list(alpha=alpha))
}

stability.selection.c060.at <- function (x, fwer, pi_thr = 0.6) 
{
  stopifnot(pi_thr > 0.5, pi_thr < 1)
  if (class(x$fit)[1] == "multnet") {
    p <- dim(x$fit$beta[[1]])[1]
  }
  else {
    p <- dim(x$fit$beta)[1]
  }
  qv <- ceiling(sqrt(fwer * (2 * pi_thr - 1) * p))
  lpos <- which(x$qs > qv)[1]
  if (!is.na(lpos)) {
    stable <- which(x$x[, lpos] >= pi_thr)
    stable.p <- x$x[stable, lpos]
    lpos.order <- order(x$x[, lpos],decreasing = TRUE)
    lpos.order.p <- x$x[lpos.order,lpos]
  }
  else {
    stable <- NULL
    lpos.order <- NULL
  }
  mean.p = rowMeans(x$x)
  mean.order = order(mean.p,decreasing = TRUE)
  mean.order.p = mean.p[mean.order]
  #Can use factor(z,levels=labels[mean.order],ordered=T)
  #to order some z by the ranking of variables
  #returned by this method
  out <- list(stable = stable, 
              stable.p = stable.p,
              lpos.order = lpos.order,
              lpos.order.p = lpos.order.p,
              mean.order = mean.order,
              mean.order.p = mean.order.p,
              labels = rownames(x$x),
              lambda = x$fit$lambda[lpos], 
              lpos = lpos, fwer = fwer,
              stab.path=x)
  return(out)
}

plot.stability.selection.c060.at <- function (y, 
                                              annot = TRUE, 
                                              main = "Stability path", 
                                              log.scale = TRUE, 
                                              nvar = 8,
                                              xvar="lambda",
                                              rank="lpos") 
{
  x = y$stab.path
  prob = t(x$x)
  nzeros <- which(colSums(prob) != 0)  
  p = ncol(prob)
  
  rank.order = switch(rank,lpos=y$lpos.order,mean=y$mean.order)
  selection <- rep("unselected", ncol(prob))
  names(selection) = colnames(prob)
  labels = names(selection)
  rank.sel = paste("Top",nvar,"by rank")
  fwer.sel = paste("FWER <",y$fwer)
  selection[rank.order[1:nvar]] = rank.sel
  if(!is.null(y$stable)) {
    selection[y$stable] = fwer.sel
  }
  
  xv <- switch(xvar, fraction = x$fit$lambda/max(x$fit$lambda), 
               x$fit$lambda)
  if (xvar == "lambda") {
    if (log.scale == T) {
      xv <- log10(xv)
    }
  }
  data.coef <- melt(data.frame(xvar = xv, prob = t(x$x)), 
                    id="xvar")
  data.coef$selection <- factor(rep(selection, each = length(xv)))
  data.coef$labels <- factor(rep(names(selection), each = length(xv)),
                             levels=labels[rank.order],
                             ordered=T)
  colnames(data.coef) <- c("xvar", "var", "prob", "Selection", 
                           "Variables")
  linetype.values = c(unselected = "dotted")
  linetype.values[rank.sel] = "dashed"
  linetype.values[fwer.sel] = "solid"
  #data.coef$prob = log10(data.coef$prob)
  ##make.global(data.coef)
  d <- ggplot() + 
    #geom_line(aes(x = xvar, 
    #              y = prob)) + 
    geom_line(data=data.coef[data.coef$Selection == "unselected",],
              aes(x = xvar, y = prob, group=var),
              color="black")+    
    geom_line(data=data.coef[data.coef$Selection != "unselected",],
              aes(x = xvar, y = prob, group=var,
                  linetype = Selection, 
                  colour = Variables),
              size=1) +
    scale_linetype_manual(values = linetype.values)+
    labs(x = switch(xvar, fraction = expression(lambda[1]/max[lambda]), 
                    ifelse(log.scale, expression(log[10](lambda)), 
                           expression(lambda[1]))), y = "Selection probabilities") + 
    ggtitle(main) +
    coord_fixed()
  #d <- d + scale_x_reverse()
  if (is.null(labels)) {
    d <- d + theme(legend.position = "none")
  }
  else {
    if (length(labels[nzeros]) != length(nzeros)) {
      d <- d + theme(legend.position = "none")
    }
  }
  return (d)
}


new_deseq2 <- function(...) {
  x = new_mgsatres(...)
  class(x) <- append(class(x),"deseq2",0)
  return(x)
}

get.deseq2.result.description <- function(x) {
  paste(capture.output(show(x))[1:1],collapse=";")
}

deseq2.report.results <- function(res,formula.rhs,result.task,wrap.vals=T) {
  formula.rhs = paste(formula.rhs,collapse = " ")
  #res.descr = get.deseq2.result.description(res)
  res.df = cbind(feature = rownames(res),as.data.frame(res))
  caption=paste(formula.rhs,
                arg.list.as.str(result.task),
  #              res.descr,
                sep=";")  
  #report$add.printed(summary(res),caption=paste("DESeq2 summary for task:",caption))
  report$add.table(as.data.frame(res.df),
                   caption=paste("DESeq2 results for task:",caption),
                   show.row.names=F,
                   wrap.vals=wrap.vals)
  caption
}

as_data_table_with_rownames <- function(df) {
  id_rownames = rownames(df)
  if(is.null(id_rownames)) {
    stop("Need defined rownames for conversion to data.table")
  }
  df = as.data.table(df)
  df[,.id_rownames:=id_rownames]
  df
}

deseq2.join_results_with_design <- function(dds,res) {
  library(DESeq2)
  library(data.table)
  rd = rowData(dds)
  rownames(rd) = rownames(dds)
  rd = as_data_table_with_rownames(rd)
  res = as_data_table_with_rownames(res)
  res = res[rd,on=".id_rownames"]
  res = res[order(abs(stat),decreasing = T,na.last = T)]
  res = setDF(res)
  rownames(res) = res[,".id_rownames"]
  res[[".id_rownames"]] = NULL
  res
}

#' m_a can be m_a or DESeqDataSet type
#' If m_a is DESeqDataSet, and formula.rhs is not NULL,
#' if will override formula in the DESeqDataSet
deseq2.report <- function(m_a,
                          formula.rhs,
                          test.task=list(),
                          result.tasks=list(list()),
                          round.to.int=T,
                          alpha=0.05,
                          do.report = T
) {
  library(DESeq2)
  library(foreach)
  if(do.report) {
    report.section = report$add.header("DESeq2 tests and data normalization",section.action="push")
    report$add.package.citation("DESeq2")
  }
  if(!inherits(m_a,"DESeqDataSet")) {
    dds = as.dds.m_a(m_a=m_a,
                     formula.rhs=formula.rhs,
                     force.lib.size=T,
                     round.to.int=round.to.int)
  }
  else {
    dds = m_a
    if(!is.null(formula.rhs)) {
      design(dds) = as.formula(sprintf("~ %s",formula.rhs))
    }
  }
  dds = do.call(DESeq,
                c(list(object=dds),
                  test.task)
  )
  formula.rhs = as.character(design(dds))[2] #cuts the ~
  res.all = foreach(result.task=result.tasks) %do% {
    if(is.null(result.task$alpha)) {
      result.task$alpha = alpha
    }
    res = do.call(results,
                  c(list(object=dds),
                    result.task)
    )
    res = deseq2.join_results_with_design(dds,res)    
    if(do.report) {
      deseq2.report.results(res=res,formula.rhs=formula.rhs,result.task=result.task)      
    }
    res
  }
  return (new_deseq2(dds=dds,results=res.all,result.tasks=result.tasks,formula.rhs=formula.rhs))
}

glmnet.stabpath.c060.report <- function(m_a,
                                        resp.attr,
                                        family="gaussian",
                                        steps=600,
                                        weakness=0.8,
                                        standardize.count=T,
                                        transform.count="ident",
                                        pred.attr=NULL,
                                        fwer.alpha=0.05) {
  
  report$add.header(paste(
    "Glmnet stability path analysis for response (",
    resp.attr,
    ")"
  )
  )
  
  m = m_a$count
  count = switch(transform.count,
                 boxcox=norm.boxcox(m),
                 ihs=ihs(m,1),
                 ident=m,
                 binary=(m > 0),
                 clr=norm.clr(m))
  
  
  attr = m_a$attr
  if(standardize.count) {
    count = decostand(count,method="standardize",MARGIN=2)
  }
  
  cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
  registerDoSNOW(cl)  
  cv.res = cv.glmnet.alpha(attr[,resp.attr],count,family=family,standardize=F)
  stopCluster(cl)
  
  penalty.alpha = cv.res$alpha
  #alpha = 0.8
  stab.res.c060 = stabpath(attr[,resp.attr],count,weakness=weakness,
                           family=family,steps=steps,
                           alpha=penalty.alpha,standardize=F)
  #stab.res.c060 = stabpath(m_a$attr$A1C[m_a$attr$T1D=="T1D"],
  #                               count[m_a$attr$T1D=="T1D",],
  #                               weakness=0.9,
  #                               family="gaussian",steps=600,
  #                               alpha=penalty.alpha,standardize=F)
  
  fwer = fwer.alpha
  pi_thr = 0.6
  stab.feat.c060 = stability.selection.c060.at(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
  #stab.feat.c060 = stability.selection(stab.res.c060,fwer=fwer,pi_thr=pi_thr)
  stab.feat.c060$penalty.alpha = penalty.alpha
  
  pl.stab = tryCatchAndWarn({
    #p = plot.stability.selection.c060.at(stab.feat.c060,rank="mean")
    p = plot.stability.selection.c060.at(stab.feat.c060,xvar="fraction",rank="lpos")
    #ggsave(stab.path.file)
    p
  })
  
  #report$add.printed(stab.res.c060$fit,
  #                   caption=paste(
  #                     "Glmnet stability path analysis for response (",
  #                     resp.attr,
  #                     ")"
  #                     )
  #)
  
  report$add.printed(summary(m_a$attr[,resp.attr]),
                     caption="Summary of response variable")
  
  report$add.package.citation("c060")
  report$add.descr("This multivariate feature selection method implements 
                  stability selection procedure by Meinshausen and Buehlmann (2010) 
                  The features (e.g. taxonomic features)
                   are ranked according to their probability to be selected
                   by models built on multiple random subsamples of the input dataset.")
  
  report$add.descr(paste("Alpha penalty parameter was chosen in cross-validation as",
                         penalty.alpha))
  
  report$add.vector(stab.feat.c060$stable,
                    caption="Features that passed FWER in stability analysis")
  
  if(!is.null(pl.stab)) {
    report$add(pl.stab,caption="Stability path. Probability of each variable 
               to be included in the model as a function of L1 regularization
               strength. Paths for top ranked variables are colored. The variables
               (if any) that passed family wide error rate (FWER) cutoff are
               plotted as solid lines.")
  }
  
  return (stab.feat.c060)
}


stabsel.report <- function(m_a,
                           resp.attr,
                           fitfun="glmnet.lasso",
                           parfitfun="cv.glmnet.alpha",
                           args.fitfun = list(
                             family="gaussian",
                             standardize=T                                        
                           ),
                           args.stabsel = list(
                             PFER=0.05,
                             sampling.type="SS",
                             assumption="r-concave",
                             q=NULL,
                             B=200
                           )
) {
  
  require(stabs)
  report$add.header(paste(
    "Stability selection analysis for response (",
    paste(resp.attr,collapse=","),
    ")"
  )
  )
  
  resp = m_a$attr[,resp.attr]
  if(is.factor(resp)) {
    resp = factor(resp) #remove levels without observations
  }
  count = m_a$count
  
  if(is.null(args.stabsel$q)) {
    #MB's paper
    q = ceiling(sqrt(0.8*nrow(count)))
    ##alternative is nrow(count)/log(ncol(count))
    ##as the condition on the size of support for Theta
    ##in refs [CT07, BRT09, BvdG11] from
    ##{Confidence Intervals and Hypothesis Testing for
    ##High-Dimensional Regression,
    ##Adel Javanmard and Andrea Montanari}
    ##http://web.stanford.edu/~montanar/sslasso/
    
    q = min(q,ncol(count))
    args.stabsel$q = q
  }
  else {
    q = args.stabsel$q
  }
  
  param.fit.res = do.call(parfitfun,
                          c(list(x=count,y=resp,q=q),
                            args.fitfun
                          )
  )
  param.fit = param.fit.res$param
  args.fitfun = c(args.fitfun,param.fit)
  stab.res = do.call(stabsel,c(
    list(x=count,y=resp,
         fitfun=fitfun,
         args.fitfun=args.fitfun),
    args.stabsel
  )
  )
  
  report$add.printed(summary(resp),
                     caption=paste("Summary of response variable",resp.attr))
  
  report$add.package.citation("stabs")
  report$add.descr("This multivariate feature selection method implements 
                  stability selection procedure by Meinshausen and Buehlmann (2010) 
                  and the improved error bounds by Shah and Samworth (2013). 
                  The features (e.g. taxonomic features)
                   are ranked according to their probability to be selected
                   by models built on multiple random subsamples of the input dataset.")
  
  report$add.descr(paste("Base selection method parameters that were chosen based on
                         cross-validation are:",
                         arg.list.as.str(param.fit)))
  
  report$add.descr(paste("All base selection method parameters are:",
                         arg.list.as.str(args.fitfun)))  
  
  report$add.descr(paste("Stability selection method parameters are:",
                         arg.list.as.str(args.stabsel)))
  
  #report$add.printed(stab.res)
  
  cutoff.descr = sprintf("Probability cutoff=%s corresponds to per family error rate PFER=%s",
                         stab.res$cutoff,
                         stab.res$PFER)  
  
  report$add.vector(get.feature.ranking(stab.res,only.names=F)$ranked,
                    name="Prob(selection)",
                    caption=paste("Selection probability for the variables.",
                                  cutoff.descr)
  )
  
  report$add(plot(stab.res, main = NULL, type = "maxsel",
                  ymargin = 16, np = 20),
             caption=paste("Selection probability for the top ranked variables.",
                           cutoff.descr,
                           "(vertical line).")
  )
  
  return (stab.res)
}

num.levels <- function(x) {
  length(levels(factor(x)))
}

num.levels.data.frame <- function(x) {
  nrow(as.data.frame(table(x)))
}


genesel.stability.report <- function(m_a,group.attr,
                                     genesel.param=list(),
                                     do.nmds=F,
                                     do.plot.profiles = T,
                                     plot.profiles.task=list(),
                                     norm.count.task=NULL) {
  
  ## m_a$count passed here should be normalized for library size, because we perform Wilcox test
  ## inside, and could find false significance due to different depth of sequencing
  
  report$add.header("GeneSelector stability ranking")
  if(genesel.param$type %in% c("unpaired","paired") &&
     num.levels(m_a$attr[,group.attr]) != 2) {
    report$add.descr(sprintf("GeneSelector analysis is skipped 
                             because grouping factor %s does not have two levels",group.attr))
    return(NULL)
  }
  report$add.package.citation("GeneSelector")
  if(!is.null(norm.count.task)) {
    m_a <- norm.count.report(m_a,
                             descr="GeneSelector",
                             norm.count.task=norm.count.task)
  }  
  report$add.descr(sprintf("Wilcoxon test (rank-sum for independent samples and signed-rank for paired samples) 
                   is applied to each feature (feature, gene) on random
                   subsamples of the data. Consensus ranking is found with a
                   Monte Carlo procedure ((method AggregateMC in GeneSelector package). 
                   features ordered according to the consensus ranking
                   are returned, along with the p-values, statistic and effect size 
                   computed on the full
                   original dataset. In a special case when no replications are requested,
                   features are ordered by the adjuested p-value. 
                   P-values are reported with and without the 
                   multiple testing correction of Benjamini & Hochberg. The effect sizes
                   for Wilcoxon tests are reported as: common-language effect
                   size (proportion of pairs where observations from the second group
                   are larger than observations from the first group; no effect
                   corresponds to 0.5); rank-biserial
                   correlation (common language effect size minus its complement; no
                   effect corresponds to 0; range is [-1;1]) and
                   absolute value of r (as defined in Cohen, J. (1988). Statistical power 
                   analysis for the behavioral sciences (2nd ed.). Hillsdale, NJ: Erlbaum.).
                   For paired samples, when calculating the common language effect size,
                   only paired observations are used, and one half of the number of ties is 
                   added to the numerator (Grissom, R. J., and J. J. Kim. \"Effect Sizes for Research: Univariate 
                   and Multivariate Applications, 2nd Edn New York.\" NY: Taylor and Francis (2012)).
                   Logarithm in base 2 of the fold change (l2fc) is also reported if requested.
                   For independent samples, the fold change is computed between the sample means of
                   the groups (last to first). For paired samples - as the sample median of the logfold change
                   in each matched pair of observations."))
  
  res.genesel.stability = do.call(
    genesel.stability,
    c(list(m_a,group.attr=group.attr),
      genesel.param)
  )
  
  ord = res.genesel.stability$stab_feat$name
  feature.order = list(list(ord=ord,ord_descr="GeneSelector paired test ranking",sfx="gsp"))
  
  report$add.descr(sprintf("Stability selection parameters are: %s",arg.list.as.str(genesel.param)))
  
  report$add.printed(summary(m_a$attr[,group.attr]),
                     caption=paste("Summary of response variable (unpaired samples)",group.attr))
  
  caption.descr = ""
  if(genesel.param$type=="paired" && !is.null(genesel.param$block.attr)) {
    caption.descr = sprintf("Samples are paired according to attribute %s, resulting in %s pairs.",
                            genesel.param$block.attr,res.genesel.stability$n.samp)
  }
  
  caption.descr = sprintf("%s When fold change or difference is computed, this is done as '%s'.",caption.descr,
                          paste(res.genesel.stability$levels.last.first,collapse=" by ")
  )
  
  report$add.table(res.genesel.stability$stab_feat,
                   caption=
                     paste(sprintf("GeneSelector stability ranking for response %s.",group.attr),
                           caption.descr)
  )
  
  if(do.plot.profiles && genesel.param$type!="unpaired") {
    
    id.vars.list = plot.profiles.task$id.vars.list
    ## add group.attr if it is not already where in
    ## order to get empty element in the next step
    if(! group.attr %in% id.vars.list) {
      id.vars.list = c(group.attr,id.vars.list)
    }
    ## remove group.attr from all elements of id.vars.list
    plot.profiles.task$id.vars.list = sapply(id.vars.list,function(y) y[y != group.attr])
    ## remove group.attr from feature.meta.x.vars
    plot.profiles.task$feature.meta.x.vars = 
      plot.profiles.task$feature.meta.x.vars[plot.profiles.task$feature.meta.x.vars != group.attr]
    plot.profiles.task = within(plot.profiles.task, {
      do.profile=T
    })
    
    if(!is.null(res.genesel.stability$m_a.contrasts)) {
      m_a.c = res.genesel.stability$m_a.contrasts
      tryCatchAndWarn({
        
        plot.profiles.task$value.name = "abundance.diff"
        plot.profiles.task$feature.descr = paste("Abundance difference between paired samples.",
                                                 caption.descr)
        
        do.call(plot.profiles,
                c(list(m_a=m_a.c,
                       feature.order=feature.order),
                  plot.profiles.task
                )
        )
        
      })
      
    }
    
    if(!is.null(res.genesel.stability$m_a.lfc.paired)) {
      m_a.c = res.genesel.stability$m_a.lfc.paired
      tryCatchAndWarn({
        
        plot.profiles.task$value.name = "l2fc"
        plot.profiles.task$feature.descr = paste("Log2 fold change in abundance between paired samples.",
                                                 caption.descr)
        
        if(is.null(plot.profiles.task$show.profile.task)) {
          plot.profiles.task$show.profile.task = list()
        }
        plot.profiles.task$show.profile.task$stat_summary.fun.y="median"
        do.call(plot.profiles,
                c(list(m_a=m_a.c,
                       feature.order=feature.order),
                  plot.profiles.task
                )
        )
        
      })
      
    }
    
  }
  
  if(do.nmds) {
    report$add.package.citation("vegan")
    m_a.mds = m_a
    m_a.mds$count = m_a.mds$count[,res.genesel.stability$stab_feat$name]
    report$add(
      plot.features.mds(m_a.mds,sample.group=m_a$attr[,group.attr]),
      caption="metaMDS plot of features selected by GeneSelector. 'x' marks features, 'o' marks samples"
    )
  }
  
  return (res.genesel.stability)
  
}

test.counts.glmer.report <- function(m_a,
                                     tasks,
                                     alpha=0.05) {
  descr.tpl = "
  Mixed model analysis of count data.
  The binomial family
  is used to build a set of univariate models, with each
  model describing the counts for one feature.
  We add a random effect for each sample to account
  for the overdispersion;
  P-values are estimated from the model under a null hypothesis
  of zero coefficients and a two-sided alternative. 
  Benjamini & Hochberg (1995) method is used 
  for multiple testing correction, and the significant features
  are reported.
  %s
  "
  
  report$add.header("Binomial mixed model analysis")
  report$add.package.citation("lme4")
  res = lapply(tasks,function(task) {
    with(task,{
      res.glmer = test.counts.glmer(m_a,alpha=alpha,
                                    formula_rhs=formula.rhs,
                                    linfct=linfct)
      
      descr = sprintf(descr.tpl,descr.extra)
      report$add.descr(descr)
      report.counts.glmer(report,res.glmer)
      res.glmer
    })
  })
  return(res)
}

test.counts.adonis.report <- function(m_a,
                                      tasks,
                                      n.perm=4000,
                                      dist.metr="bray",
                                      col.trans="range",
                                      data.descr="proportions of counts",
                                      norm.count.task=NULL) {
  library(vegan)
  report$add.header(paste("PermANOVA (adonis) analysis of ",data.descr))
  report$add.package.citation("vegan")  
  is.dist = inherits(m_a$count,"dist")
  if(!is.dist && !is.null(norm.count.task)) {
    m_a <- norm.count.report(m_a,
                             descr="Adonis",
                             norm.count.task=norm.count.task)
  }
  
  ##Negative values break bray-curtis and jaccard distances; we standardize to "range" to reduce
  ##the influence of very abundant species:
  count = m_a$count
  if(!is.dist && !is.null(col.trans) && col.trans != "ident") {
    count = decostand(count,method=col.trans,MARGIN=2)
    col.trans.descr = sprintf(" Profile columns are normalized with %s method of decostand function.",col.trans)
  }
  else {
    col.trans.descr = ""
  }
  
  if(!is.dist) {
    dist.metr.descr = sprintf(" Dissimilarity index is %s.",dist.metr)
  }
  else {
    dist.metr.descr = " Using supplied distance matrix."
  }
  
  #print(ad.res)
  report$add.descr(sprintf("Non-parametric multivariate test for association between
                           %s and meta-data variables.%s%s",
                           data.descr,
                           col.trans.descr,
                           dist.metr.descr))
  
  res = lapply(tasks,function(task) {
    strata = task$strata #implicitely defined here even if undefined in task
    if(is.null(strata)) {
      strata.descr = ""
    }
    else {
      strata.descr = paste("with strata = ",strata)
    }
    with(task,{
      formula_str = paste("count",formula.rhs,sep="~")
      perm = permute::how(nperm=n.perm)
      if(!is.null(strata)) {
        setBlocks(perm) = m_a$attr[,strata]
      }
      ad.res = adonis2(
        as.formula(formula_str),
        data=m_a$attr,
        permutations=perm,
        method=dist.metr)
      
      #report$add.p(pandoc.formula.return(as.formula(formula_str),caption=descr)
      report$add.printed(ad.res,caption=paste(descr,
                                              "with formula",
                                              pandoc.escape.special(formula_str),
                                              strata.descr
      )
      )
      report$add(ad.res,
                 caption=paste(descr,"Adonis summary"))
      ad.res
    })
  })
  
  return (res)
  
}

select.samples.report <- function(m_a,
                                  feature.ranking,
                                  resp.attr,
                                  n.species,
                                  n.samples) {
  
  species.sel=take_first(feature.ranking,
                         n.species)
  select.samples(m_a=m_a,
                 species.sel=species.sel,
                 sample.group.name=resp.attr,
                 n.select=n.samples)
  
}


plot.profiles.abund <- function(m_a,
                                label,
                                plot.profiles.task,                                
                                res.tests=NULL,
                                norm.count.task=NULL,
                                feature.ranking="stabsel") {
  
  if(!is.null(norm.count.task)) {
    m_a <- norm.count.report(m_a,
                             res.tests=res.tests,
                             descr="abundance plots",
                             norm.count.task=norm.count.task)
  }
  
  tryCatchAndWarn({
    
    feature.order = list(list(ord=NULL,ord_descr="average abundance",sfx="ab"))
    if (!is.null(res.tests)) {
      ord=get.feature.ranking(res.tests,feature.ranking)$ranked
      if(!is.null(ord)) {
        feature.order[[2]] = list(ord=ord,
                                  ord_descr=sprintf("Ranking by '%s' method",feature.ranking),
                                  sfx=feature.ranking)
      }
    }
    
    do.call(plot.profiles,
            c(list(m_a=m_a,
                   feature.order=feature.order),
              plot.profiles.task
            )
    )
    
  })
  
}

##pull necessary data (such as DESeq2 normalization) from previously
##done tests into the norm.count.task object
update.norm.count.task <- function(norm.count.task,res.tests=NULL) {
  ## if dds is required but not defined (set to NA)
  if(!is.null(norm.count.task$method.args$dds) && 
     is.na(norm.count.task$method.args$dds)) {
    if(!is.null(res.tests$deseq2$dds)) {
      norm.count.task$method.args$dds = res.tests$deseq2$dds
    }
  }
  return (norm.count.task)  
}

norm.count.report <- function(m_a,norm.count.task=NULL,res.tests=NULL,descr=NULL) {
  
  if(is.null(descr)) {
    descr = ""
  }
  else {
    descr = paste("for",descr)
  }
  
  if(!is.null(norm.count.task)) {
    
    descr = paste("Count normalization method",descr,":",arg.list.as.str(norm.count.task))
    
    norm.count.task = update.norm.count.task(norm.count.task,res.tests=res.tests)
    
    m_a.norm <- do.call(norm.count.m_a,
                        c(list(m_a),
                          norm.count.task)
    )
    
  }
  else {
    
    descr = paste("Counts are not normalized and no features are dropped",descr)
    
    m_a.norm = m_a
    
  }
  
  report$add.descr(descr)
  
  return(m_a.norm)
}

## Either returns default m_a.norm or updates norm.count.task within task if defined
pull.norm.count.task <- function(m_a,m_a.norm,task,res.tests) {
  if(!is.null(task$norm.count.task)) {
    task$norm.count.task = 
      update.norm.count.task(task$norm.count.task,res.tests=res.tests)
    m_a.task = m_a
  }
  else {
    m_a.task = m_a.norm
  }
  return(list(m_a.task=m_a.task,task=task))
}

test.counts.project <- function(m_a,
                                label,
                                do.genesel=T,
                                do.stabsel=T,
                                do.glmer=T,
                                do.adonis=T,
                                do.select.samples=F,
                                do.deseq2=T,
                                do.divrich=T,
                                do.divrich.pre.filter=T,
                                do.divrich.post.filter=T,
                                do.plot.profiles.abund=T,
                                do.heatmap.abund=T,
                                do.ordination=T,
                                do.network.features.combined=T,
                                do.extra.method=F,
                                do.aggr.after.norm=F,
                                count.filter.feature.options=NULL,
                                norm.count.task=NULL,
                                stabsel.task=NULL,
                                genesel.task=NULL,
                                adonis.task=NULL,
                                glmer.task=NULL,
                                select.samples.task=NULL,
                                deseq2.task=NULL,
                                divrich.task=NULL,
                                plot.profiles.task=NULL,
                                plot.profiles.abund.task=NULL,
                                heatmap.abund.task=NULL,
                                heatmap.combined.task=NULL,
                                ordination.task=NULL,
                                network.features.combined.task=NULL,
                                alpha=0.05,
                                do.return.data=T,
                                feature.ranking="stabsel",
                                extra.method.task = NULL,
                                aggr.after.norm.task = NULL
) {
  
  report.section = report$add.header("Data analysis",section.action="push")
  
  res = new_mgsatres()
  
  ##drop feature "other"
  m_a.abs = norm.count.m_a(m_a,method="ident")
  
  if(!is.null(genesel.task)) {
    if(is.null(genesel.task$plot.profiles.task)) {
      genesel.task$plot.profiles.task=plot.profiles.task
    }
  }
  
  if(do.divrich) {
    if(is.null(divrich.task$beta.task$adonis.task)) {
      if(!is.null(divrich.task$beta.task)) 
      { divrich.task$beta.task$adonis.task = adonis.task }
      if(!do.adonis) {
        divrich.task$beta.task$adonis.task = NULL
      }
    }
    if(is.null(divrich.task$counts.genesel.task)) {
      divrich.task$counts.genesel.task=genesel.task
    }
    if(!do.genesel) {
      divrich.task$counts.genesel.task = NULL
    }
  }
  if(do.divrich && do.divrich.pre.filter) {
    tryCatchAndWarn({ 
      res$divrich <- do.call(mgsat.divrich.report,
                             c(list(m_a=m_a.abs,
                                    plot.profiles.task=plot.profiles.task,
                                    extra.header="Before count filtering"),
                               divrich.task)
      )
    })
  }
  
  m_a = report.count.filter.m_a(m_a,
                                count.filter.options=count.filter.feature.options,
                                "Filtering features")
  export.taxa.meta(m_a,
                   label=label,
                   descr="After final feature filtering",
                   row.proportions=T,
                   row.names=F)    
  
  m_a.abs = norm.count.m_a(m_a,method="ident")
  
  if(is.null(count.filter.feature.options) || length(count.filter.feature.options)==0) {
    divrich.post.filter=F
  }
  
  #make.global(m_a)
  
  if(do.divrich && do.divrich.post.filter) {
    tryCatchAndWarn({ 
      res$divrich <- do.call(mgsat.divrich.report,
                             c(list(m_a=m_a.abs,
                                    plot.profiles.task=plot.profiles.task,
                                    extra.header="After count filtering"),
                               divrich.task)
      )
    })
  }
  
  if(do.deseq2) {
    tryCatchAndWarn({ 
      ##drop "other" category with norm method "ident"
      if(is.null(deseq2.task$alpha)) {
        deseq2.task$alpha = alpha
      }
      res$deseq2 = do.call(deseq2.report,c(list(m_a=m_a.abs),deseq2.task))
    })
  }
  
  report$add.header("Default transformations for further data analysis",section.action="push")
  report$add.descr("Specific methods can override these and use their own normalization.")
  ## this is done after an optional call to deseq2 in case the norm.method
  ## wants deseq2 normalization
  
  m_a.norm <- norm.count.report(m_a,
                                res.tests=res,
                                descr="data analysis (unless modified by specific methods)",
                                norm.count.task)
  if(do.aggr.after.norm && !is.null(aggr.after.norm.task)) {
    aggr.after.norm.func = aggr.after.norm.task$func
    aggr.after.norm.task$func = NULL
    res.aan = do.call(aggr.after.norm.func,
                      c(list(m_a=m_a,
                             m_a.norm=m_a.norm,
                             m_a.abs=m_a.abs,
                             res.tests=res),
                        aggr.after.norm.task
                      )
    )
    m_a = res.aan$m_a
    m_a.norm = res.aan$m_a.norm
    m_a.abs = res.aan$m_a.abs
    
    export.taxa.meta(m_a,
                     label=label,
                     descr=paste("Raw after default transformations",sep="."),
                     row.proportions=T,
                     row.names=F)
  }
  
  export.taxa.meta(m_a.norm,
                   label=label,
                   descr=paste("Normalized after default transformations",sep="."),
                   row.proportions=F,
                   row.names=F)
  
  report$pop.section()
  
  if(do.genesel) {
    genesel.norm.t = pull.norm.count.task(m_a=m_a,m_a.norm=m_a.norm,
                                          task=genesel.task,res.tests=res)
    tryCatchAndWarn({ 
      res$genesel = do.call(genesel.stability.report,
                            c(list(m_a=genesel.norm.t$m_a.task),
                              genesel.norm.t$task))
    })
  }
  
  #DEBUG:
  #make.global()
  
  if(do.stabsel) {
    tryCatchAndWarn({ 
      res$stabsel = do.call(stabsel.report,c(list(m_a=m_a.norm),
                                             stabsel.task
      )
      )
    })
  }
  
  if(do.glmer) {  
    tryCatchAndWarn({ 
      res$glmer = do.call(test.counts.glmer.report,
                          c(
                            list(m_a=m_a.abs,
                                 alpha=alpha),
                            glmer.task
                          )
      )
    })
  }
  
  if(do.adonis) {
    tryCatchAndWarn({ 
      adonis.norm.t = pull.norm.count.task(m_a=m_a,m_a.norm=m_a.norm,
                                           task=adonis.task,res.tests=res)
      res$adonis = do.call(test.counts.adonis.report,
                           c(
                             list(m_a=adonis.norm.t$m_a.task),
                             adonis.norm.t$task
                           )
      )
    })
  }
  
  if( do.select.samples ) {
    
    tryCatchAndWarn({ 
      do.call(select.samples.report,
              c(
                list(m_a=m_a.norm,
                     feature.ranking=get.feature.ranking(res,feature.ranking)$ranked
                ),
                select.samples.task
              )
      )
    })
    
  }
  
  if( do.plot.profiles.abund ) {
    
    tryCatchAndWarn({ 
      
      plot.profiles.abund.norm.t = pull.norm.count.task(m_a=m_a,m_a.norm=m_a.norm,
                                                        task=plot.profiles.abund.task,res.tests=res)
      
      do.call(plot.profiles.abund,
              c(
                list(m_a=plot.profiles.abund.norm.t$m_a.task,
                     label=label,
                     res.tests=res,
                     plot.profiles.task=plot.profiles.task,
                     feature.ranking=feature.ranking),
                plot.profiles.abund.norm.t$task
              )
      )
    })
  }
  
  if (do.heatmap.abund) {
    
    if(F) tryCatchAndWarn({
      do.call(heatmap.counts,
              c(list(m_a=m_a.norm),
                heatmap.abund.task
              )
      )
      
    })
    
    tryCatchAndWarn({
      do.call(heatmap.combined.report,
              c(list(m_a=m_a,
                     m_a.norm=m_a.norm,
                     res.tests=res),
                heatmap.combined.task
              )
      )
      
    })
    
  }
  
  if (do.ordination) {
    tryCatchAndWarn({
      do.call(ordination.report,
              c(list(m_a=m_a.norm,res=res),
                ordination.task
              )
      )
      
    })
  }  
  
  if (do.network.features.combined) {
    
    tryCatchAndWarn({
      do.call(network.features.combined.report,
              c(list(m_a=m_a.abs,res.tests=res),
                network.features.combined.task
              )
      )
      
    })
  }    
  
  if(do.extra.method && !is.null(extra.method.task)) {
    extra.method.func = extra.method.task$func
    extra.method.task$func = NULL
    res$extra.method = do.call(extra.method.func,
                               c(list(m_a=m_a,
                                      m_a.norm=m_a.norm,
                                      res.tests=res),
                                 extra.method.task
                               )
    )
  }
  
  if(do.return.data) {
    res$m_a = m_a
  }
  return (res)
}

## annHeatmap2 from Heatplus, with fixed reordering of labels
annHeatmap2AT <-
  function (x, dendrogram, annotation, cluster, labels, scale = c("row", 
                                                                  "col", "none"), breaks = 256, col = g2r.colors, legend = FALSE) 
  {
    if (!is.matrix(x) | !is.numeric(x)) 
      stop("x must be a numeric matrix")
    nc = ncol(x)
    nr = nrow(x)
    if (nc < 2 | nr < 2) 
      stop("x must have at least two rows/columns")
    def = list(clustfun = hclust, distfun = dist, status = "yes", 
               lwd = 3, dendro = NULL)
    dendrogram = extractArg(dendrogram, def)
    def = list(data = NULL, control = list(), asIs = FALSE, inclRef = TRUE)
    annotation = extractArg(annotation, def)
    def = list(cuth = NULL, grp = NULL, label = NULL, col = RainbowPastel)
    cluster = extractArg(cluster, def)
    def = list(cex = NULL, nrow = 3, side = NULL, labels = NULL, status = "yes")
    labels = extractArg(labels, def)
    if (is.logical(legend)) {
      if (legend) 
        leg = NULL
      else leg = 0
    }
    else {
      if (!(legend %in% 1:4)) 
        stop("invalid value for legend: ", legend)
      else leg = legend
    }
    layout = heatmapLayout(dendrogram, annotation, leg.side = leg)
    x2 = x
    scale = match.arg(scale)
    if (scale == "row") {
      x2 = sweep(x2, 1, rowMeans(x, na.rm = TRUE))
      sd = apply(x2, 1, sd, na.rm = TRUE)
      x2 = sweep(x2, 1, sd, "/")
    }
    else if (scale == "column") {
      x2 = sweep(x2, 2, colMeans(x, na.rm = TRUE))
      sd = apply(x2, 2, sd, na.rm = TRUE)
      x2 = sweep(x2, 2, sd, "/")
    }
    breaks = niceBreaks(range(x2, na.rm = TRUE), breaks)
    col = breakColors(breaks, col)
    dendrogram$Row = within(dendrogram$Row, if (!inherits(dendro, 
                                                          "dendrogram")) {
      dendro = clustfun(distfun(x))
      dendro = reorder(as.dendrogram(dendro), rowMeans(x, na.rm = TRUE))
    })
    dendrogram$Col = within(dendrogram$Col, if (!inherits(dendro, 
                                                          "dendrogram")) {
      dendro = clustfun(distfun(t(x)))
      dendro = reorder(as.dendrogram(dendro), colMeans(x, na.rm = TRUE))
    })
    rowInd = with(dendrogram$Row, if (status != "no") 
      order.dendrogram(dendro)
      else 1:nr)
    colInd = with(dendrogram$Col, if (status != "no") 
      order.dendrogram(dendro)
      else 1:nc)
    x2 = x2[rowInd, colInd]
    labels$Row = within(labels$Row, {
      if (is.null(cex)) 
        cex = 0.2 + 1/log10(nr)
      if (is.null(side)) 
        side = if (is.null(annotation$Row$data)) 
          4
      else 2
      if(status!="no") {
        if (is.null(labels)) 
          labels = rownames(x2)
        else
          ## need to reorder - plot.annHeatmap does not do it
          labels = labels[rowInd]
      }
      else {
        labels = NULL
      }
    })
    labels$Col = within(labels$Col, {
      if (is.null(cex)) 
        cex = 0.2 + 1/log10(nc)
      if (is.null(side)) 
        side = if (is.null(annotation$Col$data)) 
          1
      else 3
      if(status!="no") {
        if (is.null(labels)) 
          labels = colnames(x2)
        else
          ## need to reorder - plot.annHeatmap does not do it
          labels = labels[colInd]
      }
      else {
        labels = NULL
      }
    })
    cluster$Row = within(cluster$Row, if (!is.null(cuth) && (cuth > 
                                                             0)) {
      grp = cutree.dendrogram(dendrogram$Row$dendro, cuth)[rowInd]
    })
    cluster$Col = within(cluster$Col, if (!is.null(cuth) && (cuth > 
                                                             0)) {
      grp = cutree.dendrogram(dendrogram$Col$dendro, cuth)[colInd]
    })
    annotation$Row = within(annotation$Row, {
      data = convAnnData(data, asIs = asIs, inclRef = inclRef)
    })
    annotation$Col = within(annotation$Col, {
      data = convAnnData(data, asIs = asIs, inclRef = inclRef)
    })
    ret = list(data = list(x = x, x2 = x2, rowInd = rowInd, colInd = colInd, 
                           breaks = breaks, col = col), dendrogram = dendrogram, 
               cluster = cluster, annotation = annotation, labels = labels, 
               layout = layout, legend = legend)
    class(ret) = "annHeatmap2AT"
    ret
  }

library(Heatplus)
environment(annHeatmap2AT) <- asNamespace('Heatplus')

plot.annHeatmap2AT <-
  function (x, widths, heights, ...) 
  {
    if (!missing(widths)) 
      x$layout$width = widths
    if (!missing(heights)) 
      x$layout$height = heights
    with(x$layout, layout(plot, width, height, respect = TRUE))
    nc = ncol(x$data$x2)
    nr = nrow(x$data$x2)
    doRlab = !is.null(x$labels$Row$labels)
    doClab = !is.null(x$labels$Col$labels)
    mmar = c(1, 0, 0, 2)
    if (doRlab) 
      mmar[x$labels$Row$side] = x$labels$Row$nrow
    if (doClab) 
      mmar[x$labels$Col$side] = x$labels$Col$nrow
    with(x$data, {
      par(mar = mmar)
      image(1:nc, 1:nr, t(x2), axes = FALSE, xlim = c(0.5, 
                                                      nc + 0.5), ylim = c(0.5, nr + 0.5), xlab = "", ylab = "", 
            col = col, breaks = breaks, ...)
    })
    with(x$labels, {
      if (doRlab) 
        axis(Row$side, 1:nr, las = 2, line = -0.5, tick = 0, 
             labels = Row$labels, cex.axis = Row$cex)
      if (doClab) 
        axis(Col$side, 1:nc, las = 2, line = -0.5, tick = 0, 
             labels = Col$labels, cex.axis = Col$cex)
    })
    with(x$dendrogram$Col, if (status == "yes") {
      par(mar = c(0, mmar[2], 3, mmar[4]))
      cutplot.dendrogram(dendro, h = x$cluster$Col$cuth, cluscol = x$cluster$Col$col, 
                         horiz = FALSE, axes = FALSE, xaxs = "i", leaflab = "none", 
                         lwd = x$dendrogram$Col$lwd)
    })
    with(x$dendrogram$Row, if (status == "yes") {
      par(mar = c(mmar[1], 3, mmar[3], 0))
      cutplot.dendrogram(dendro, h = x$cluster$Row$cuth, cluscol = x$cluster$Row$col, 
                         horiz = TRUE, axes = FALSE, yaxs = "i", leaflab = "none", 
                         lwd = x$dendrogram$Row$lwd)
    })
    if (!is.null(x$annotation$Col$data)) {
      par(mar = c(1, mmar[2], 0, mmar[4]), xaxs = "i", yaxs = "i")
      picketPlot(x$annotation$Col$data[x$data$colInd, , drop = FALSE], 
                 grp = x$cluster$Col$grp, grpcol = x$cluster$Col$col, 
                 control = x$annotation$Col$control, asIs = TRUE)
    }
    if (!is.null(x$annotation$Row$data)) {
      par(mar = c(mmar[1], 0, mmar[3], 1), xaxs = "i", yaxs = "i")
      picketPlot(x$annotation$Row$data[x$data$rowInd, , drop = FALSE], 
                 grp = x$cluster$Row$grp, grpcol = x$cluster$Row$col, 
                 control = x$annotation$Row$control, asIs = TRUE, 
                 horizontal = FALSE)
    }
    if (x$legend) {
      if (x$layout$legend.side %in% c(1, 3)) {
        par(mar = c(2, mmar[2] + 2, 2, mmar[4] + 2))
      }
      else {
        par(mar = c(mmar[1] + 2, 2, mmar[3] + 2, 2))
      }
      doLegend(x$data$breaks, col = x$data$col, x$layout$legend.side)
    }
    invisible(x)
  }

environment(plot.annHeatmap2AT) <- asNamespace('Heatplus')

## taxa count columns in meta.data must be already sorted by some abundance metrics
heatmap.counts <- function(m_a,
                           attr.annot.names=NULL,
                           attr.row.labels=NULL,
                           caption="Heatmap",
                           max.species.show=30,
                           do.row.clust=T,
                           do.col.clust=T,
                           stand.clust=NULL,
                           trans.clust=NULL,
                           dist.metr="euclidean",
                           stand.show="range",
                           trans.show=NULL,
                           attr.order=NULL,
                           agglo.fun.order=sum,
                           cluster.row.cuth=2) {
  
  library(RColorBrewer)
  library(vegan)
  
  ##permute samples to make sure that our dendrogram
  ##clustering is not influenced by the original order
  perm.ind = sample(nrow(m_a$count))
  
  count.src = m_a$count[perm.ind,]
  count = count.src
  
  if(!is.null(trans.clust)) {
    count = do.call(trans.clust,list(count))
  }
  if(!is.null(stand.clust)) {
    count = decostand(count,method=stand.clust,MARGIN=2)
  }
  
  attr = m_a$attr[perm.ind,]
  if(do.row.clust) {
    data.dist.samp <- vegdist(count, method = dist.metr)
    row.clus <- hclust(data.dist.samp, "ward.D2")
    row.dendro = as.dendrogram(row.clus)
    
    if(is.null(attr.order)) {
      wgts = rowMeans(count, na.rm = TRUE)
    }
    else {
      wgts = attr[,c(attr.order)]
    }
    row.dendro = reorder(row.dendro, wgts, agglo.FUN=agglo.fun.order)
    #row.ind = order.dendrogram(row.dendro)
  }
  else {
    row.dendro=NA
  }
  if(!is.null(attr.annot.names)) {
    attr.annot = attr[,attr.annot.names,drop=F]
  }
  else {
    attr.annot = NULL
  }
  
  if(do.col.clust) {
    ## cluster on normalized columns
    count.sub = count[,1:min(ncol(count),max.species.show)]
    # you have to transpose the dataset to get the taxa as rows
    data.dist.taxa <- vegdist(t(count.sub), method = dist.metr)
    col.clus <- hclust(data.dist.taxa, "ward.D2")
    col.dendro = as.dendrogram(col.clus)
  }
  else {
    col.dendro=NA
  }
  ## go back to un-normalized
  count.sub = count.src[,1:min(ncol(count),max.species.show)]
  if(!is.null(trans.show)) {
    count.sub = do.call(trans.show,list(count.sub))
  }
  if(!is.null(stand.show)) {
    count.sub = decostand(count.sub,method=stand.show,MARGIN=2)
  }
  
  pl.heat = annHeatmap2AT(count.sub,
                          col = colorRampPalette(c("lightyellow", "red"), space = "Lab")(50),
                          #col = heat.colors(50),
                          breaks = niceBreaks(c(min(count.sub),max(count.sub)),50),
                          scale = "none", # we get false bands with default standartization
                          legend = F,
                          dendrogram = list(status="yes",
                                            Row = list(dendro = row.dendro, status = ifelse(do.row.clust,"yes","no")), 
                                            Col = list(dendro = col.dendro, status = ifelse(do.col.clust,"yes","no"))),
                          labels = list(Col = list(nrow = 20),
                                        Row = if(is.null(attr.row.labels)) list(status="no") 
                                        else list(nrow=6,labels=attr[,attr.row.labels])), #cex=0.8
                          ann = list(Row = list(data = attr.annot,control=list(hbuff=0,vbuff=0))),
                          cluster = list(Row = list(cuth = cluster.row.cuth, 
                                                    col = function(n) {brewer.pal(n, "Set2")})) 
                          # cuth gives the height at which the dedrogram should be cut to form 
                          # clusters, and col specifies the colours for the clusters
  )
  report$add(plot(pl.heat),caption=caption)
}


ComplexHeatmap.add.attr <- function(attr.names,data,show_row_names = TRUE,width=NA,...) {
  n.h = length(attr.names)
  if(!is.null(width)) {
    if(is.na(width)) {
      width = grid::unit(1.0,"lines")
    }
  }
  get_pal = function(attr.name) {
    generate.colors.mgsat(data[,attr.name,drop=T],value="palette")
  }
  attr.name = attr.names[1]
  h = ComplexHeatmap::Heatmap(data[,attr.name,drop=F],name=attr.name,
                              show_row_names = ifelse(n.h>1,F,show_row_names),
                              width = width,
                              col = get_pal(attr.name),
                              ...)
  n.h = n.h - 1
  if (n.h > 0) {
    for(attr.name in attr.names[2:length(attr.names)]) {
      h = h + ComplexHeatmap::Heatmap(data[,attr.name,drop=F],name=attr.name,
                                      show_row_names = ifelse(n.h>1,F,show_row_names),
                                      width = width,
                                      col = get_pal(attr.name),
                                      ...) 
      n.h = n.h - 1
    }
  }
  return (h)
}

make.plot.ylim <- function(y,margin.rel=0.1,margin.abs=1e-16) {
  if(length(margin.abs)==1) {
    margin.abs = c(margin.abs,margin.abs)
  }
  ymin = min(y)
  ymax = max(y)
  ymin = ymin - margin.rel*(ymax-ymin)
  ymax = ymax + margin.rel*(ymax-ymin)
  if(ymax==ymin) {
    ymin = ymin - margin.abs[1]
    ymax = ymax + margin.abs[2]
  }
  c(ymin,ymax)
}

heatmap.feature.test.annot <- function(test.res,
                                       base=NULL,
                                       effect=NULL,
                                       p.val=NULL,
                                       p.val.adj=NULL,
                                       effect_baseline=0,
                                       p.val.adj.alpha=0.05,
                                       annot.which = "column",
                                       fontsize=NULL,
                                       fontsize_leg=NULL,
                                       cex_axis=0.6,
                                       cex_label=0.8) {
  library(ComplexHeatmap)
  gp_label = grid::gpar(fontsize=fontsize,cex=cex_label)
  annot_label = function(name,label=name) decorate_annotation(name, {grid.text(label, unit(1,"npc") + unit(0.5, "char"), 
                                                                               just = "left",
                                                                               gp = gp_label)})
  annot_labels = c()
  make_annot_label = function(name,label=name) {
    annot_labels <<- c(annot_labels,label)
    function() annot_label(name,label=label)
  }
  test.res = as.data.frame(test.res)
  n.feat = nrow(test.res)
  annots = list()
  decorations = list()
  if(!is.null(base)) {
    annots[[base]] = ComplexHeatmap::anno_points(test.res[[base]],
                                                 axis = T,
                                                 which=annot.which,
                                                 ylim = make.plot.ylim(test.res[[base]]),
                                                 axis_gp = gpar(fontsize=fontsize,cex=cex_axis))
    decorations[[base]] = make_annot_label(base)
  }
  if(!is.null(effect)) {
    y = test.res[[effect]]
    annots[[effect]] = ComplexHeatmap::anno_barplot(y,
                                                    baseline = effect_baseline,
                                                    axis = T,
                                                    which=annot.which,
                                                    ylim = make.plot.ylim(y),
                                                    gp = gpar(fill = ifelse(y > 0, "red", "blue")),
                                                    axis_gp = gpar(fontsize=fontsize,cex=cex_axis))
    decorations[[effect]] = make_annot_label(effect)
  }
  if(!is.null(p.val)) {
    y_orig = test.res[[p.val]]
    y = -log10(y_orig+1e-16)
    annots[[p.val]] = ComplexHeatmap::anno_points(y,
                                                  axis = T,
                                                  which=annot.which,
                                                  ylim = make.plot.ylim(y,margin.abs = c(0,1)),
                                                  gp = gpar(col = ifelse(y_orig <= p.val.adj.alpha, "green", "black")),
                                                  axis_gp = gpar(fontsize=fontsize,cex=cex_axis))
    decorations[[p.val]] = make_annot_label(p.val,label=sprintf("%s, -log10",p.val))
  }
  if(!is.null(p.val.adj)) {
    y = -log10(test.res[[p.val.adj]]+1e-16)
    annots[[p.val.adj]] = ComplexHeatmap::anno_points(y,
                                                      axis = T,
                                                      which=annot.which,
                                                      ylim = make.plot.ylim(y,margin.abs = c(0,1)),
                                                      axis_gp = gpar(fontsize=fontsize,cex=cex_axis))
    decorations[[p.val.adj]] = make_annot_label(p.val.adj,label=sprintf("%s, -log10",p.val.adj))
  }
  pad_annot_label = unit(0,"mm")
  for(a_l in annot_labels) {
    pad_annot_label = max(grid::grobWidth(grid::textGrob(a_l,gp = gp_label)),pad_annot_label)
  }
  ## because the default width for row cluster is 1cm, see
  ## https://bioconductor.org/packages/release/bioc/vignettes/ComplexHeatmap/inst/doc/s9.examples.html#toc_6
  padding = unit.c(unit(c(2, 2), "mm"),
                   unit(2, "mm"), max(pad_annot_label - unit(1, "cm"),unit(2, "mm")))
  list(
    plot=do.call(ComplexHeatmap::HeatmapAnnotation,
                 c(annots,
                   list(gap=grid::unit(0.5,"char")))),
    finalize=function() { for(decor in decorations) decor() },
    padding = padding
  )
}

heatmap.cluster.rows <- function(m_a,main.meta.var,clustering_distance_rows,km) {
  n.obs = nrow(m_a$count)
  ## "pearson", "spearman" and "kendall" are only understood by this internal function from ComplexHeatmap
  ## pam will silently use "euclidean"
  diss = ComplexHeatmap:::get_dist(m_a$count,clustering_distance_rows)
  
  if(n.obs >= 6) {
    if(km<1) {
      split = fpc::pamk(diss, krange = 1:min(n.obs-2,10))$pamobject$clustering
      split.descr = "Number of cluster splits is determined automatically with method `fpc::pamk`"
    }
    else {
      split = cluster::pam(diss, k=km)$clustering
      split.descr = "Number of cluster splits is set to a fixed value that is passed to method `cluster::pam`"
    }
  }
  else {
    split = NULL
    split.descr = "Not splitting clusters due to low number of observations"
  }
  if(!is.null(split)){
    if(num.levels(split)<2){
      split = NULL
      split.descr = "Splitting clusters resulted in a single partition"
    }
  }
  caption.g.test=sprintf("G-test of independence between automatic cluster splits and attribute '%s'. %s.",
                         main.meta.var,split.descr)
  g.t = NULL
  if(!is.null(split)) {
    if(num.levels(split)>1 && !is.null(main.meta.var) && num.levels(m_a$attr[,main.meta.var])>1) {
      g.t = g.test(m_a$attr[,main.meta.var],split)
    }
    m_a$attr$.Heatmap.Cluster.Split = split
  }
  list(split=split,m_a=m_a,g.t=g.t,main.meta.var=main.meta.var,caption.g.test=caption.g.test,split.descr=split.descr)  
}

heatmap.diff.abund <- function(m_a,
                               res.test.df=NULL,
                               base=NULL,
                               effect=NULL,
                               p.val=NULL,
                               p.val.adj=NULL,
                               effect_baseline=0,
                               p.val.adj.alpha=0.05,
                               hmap.label="Abundance",
                               hmap.width=1000,
                               hmap.height=hmap.width*0.8,
                               attr.annot.names=c(),
                               clustering_distance_rows="pearson",
                               cluster_columns=T,
                               km=0,
                               show_row_names = F,
                               show_column_names = T,
                               max.n.columns=NULL,
                               top_annotation_height = unit(10,"lines"),
                               column_names=NULL,
                               column_names_max_height = NULL,
                               column_names_gp = NULL,
                               max_column_names_symbols=30) {
  if(!is.null(res.test.df)) {
    res.test.df = as.data.table(res.test.df)
    res.test.df = res.test.df[data.table(ID=colnames(m_a$count)),on="ID"]
    stopifnot(all(res.test.df$ID==colnames(m_a$count)))
  }  
  if(!is.null(max.n.columns)) {
    max.n.columns = min(max.n.columns,ncol(m_a$count))
    m_a$count = m_a$count[,1:max.n.columns,drop=F]
  }
  main.meta.var = NULL
  if(length(attr.annot.names)>0) {
    main.meta.var = attr.annot.names[[1]]
  }
  fontsize = ggplot2::theme_get()$text$size
  fontsize_leg = fontsize*0.8
  #fontsize = grid::unit(fontsize,"points")
  #fontsize_leg = grid::unit(fontsize_leg,"points")
  library(ComplexHeatmap) # need it for `+`
  count = m_a$count
  colnames_count = colnames(count)
  if(!is.null(column_names)) {
    colnames_count = column_names
    if(!is.null(max.n.columns)) {
      colnames_count = colnames_count[1:max.n.columns]
    }
  }
  colnames_count = paste0(substr(colnames_count,1,max_column_names_symbols),
                          ifelse(stringr::str_length(colnames_count)>max_column_names_symbols,
                                 "...",
                                 ""))
  ## add a running index if any column names become identical after trimming
  if(anyDuplicated(colnames_count)) {
    colnames_count = paste0(seq_along(colnames_count),". ",colnames_count)
  }
  colnames(count) = colnames_count
  #column_names_max_height = grid::unit(1, "strwidth",colnames_count[which.max(stringr::str_length(colnames_count))[1]])
  if(is.null(column_names_gp)) {
    column_names_gp = grid::gpar(fontsize=fontsize,cex=0.8)
  }
  if(is.null(column_names_max_height)) {
    column_names_max_height = ComplexHeatmap::max_text_width(colnames_count, gp = column_names_gp)
  }
  #column_names_max_height = unit(10, "cm")
  rows.cluster = heatmap.cluster.rows(m_a=m_a,
                                      main.meta.var=main.meta.var,
                                      clustering_distance_rows=clustering_distance_rows,
                                      km=km)  
  top_annot = NULL
  if(!is.null(res.test.df)) {
    top_annot = heatmap.feature.test.annot(test.res=res.test.df,
                                           base=base,
                                           effect=effect,
                                           p.val=p.val,
                                           p.val.adj=p.val.adj,
                                           effect_baseline=effect_baseline,
                                           p.val.adj.alpha=p.val.adj.alpha,
                                           annot.which = "column",
                                           fontsize = fontsize)    
  }
  labels_gp = grid::gpar(fontsize = fontsize_leg)
  h = ComplexHeatmap::Heatmap(count,name=hmap.label,
                              cluster_columns=cluster_columns,
                              show_row_names = show_row_names,
                              show_column_names = show_column_names,
                              clustering_distance_rows = clustering_distance_rows, 
                              split=rows.cluster$split,
                              column_names_max_height = column_names_max_height,
                              column_names_gp = column_names_gp,
                              row_names_gp = grid::gpar(fontsize = fontsize),
                              heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
                                                          labels_gp = labels_gp,
                                                          grid_height = max(ComplexHeatmap::max_text_height("0",
                                                                                                            gp = labels_gp)+unit(0.1,"char"),
                                                                            unit(4, "mm"))),
                              top_annotation = top_annot$plot,
                              top_annotation_height = top_annotation_height)
  if(length(attr.annot.names)>0) {
    h = h + ComplexHeatmap.add.attr(attr.annot.names,
                                    m_a$attr,
                                    show_row_names=F,
                                    column_names_gp = column_names_gp,
                                    row_names_gp = grid::gpar(fontsize = fontsize),
                                    heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
                                                                labels_gp = labels_gp,
                                                                grid_height = max(ComplexHeatmap::max_text_height(attr.annot.names,
                                                                                                                  gp = labels_gp)+unit(0.1,"char"),
                                                                                  unit(4, "mm"))))
  }
  report$add({draw(h,padding=top_annot$padding); top_annot$finalize()},
             caption=sprintf("Clustered heatmap of normalized abundance values. %s.",rows.cluster$split.descr),
             width=hmap.width,height=hmap.height,hi.res.width = hmap.width, hi.res.height=hmap.height)
  if(!is.null(rows.cluster$split)) {
    report$add(rows.cluster$g.t,caption = rows.cluster$caption.g.test)
  }
  export.taxa.meta(rows.cluster$m_a,
                   label="htmap",
                   descr="Data used for heatmap with added row cluster splits",
                   row.proportions=F,
                   row.names=F)
}

heatmap.combined.report <- function(m_a,
                                    m_a.norm,
                                    res.tests,
                                    norm.count.task=NULL,
                                    hmap.width,
                                    hmap.height=hmap.width*0.8,
                                    attr.annot.names=c(),
                                    clustering_distance_rows="pearson",
                                    cluster_columns=T,
                                    km.abund=0,
                                    km.diversity=0,
                                    show_row_names = F,
                                    show_column_names = T,
                                    max.n.columns=NULL) {
  if(!is.null(max.n.columns)) {
    max.n.columns = min(max.n.columns,ncol(m_a.norm$count))
    m_a.norm$count = m_a.norm$count[,1:max.n.columns,drop=F]
  }
  main.meta.var = NULL
  if(length(attr.annot.names)>0) {
    main.meta.var = attr.annot.names[[1]]
  }
  fontsize = ggplot2::theme_get()$text$size
  fontsize_leg = fontsize*0.8
  #fontsize = grid::unit(fontsize,"points")
  #fontsize_leg = grid::unit(fontsize_leg,"points")
  library(ComplexHeatmap) # need it for `+`
  count = m_a.norm$count
  colnames_count = colnames(count)
  #column_names_max_height = grid::unit(1, "strwidth",colnames_count[which.max(stringr::str_length(colnames_count))[1]])
  column_names_max_height = ComplexHeatmap::max_text_width(colnames_count, gp = gpar(fontsize = fontsize))
  rows.cluster = heatmap.cluster.rows(m_a=m_a.norm,
                                      main.meta.var=main.meta.var,
                                      clustering_distance_rows=clustering_distance_rows,
                                      km=km.abund)  
  h = ComplexHeatmap::Heatmap(count,name="Abundance",
                              cluster_columns=cluster_columns,
                              show_row_names = show_row_names,
                              show_column_names = show_column_names,
                              clustering_distance_rows = clustering_distance_rows, 
                              split=rows.cluster$split,
                              column_names_max_height = column_names_max_height,
                              column_names_gp = grid::gpar(fontsize = fontsize,cex=0.8),
                              row_names_gp = grid::gpar(fontsize = fontsize),
                              heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
                                                          labels_gp = grid::gpar(fontsize = fontsize_leg))) +
    #     ComplexHeatmap::HeatmapAnnotation(df=m_a.norm$attr[,attr.annot.names,drop=F],
    #                             which="row",
    #                             annotation_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
    #                                                         labels_gp = grid::gpar(fontsize = fontsize)),
    #                             show_annotation_name = T,
    #                             annotation_name_gp = grid::gpar(fontsize = fontsize))
    
    ComplexHeatmap.add.attr(attr.annot.names,
                            m_a.norm$attr,
                            show_row_names=F,
                            column_names_gp = grid::gpar(fontsize = fontsize,cex=0.8),
                            row_names_gp = grid::gpar(fontsize = fontsize),
                            heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
                                                        labels_gp = grid::gpar(fontsize = fontsize_leg)))
  
  report$add(h,caption=sprintf("Clustered heatmap of normalized abundance values. %s.",rows.cluster$split.descr),
             width=hmap.width,height=hmap.height,hi.res.width = hmap.width, hi.res.height=hmap.height)
  if(!is.null(rows.cluster$split)) {
    report$add(rows.cluster$g.t,caption = rows.cluster$caption.g.test)
  }
  mor_rowAnnotations <- m_a.norm$attr[,attr.annot.names,drop=F]
  mor_p = morpheus::morpheus(m_a.norm$count,
           rowAnnotations=mor_rowAnnotations,
           rows=c(list(list(field="id", display=list("text"))),
                  lapply(attr.annot.names,function(x) list(field=x,display=list("text","color")))
                  )
           )
  report$add.package.citation("morpheus")
  report$add.widget(mor_p,show.inline=F,
                    caption = "Dynamic Morpheus heatmap of normalized abundance values. 
                    It is available here through the link only because it can take a while to render for large datasets.
                    This is very customizable. What you will see initially is just a default starting configuration. Explore its menus.")
  
  export.taxa.meta(rows.cluster$m_a,
                   label="htmap",
                   descr="Data used for heatmap with added row cluster splits (clustering by abundance profile)",
                   row.proportions=F,
                   row.names=F)
  
  if(!is.null(get.diversity(res.tests,type="diversity"))) {
    div = log(get.diversity(res.tests,type="diversity")$e)
    m_a.norm$count = div
    clustering_distance_rows = "pearson"
    rows.cluster = heatmap.cluster.rows(m_a=m_a.norm,
                                        main.meta.var=main.meta.var,
                                        clustering_distance_rows=clustering_distance_rows,
                                        km=km.diversity)  
    h.d = ComplexHeatmap::Heatmap(div,name="Renyi diversity indices",
                                  cluster_columns=F,
                                  show_row_names = F, 
                                  clustering_distance_rows = clustering_distance_rows, 
                                  split=rows.cluster$split,
                                  column_names_gp = grid::gpar(fontsize = fontsize,cex=0.8),
                                  width=grid::unit(ncol(m_a.norm$count),"lines"),
                                  heatmap_legend_param = list(title_gp = grid::gpar(fontsize = fontsize), 
                                                              labels_gp = grid::gpar(fontsize = fontsize_leg)))
    h = h.d + h
    report$add(h,caption=sprintf("Clustered heatmap of diversity and normalized abundance values. %s.",rows.cluster$split.descr),
               width=hmap.width,height=hmap.height,hi.res.width = hmap.width, hi.res.height=hmap.height)
    if(!is.null(rows.cluster$split)) {
      report$add(rows.cluster$g.t,caption = rows.cluster$caption.g.test)
    }
    export.taxa.meta(rows.cluster$m_a,
                     label="htmap",
                     descr="Data used for heatmap with added row cluster splits (clustering by Renyi diversity indices)",
                     row.proportions=F,
                     row.names=F)
  }
}

#' turn ordered quantiiles (as returned by quantcut.ordered()) into a color gradient
#' x - original data vector
#' x_q - quantiles vector same length as x
ordered.color <- function(x,x_q,as.rgb=F) {
  names(x_q) = names(x)
  pal_tbl = data.table(x_q=levels(x_q))[,x_q:=factor(x_q,levels=x_q,ordered = T)]
  pal_tbl[,col:=scales::seq_gradient_pal("blue", "red", "Lab")(vegan::decostand(.I,method = "range"))]
  pal_tbl[,col:=factor(col,levels = col, ordered = T)]
  if(as.rgb) {
    pal_tbl = cbind(pal_tbl,t(col2rgb(pal_tbl$col)))
  }
  y = data.table(x=x,x_q=x_q)
  y = pal_tbl[y,on="x_q"]
  if(!is.null(names(x_q))) {
    y[,rn:=names(x_q)]
  }
  list(data=y,pal_tbl=pal_tbl)
}

quantcut.ordered.color <- function(x,q=4,as.rgb=F,...) {
  x_q = quantcut.ordered(x,q=q,...)
  ordered.color(x,x_q,as.rgb = as.rgb)
}

#' qcol as returned by quantcut.ordered.color
ordered.color.legend <- function(qcol,title) {
  library(data.table)
  library(ggplot2)
  tbl_col = qcol$data
  tbl_col = tbl_col[,.(count=.N),by=col]
  tbl_col = tbl_col[qcol$pal_tbl,on="col"][order(col)]
  tbl_col[,right := as.numeric(stringr::str_match(x_q,",([.0-9eE]+)")[,2])]
  tbl_col[,left := as.numeric(stringr::str_match(x_q,"([.0-9eE]+),")[,2])]
  ticks = c(tbl_col[1,left],tbl_col[,right])
  ticks = sprintf("%s%%",ticks*100)
  ticks = c("",ticks[2:(length(ticks)-1)],"")
  col_pal = as.character(tbl_col$col)
  names(col_pal) = tbl_col$x_q
  #ticks = rev(ticks)
  tbl_col[,x_q:=ordered(x_q,levels=rev(levels(x_q)))]
  theme_font_size_abs = ggplot2::theme_get()$text$size
  fontsize = 3
  fontsize_final = fontsize*theme_font_size_abs
  pl = ggplot(tbl_col,aes(x=factor(1),fill=x_q)) + 
    geom_bar(position="stack",color="black") + 
    scale_fill_manual(values=col_pal) +
    scale_y_continuous(name = title, 
                       breaks=(seq_along(ticks)-1),
                       labels = ticks) +
    ggtitle(title) +
    #coord_flip() +
    theme(text=element_text(color=c("black","black"),size = fontsize_final),
          plot.title = element_text(size=fontsize_final),
          axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          legend.position="none",
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank())
  return(list(tbl_palette=tbl_col,plot_palette=pl))
}

make.color.legend.scatter.js3d <- function(xyz,color_val,color) {
  xyz_max = apply(xyz,2,min)
}

plot.scatter.js3d <- function(xyz,data,color=NULL,
                              labels=NULL,size=NULL,pch=NULL,
                              renderer="auto",show.color.legend=T,
                              num.ticks=c(6,6,6),
                              axis.scale=NA,
                              lines.args=NULL,
                              cex.lab=0.8,
                              ...) {
  library(threejs)
  args = list(color=color,labels=labels,size=size,pch=pch)
  args = interpret.args.in.df(args,data)
  
  color_val = args$color
  if(!is.null(args$color) && length(args$color) > 1) {
    args$color = generate.colors.mgsat(args$color)
  }
  else {
    show.color.legend = F
  }
  if(!is.null(args$size)) {
    args$size = vegan::decostand(args$size,method="range")*2+1
  }
  pl = do.call(threejs::scatterplot3js,
               c(list(as.matrix(xyz)),
                 args,
                 list(renderer=renderer,
                      num.ticks=num.ticks,
                      axis.scale=axis.scale,
                      cex.lab=cex.lab),
                 list(...)
               )
  )
  if(F && show.color.legend) {
    points3d_args = make.color.legend.scatter.js3d(xyz,color_val,color)
    pl = do.call(points3d,c(list(pl),
                            points3d_args))
  }
  if(!is.null(lines.args)) {
    pl = do.call(lines3d,c(list(pl),
                           lines.args))
  }
  #scatterplot3js(x, y, z, pch="@") %>%
  #points3d(pl,x + 0.1, y + 0.1, z, color="red", pch=paste("point", 1:5))
  pl
}

## Wrapper around phyloseq::ordinate
phyloseq.ordinate <- function(physeq, method = "DCA", distance = "bray", formula = NULL, ...) {
  args = list(...)
  
  library(phyloseq)
  ## Call metaMDS directly when distance matrix is provided because phyloseq::ordinate
  ## fails to pass k parameter to metaMDS
  
  ord.func = "ordinate"
  ord.args = c(list(physeq,method=method,distance=distance),args)
  
  if (method == "NMDS") {
    if (inherits(distance, "dist")) {
      ord.func = vegan::metaMDS
      ord.args = c(list(distance),args)
    }
  }
  
  ord = do.call(ord.func,
                ord.args
  )
  return (ord)  
}

plot_ordination.2d <- function(physeq,ordination,
                               type = "samples", axes = 1:2,
                               legend.point.size = rel(4),
                               legend.position="right",
                               border=T,
                               ggplot.extra=list(),
                               ...) {
  library(phyloseq)
  pt = list(...)
  
  df.plot =  plot_ordination(physeq,ordination,type=type,axes=axes,justDF = T)
  
  if(is.null(pt$alpha) && nrow(df.plot)>100) {
    ##helps with overplotting
    pt$alpha = 0.5
  }
  
  df.names = colnames(df.plot)
  pl = ggplot(df.plot,aes_string(x=df.names[1],y=df.names[2]))
  
  if(!is.null(pt$color)) {
    args = list(color=pt$color)
    args = interpret.args.in.df(args,df.plot)
    
    if(!is.null(args$color) && length(args$color) > 1) {
      palette = generate.colors.mgsat(args$color,value="palette")
      guide_ovr = guide_legend(override.aes = list(size=legend.point.size))
      pl = pl + 
        scale_size(range = c(4,10)) +
        scale_fill_manual(values = palette) +
        #scale_color_manual(values = palette) +
        #guides(colour = guide_ovr,fill = guide_ovr)
        guides(fill = guide_ovr)
      ##TODO: the above size override will probably break color bar because
      ##it always forces guide_legend, which is an alternative to color bar.
      ##Need to figure out when color is continuous scale and use alternative
      ##override. As of 2017/07, there seems to be no way to increase legend
      ##point size in a generic way.
    }
  }
  
  options.fixed = list()
  if(border) {
    pt$fill = pt$color
    pt$color = NULL
    options.fixed$shape=21
  }
  if(is.null(pt$size)) {
    options.fixed$size=rel(4)
  }
  
  pl = pl + make.geom(geom_point,data=df.plot,
                      options=pt,options.fixed = options.fixed)
  
  pl = pl + theme_bw(base_size = 18) + theme(legend.position=legend.position)
  if(!is.null(ggplot.extra)) {
    for(g in ggplot.extra) {
      pl = pl + g
    }
  }
  return (pl)
}

interpret.args.point_lines <- function(dt,
                                       axes = 1:2,
                                       line.group=NULL,
                                       line.order=NULL,
                                       line.order.cyclic=F,
                                       line.color.order=T,
                                       find.hull=F,
                                       line.width=5,
                                       line.color=NULL,
                                       line.alpha=NULL) {
  ret = NULL
  dt = copy(data.table::as.data.table(dt))
  dt[,.ind_point:=.I]
  if (!is.null(line.color)) {
    dt[,.line_color:=line.color]
  }
  dt[[".line_order"]] = NULL
  if(!is.null(line.order)) {
    setorderv(dt,c(line.group,line.order))
    dt[[".line_order"]] = dt[[line.order]]
  }
  if(line.color.order) {
    dt[,.line_color := quantcut.ordered.color(as.numeric(.line_order))$data$col]
  }
  
  if(!is.null(line.group)) {
    dt[, `:=`(.ind_point_to = 
                c(if(.N>1) .ind_point[2:.N] else as.integer(c()),
                  if(line.order.cyclic) .ind_point[1] else NA),
              .ind_point_group = seq(.N)), by = line.group]
  }
  dt = dt[!is.na(.ind_point_to)]
  ind_point_from = dt[,.ind_point]
  ind_point_to = dt[,.ind_point_to]
  line.color = dt[,.line_color]
  if(length(axes)==3) {
    ret = list(from=ind_point_from,
               to=ind_point_to,
               lwd=line.width,
               alpha=line.alpha,
               color=line.color)
    ret = plyr::compact(ret)
  }
  else {
    stop("Not implemented yet for dimensions other than 3")
  }
  ret
}

plot_ordination.3d <- function(physeq,ordination,
                               type = "samples", axes = 1:3,
                               lines.args = NULL,
                               ...) {
  library(phyloseq)
  pt = list(...)
  
  df.plot =  plot_ordination(physeq,ordination,type=type,axes=axes,justDF = T)
  if(!is.null(lines.args)) {
    lines.args = do.call(interpret.args.point_lines,
                         c(list(df.plot,axes=axes),
                           lines.args))
  }
  
  pl = do.call(plot.scatter.js3d,
               c(list(df.plot[,seq_along(axes),drop=F],df.plot,
                      lines.args=lines.args),
                 pt))
  
  
  return (pl)
}


ordination.report <- function(m_a,res=NULL,distance="bray",ord.tasks,sub.report=T,descr="",
                              legend.position="bottom") {
  require(phyloseq)
  report.section = report$add.header(paste("Ordinations",descr,sep = ", "),section.action="push",sub=sub.report)  
  report$add.package.citation("phyloseq")  
  report$add.package.citation("vegan")  
  ph = m_a.to.phyloseq(m_a)
  for(ord.task in ord.tasks) {
    if(!is.null(ord.task$ordinate.task$formula)) {
      ord.task$ordinate.task$formula = as.formula(sprintf("~%s",ord.task$ordinate.task$formula))
    }
    
    ## Call metaMDS directly when distance matrix is provided because phyloseq::ordinate
    ## fails to pass k parameter to metaMDS
    
    ord.func = "phyloseq.ordinate"
    ord.args = c(list(ph,distance=distance),ord.task$ordinate.task)
    
    
    ord = do.call(ord.func,
                  ord.args
    )
    
    pt.orig = ord.task$plot.task
    if(is.null(pt.orig$axes)) {
      pt.orig$axes = 1:2
    }
    if(is.null(pt.orig$type)) {
      pt.orig$type = "samples"
    }
    
    pt = pt.orig
    pt$axes = NULL
    pt$type = NULL
    pt$pch = NULL
    pt$lines.args = NULL
    #pt.ggplot.extra = pt$ggplot.extra
    #pt.legend.point.size = pt$legend.point.size
    #if(is.null(pt.legend.point.size)) {
    #  pt.legend.point.size = rel(4)
    #}
    #pt$ggplot.extra = NULL
    #pt$legend.point.size= NULL
    
    pt = plyr::compact(pt)
    pl =  do.call(plot_ordination.2d,
                  c(list(ph,ord,type=pt.orig$type,
                         axes=pt.orig$axes,
                         #legend.point.size=pt.legend.point.size,
                         legend.position=legend.position
                         #ggplot.extra=pt.ggplot.extra
                  ),
                  pt
                  ))
    
    caption=sprintf("Ordination plot. Ordination performed with parameters %s. 
               Plot used parameters %s.",
                    arg.list.as.str(ord.task$ordinate.task),
                    arg.list.as.str(pt))
    report$add(pl,caption = caption)
    
    if(inherits(distance,"dist") || nrow(m_a$count) > 3) {
      
      ##We have to redo NMDS for 3D plot if the original number of requested ordination dimensions
      ##was less than 3
      if(ord.task$ordinate.task$method == "NMDS") {
        if(is.null(ord.task$ordinate.task$k) || ord.task$ordinate.task$k<3) {
          ord.args$k = 3
          ord = do.call(ord.func,
                        ord.args
          )        
          
        }
      }
      
      pt = pt.orig
      pt$axes = 1:3
      
      caption=sprintf("Ordination plot in 3D. Ordination performed with parameters %s. 
               Plot used parameters %s.",
                      arg.list.as.str(ord.task$ordinate.task),
                      arg.list.as.str(pt))
      report$add.widget(plot_ordination.3d(
        ph,ord,type=pt$type,axes=pt$axes,labels=pt$label,color=pt$color,size=pt$size,
        pch=pt$pch,
        lines.args = pt$lines.args,
        axis.scale = if(!is.null(pt$axis.scale)) pt$axis.scale else NA),
        caption = caption)
    }
  }
  report$pop.section()  
  #ordination.report(m_a,res=NULL,distance=as.dist(d.direct.mat),ord.tasks=list(list(ordinate.task=list(method="NMDS"),plot.task=list(type="samples",color="Genotype",size=5)),sub.report=F))
}

make.geom <- function(geom,data,options,options.data=list(),options.fixed=list(),include.data=F) {
  data.names = colnames(data)
  options = options[!sapply(options,is.null)]
  options.data.mask = options %in% data.names
  options.data = c(options[options.data.mask],options.data)
  options.fixed = c(options[!options.data.mask],options.fixed)
  if(include.data) {
    options.fixed$data = data
  }
  #make.global(name="gm")
  do.call(geom,
          c(
            ifelse(length(options.data)>0,list(mapping=do.call(aes_string,options.data)),list(mapping=NULL)),
            options.fixed
          )
  )
}

igraph.to.d3net <- function(g,vertex.data=NULL) {
  vertex.name = as.character(as_ids(V(g)))
  vertex.ind = seq_along(vertex.name) - 1 #JS zero offset
  names(vertex.ind) = vertex.name
  edge = as_data_frame(g)
  edge$source = vertex.ind[edge$from]
  edge$target = vertex.ind[edge$to]
  if(!is.null(vertex.data)) {
    vertex.data = vertex.data[vertex.name, , drop = FALSE]
    if(!nrow(vertex.data)) {
      stop("Some graph vertices are missing from vertex data table")
    }
  }
  return (list(egde=edge,vertex=data.frame(name=vertex.name),vertex.data=vertex.data))
}

mgsat.plot.igraph.vertex.options = list(size=4,alpha=0.75)
mgsat.plot.igraph.vertex.text.options = list(label = "vertex.name", hjust = 0.5, vjust = 1.5, size = 4)
mgsat.plot.igraph.edge.options = list(size = 0.5, color = "black", alpha = 0.4)

## Plot igraph object with ggplot2.
## Based on code from phyloseq:::plot_network
## *.options must be named lists, where any specified element will override default value set in
## corresponding mgsat.plot.igraph.*.options package veriable. Provide NULL element value to remove
## default value.
## vertex.data, if set, must be a data frame with row names matching values of vertices 'name' attribute,
## which can be referenced as 'vertex.name' in aesthetic parameters.
## *.options elements must match arguments to geom_point for vertex, geom_text for vertex.text 
## and geom_line for edge. Elements will be checked against vertex.data, and if matching column names,
## will be mapped within aes_string() function, otherwise provided as fixed values outside of aes in
## geom constructor.
## Value: ggplot object.
mgsat.plot.igraph <- function (g, vertex.data = NULL, 
                               vertex.options = mgsat.plot.igraph.vertex.options,
                               vertex.text.options = mgsat.plot.igraph.vertex.text.options, 
                               edge.options = mgsat.plot.igraph.edge.options,
                               vertex.text.selection=NULL,
                               layout = "layout.fruchterman.reingold",
                               extra.plot.operands=list()) 
{
  library(igraph)
  vertex.options = update.list(mgsat.plot.igraph.vertex.options,vertex.options)
  vertex.text.options = update.list(mgsat.plot.igraph.vertex.text.options,vertex.text.options)
  edge.options = update.list(mgsat.plot.igraph.edge.options,edge.options)
  if (vcount(g) < 2) {
    stop("The graph you provided, `g`, has too few vertices.")
  }
  if(is.function(layout) || is.character(layout)) {
    vertDF <- do.call(layout,list(g))
  }
  else {
    vertDF = layout
  }
  colnames(vertDF) <- c("x", "y")
  vertDF <- data.frame(vertex.name = get.vertex.attribute(g, "name"), 
                       vertDF)
  if (!is.null(vertex.data)) {
    vertDF <- data.frame(vertDF, vertex.data[as.character(vertDF$vertex.name), 
                                             , drop = FALSE])
  }
  p <- ggplot(vertDF, aes(x, y))
  p <- p + theme_bw() + theme(panel.grid.major = element_blank(), 
                              panel.grid.minor = element_blank(), axis.text.x = element_blank(), 
                              axis.text.y = element_blank(), axis.title.x = element_blank(), 
                              axis.title.y = element_blank(), axis.ticks = element_blank(), 
                              panel.border = element_blank())
  p <- p + make.geom(geom_point,data=vertDF,
                     options=vertex.options,
                     options.fixed=list(na.rm=T))
  if (!is.null(vertex.text.options$label)) {
    text.data = vertDF
    text.include.data = F
    if(!is.null(vertex.text.selection)) {
      text.data = text.data[text.data$vertex.name %in% vertex.text.selection,,drop=F]
      text.include.data = T
    }
    if(nrow(text.data)>0) {
      p <- p + make.geom(geom_text,data=text.data,
                         options=vertex.text.options,
                         options.fixed=list(na.rm = TRUE),
                         include.data=text.include.data)
    }
  }
  edgeDF <- data.frame(get.edgelist(g))
  if(nrow(edgeDF) > 0) {
    edgeDF$id <- 1:length(edgeDF[, 1])
    graphDF <- merge(melt(edgeDF, id = "id", value.name = "vertex.name"), vertDF, by = "vertex.name")
    p <- p + make.geom(geom_line,data=graphDF,
                       options=edge.options,
                       options.data=list(group="id"),
                       options.fixed=list(na.rm=T),
                       include.data=T)
  }
  ## give bigger margins to avoid cutting off the point labels
  p = p + scale_x_continuous(expand=c(0.1,0))
  for(extra.plot.operand in extra.plot.operands) {
    p = p + extra.plot.operand
  }
  return(p)
}

scale.to.range <- function(x,quant.min=0.1,quant.max=0.9,na.rm=T) {
  require(vegan)
  q = quantile(x,c(quant.min,quant.max),na.rm=na.rm)
  x[x<q[1]] = q[1]
  x[x>q[2]] = q[2]
  return (decostand(x,method="range"))
}

mgsat.plot.igraph.d3net <- function (g, vertex.data = NULL, 
                                     vertex.options = mgsat.plot.igraph.vertex.options,
                                     vertex.text.options = mgsat.plot.igraph.vertex.text.options, 
                                     edge.options = mgsat.plot.igraph.edge.options,
                                     ...) 
{
  library(networkD3)
  vertex.options = update.list(mgsat.plot.igraph.vertex.options,vertex.options)
  vertex.text.options = update.list(mgsat.plot.igraph.vertex.text.options,vertex.text.options)
  edge.options = update.list(mgsat.plot.igraph.edge.options,edge.options)
  if (vcount(g) < 2) {
    return (list(plot.obj=NULL,msg="Graph must have more than one vertix."))
  }
  
  x = igraph.to.d3net(g,vertex.data = vertex.data)
  edge = x$egde
  if(nrow(edge) == 0) {
    return (list(plot.obj=NULL,msg="Graph must have at least one edge."))
  } 
  vertex.data = x$vertex.data
  if(is.null(vertex.data)) {
    vertex.data = x$vertex
  }
  vertex.data$vertex.name = x$vertex$name
  
  radiusCalculation=NULL
  Nodesize=NULL
  if(is.numeric(vertex.options$size)) {
    radiusCalculation=JS(vertex.options$size)
  }
  else {
    Nodesize = ".Nodesize"
    vertex.data$.Nodesize = scale.to.range(vertex.data[,vertex.options$size])
    radiusCalculation=JS(sprintf('Math.sqrt(d.nodesize)*6+3'))
  }
  Group = vertex.options$color
  colorScale = NULL
  if(!is.null(Group)) {
    if(is.numeric(vertex.data[,Group])) {
      vertex.data[,Group] = signif(vertex.data[,Group],6)
      q = quantile(vertex.data[,Group],c(0.1,0.5,0.9))
      colorScale=JS(sprintf('d3.scaleLinear()
                     .domain([%s, %s, %s])
                     .range(["blue", "grey", "red"])',q[1],q[2],q[3]))
    }
    else {
      colorScale=JS('d3.schemeCategory10')
    }
  }
  NodeID = vertex.text.options$label
  opacity = vertex.options$alpha
  p = forceNetwork(Links = edge, 
                   Nodes = vertex.data,
                   Source = "source", 
                   Target = "target",
                   NodeID = NodeID,
                   Group = Group, 
                   Nodesize =  Nodesize, 
                   opacity = opacity,
                   legend=T,
                   zoom = T, 
                   clickAction = 'd.fixed = !d.fixed',
                   radiusCalculation=radiusCalculation,
                   colourScale = colorScale,
                   fontSize = 11)
  return (list(plot.obj=p,msg="OK"))
}


network.spiec.easi.options = list(
  method='mb', 
  lambda.min.ratio=1e-2, 
  nlambda=80, #15
  icov.select.params=list(rep.num=50)
)

network.spiec.easi <- function(count,
                               ...) {
  report$add.package.citation("SpiecEasi")
  options = update.list(network.spiec.easi.options,list(...))
  #se.est = dbg.cache$se.est
  #if(is.null(se.est)) {
  ## without bringing everything to the global namespace, the package gives an error of
  ## not being able to find some function through a get() call: 
  library(SpiecEasi)
  se.est <- do.call(SpiecEasi::spiec.easi,c(list(count),options))
  #}
  rownames(se.est$refit$stars) = colnames(count)
  colnames(se.est$refit$stars) = colnames(count)
  adj = getRefit(se.est)
  gr = adj2igraph(adj,vertex.attr = list(name = rownames(adj)))
  return (list(method.res=se.est,graph=gr))
}  

## vertex.data can be also provided inside all or some elements of plot.tasks,
## in which case it takes precedence over the value provided at the function call level.
network.report <- function(m_a,
                           count.filter.options=NULL,
                           vertex.data=NULL,
                           plot.tasks,
                           sub.report=T,
                           method="network.spiec.easi",
                           method.options=list(),
                           descr="") {
  report.section = report$add.header(sprintf("Network Analysis %s",descr),section.action="push",sub=sub.report) 
  report$add.descr(sprintf("Build network of interactions between features or samples and plot it. 
                           Network method is %s. Method parameters are: %s.",
                           method,
                           arg.list.as.str(method.options)))
  m_a = report.count.filter.m_a(m_a=m_a,
                                count.filter.options=count.filter.options,
                                descr="Network Analysis",warn.richness=F)
  
  net.res = do.call(method,c(list(m_a$count),method.options))
  
  gr = net.res$gr
  layout = igraph::layout.fruchterman.reingold(gr)
  for(plot.task in plot.tasks) {
    caption = plot.task$descr
    caption = sprintf("Network analysis with method %s. %s",method,caption)
    plot.task$descr = NULL
    if(is.null(plot.task$vertex.data)) {
      plot.task$vertex.data = vertex.data
    }
    gp = do.call(mgsat.plot.igraph,
                 c(list(gr,layout=layout),
                   plot.task
                 )
    )
    report$add(gp,caption=caption)
    gp.3d.res = do.call(mgsat.plot.igraph.d3net,
                        c(list(gr),
                          plot.task
                        )
    )
    if(!is.null(gp.3d.res$plot.obj)) {
      report$add.widget(gp.3d.res$plot.obj,caption=caption)
    }
    else {
      report$add.p(sprintf("Not creating 3D network plot. %s",gp.3d.res$msg))
    }
    
  }
  report$pop.section()  
}


network.features.combined.report <- function(m_a,
                                             res.tests=NULL,
                                             count.filter.options=NULL,
                                             sub.report=T,
                                             drop.unclassified=T,
                                             max.vertex.labels=30,
                                             method="network.spiec.easi",
                                             method.options=list(),
                                             descr="") {
  descr = sprintf("Feature correlation with overlaid differential abundance results %s",descr)
  if(drop.unclassified) {
    drop.names = count.filter.options$drop.names
    drop.names = c(drop.names,"other",colnames(m_a$count)[grepl("Unclassified.*",colnames(m_a$count),ignore.case=T)])
    if(is.null(count.filter.options)) {
      count.filter.options = list()
    }
    count.filter.options$drop.names = drop.names
  }
  ds2.res = get.feature.ranking(res.tests,method="deseq2",only.names=F,id.result=0)
  plot.tasks = foreach(res = ds2.res$ranked$results) %do% {
    res.descr = get.deseq2.result.description(res)
    res.df = as.data.frame(res)
    res.df$logBaseMean = log(res.df$baseMean+0.5)
    res.df$statistics = res.df$stat
    vertex.options = list(size="logBaseMean",color = "statistics", alpha=1)
    vertex.text.options = list(size = rel(3))
    vertex.text.selection = rownames(res.df)[1:min(max.vertex.labels,nrow(res.df))]
    list(vertex.data=res.df,
         vertex.options=vertex.options,
         vertex.text.options=vertex.text.options,
         vertex.text.selection=vertex.text.selection,
         descr=sprintf("Vertices are labeled by DESeq2 results for %s. Showing names for the maximum of %s top-ranked features.",res.descr,max.vertex.labels),
         extra.plot.operands=list(
           scale_colour_gradient2(midpoint=0,mid="grey",high="red",low="blue"),
           scale_size(range = c(1, 8))
         )
    )
  }
  
  if(length(plot.tasks)>0) {
    
    network.report(m_a,
                   count.filter.options=count.filter.options,
                   vertex.data=NULL,
                   plot.tasks=plot.tasks,
                   sub.report=sub.report,
                   method=method,
                   method.options=method.options,
                   descr=descr)
  }
  
}

plot.features.mds <- function(m_a,species.sel=NULL,sample.group=NULL,show.samples=T,show.species=T) {
  ##https://stat.ethz.ch/pipermail/r-help/2008-April/159351.html
  ##http://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf
  
  m = m_a$count
  ##default distance is bray
  m = decostand(m,method="range",MARGIN=2)
  sol = metaMDS(m,autotransform = F,trymax=40)
  if(show.samples) {
    site.sc <- scores(sol, display = "sites")
    plot(site.sc)
    sample.group = factor(sample.group)
    unique_sample.group = levels(sample.group)
    points(sol,display="sites",col=sample.group)
    if(!is.null(unique_sample.group)) {
      legend(-0.5,1,unique_sample.group,col=1:length(sample.group),pch=1)  
    }
  }
  if(show.species) {
    species.sc <- scores(sol, display = "species")
    if(!is.null(species.sel)) {
      species.sel.mask=(colnames(m_a$count) %in% species.sel)
    }
    else {
      species.sel.mask=rep(T,length(colnames(m_a$count)))
    }
    points(sol,display="species",pch="x",select=species.sel.mask)
    text(sol,display="species", cex=0.7, col="blue",select=species.sel.mask)
  }
  #points(site.sc,col=m_a$attr$T1D)
  #points(species.sc)
  #points(species.sc["Streptococcus_0.1.11.1.2.6.2",1:2,drop=F],pch="+")
  sol
}

cut.top.predictions<-function(scores,labels,sample.ids,n.cut,filter.by.label=F) {
  lab.levels = levels(labels)
  stopifnot(length(lab.levels)==2)
  n = length(scores)
  stopifnot(length(labels)==n && length(sample.ids)==n)
  
  n.cut = min(n.cut,round(n/2))
  
  res = list()
  
  if (filter.by.label) {
    
    lab = lab.levels[1]
    ord = order(scores)
    cum.cnt = cumsum(labels[ord]==lab)
    i.ord = first(cum.cnt>=n.cut)
    stopifnot(length(i.ord)>0) #index found
    cut = scores[ord[i.ord]]
    ids = sample.ids[(scores <= cut) & (labels == lab)]
    res[[lab]] = list(cut=cut,ids=ids)
    
    lab = lab.levels[2]
    ord = order(-scores)
    cum.cnt = cumsum(labels[ord]==lab)
    i.ord = first(cum.cnt>=n.cut)
    stopifnot(length(i.ord)>0) #index found
    cut = scores[ord[i.ord]]
    ids = sample.ids[(scores >= cut) & (labels == lab)]
    res[[lab]] = list(cut=cut,ids=ids)
    
  }
  else {
    ord = order(scores)
    cut = scores[ord[n.cut]]
    ids = sample.ids[ord[1:n.cut]]
    res[[lab.levels[1]]] = list(cut=cut,ids=ids)
    
    ord = order(-scores)
    cut = scores[ord[n.cut]]
    ids = sample.ids[ord[1:n.cut]]
    res[[lab.levels[2]]] = list(cut=cut,ids=ids)
  }
  return (res)
  
}

show.pred.perf <- function(pred,pred.score,labels) {
  library(ROCR)
  report.section = report$add.header("Prediction performance measures",section.action="push")  
  report$add.package.citation("ROCR")
  #print(table(labels,pred.score>0))
  pred.perf = prediction(pred.score,labels)
  report$add.table(table(labels,pred),caption="Confusion table")
  # Plot ROC curve
  perf <- performance(pred.perf, measure = "tpr", x.measure = "fpr")
  report$add(plot(perf),caption="ROC curve")
  # Plot precision/recall curve
  perf <- performance(pred.perf, measure = "prec", x.measure = "rec")
  report$add(plot(perf),caption="Precision/recall")
  perf <- performance(pred.perf, measure = "acc")
  report$add(plot(perf),caption="Accuracy")
}

select.samples <- function(m_a,
                           species.sel,
                           sample.group.name,
                           n.select=20,
                           filter.by.label=T,
                           selection.file=NULL) {
  library(kernlab)
  library(caret)
  
  report.section = report$add.header("Selecting most different samples with regard to a phenotype",section.action="push")
  report$add.package.citation("kernlab")
  report$add.package.citation("caret")
  report$add.descr(paste("This procedure selects those samples that are most different with regard to a grouping variable"
                         ,sample.group.name,". This can be used to select a subset for WGS or 
                         transcriptomics sequencing based on 16S profiles. Support Vector Machine is built using previously selected features,
                         and",n.select,"samples corresponding to each of the two levels
                         of the grouping variable are picked. Samples are picked as 
                         predicted correctly by the linear SVM after applying to the frequency 
                         profiles the 
                         inverse hyperbolic sign transform and normalizing across columns, 
                         and selecting samples that are maximally
                         distant from the SVM separating plane. The nuisanse parameter
                         of the SVM is picked through a grid search maximizing
                         prediction accuracy in training with resampling.
                         Accuracy of the final model is reported both on the
                         training set and with cross-validation. Several
                         performance charts are shown for the training set.
                         After selecting samples, the same metrics are shown
                         for the selected samples only. Additionally, the
                         abundance profiles are plotted for the selected samples.
                         Typically, you might expect a fairly poor
                         accuracy for the full sample set, and very good - for the
                         selected samples assuming that the number of selected
                         samples is a small fraction of the total. For reference, a purely
                         random binary classifier would result in a 50% accuracy."))
  
  sample.group = m_a$attr[,sample.group.name]
  
  ## nasty "feature" - when species.sel is a factor,
  ## the subscripting of count seems to use the integer level of a factor
  species.sel = as.character(species.sel)
  
  m = m_a$count[,species.sel]
  
  #If building kernel from a distance:
  #S = np.exp(-D * gamma), where one heuristic for choosing gamma is 1 / num_features
  #or
  #S = 1. / (D / np.max(D))
  #m = norm.boxcox(m)
  
  m = ihs(m)
  
  m = decostand(m,method="standardize",MARGIN=2)
  
  report$add.vector(colnames(m),name="feature",caption="Using these features for sample selection.")
  report$add.header("Models trained on the full dataset and their performance")  
  ##make.global(species.sel)
  ##make.global(m)
  ##make.global(sample.group)
  bootControl <- trainControl(number = 40)
  #set.seed(2)
  scaled = F
  sigma = sigest(m,scaled=scaled)[2]
  ##print(paste("sigma",sigma))
  ##using prob.mod=T screws up the model - the plane is no longer at 0 decision
  ##threshold. Maybe bias is lost?
  if(T) {
    cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
    registerDoSNOW(cl)      
    #class.weights=1/table(sample.group),
    #tuneGrid = expand.grid(sigma=sigma,C = 4**(-2:7)),
    mod.fit <- train(m, sample.group,
                     method = "svmLinear", #"svmRadial" 
                     prob.model = F,
                     cross=3,
                     class.weights=1/table(sample.group),
                     tuneGrid = expand.grid(C = 4**(-2:7)),
                     trControl = bootControl, scaled = scaled,
                     metric="Accuracy")  #"Kappa"
    stopCluster(cl)
    report$add.printed(mod.fit,caption="Results of SVM parameter fitting")
    mod = mod.fit$finalModel
  }
  else {
    for (C in c(1)) {
      mod = ksvm(x=m, y=sample.group, kernel = "vanilladot",
                 kpar = "automatic", C = C, cross = 6, prob.model = F,
                 class.weights=1/table(sample.group),scaled=scaled)
      #print(cross(mod))
      print(mod)
    }
  }
  
  report$add.printed(mod,caption="Best model trained on the full input set")  
  
  #mod.pred.prob = predict(mod, m, type = "probabilities")
  ##make.global(mod.pred.prob)
  #print(sum(mod.pred.prob[,"Control"]>0.5))
  mod.pred.score = predict(mod,m,type="decision")
  mod.pred = predict(mod,m)
  ##make.global(mod.pred.score)
  
  report$push.section(report.section)
  show.pred.perf(mod.pred,mod.pred.score,sample.group)  
  report$pop.section()
  
  #}
  report$add.header("Performance of the models and plots for selected samples")  
  sample.ids = rownames(m)
  cut.res = cut.top.predictions(mod.pred.score,
                                sample.group,
                                sample.ids,
                                n.cut=n.select,
                                filter.by.label=filter.by.label)
  ids.sel = c(plyr::laply(cut.res,function(x) x[["ids"]]))
  mask.sel = sample.ids %in% ids.sel
  ids.sel = sample.ids[mask.sel]
  pred.score.sel = mod.pred.score[mask.sel]
  pred.sel = mod.pred[mask.sel]
  sample.group.sel = sample.group[mask.sel]
  samples.sel = data.frame(SampleID=ids.sel,Group=sample.group.sel,Score=pred.score.sel)
  samples.sel = samples.sel[order(samples.sel$Score),]
  report$add.table(samples.sel,
                   caption="Selected samples, sorted by score. You can pick subsets at the opposite score extremes.")
  samples.all = data.frame(SampleID=sample.ids,Group=sample.group,Score=mod.pred.score)  
  samples.all = samples.all[order(samples.all$Score),]
  fn = report$make.file.name("sample.selection.tsv")
  write.table(samples.all,fn,sep="\t",row.names = FALSE)
  report$add.descr(paste(
    "All sample IDs with decision score and grouping variable are saved into file ",
    fn,"that you can use to fine tune the selection."))
  report$push.section(report.section)
  show.pred.perf(pred.sel,pred.score.sel,sample.group.sel)
  
  report$add.header("Abundance Plots")
  
  m_a.sel = m_a
  
  m_a.sel$attr = m_a.sel$attr[mask.sel,]
  
  m_a.sel$count = m_a.sel$count[mask.sel,species.sel]
  
  pl.hist = plot.abund.meta(m_a=m_a.sel,
                            id.vars=c(sample.group.name),
                            geom="bar",
                            id.var.dodge=sample.group.name
  )$plot
  #env=as.environment(as.list(environment(), all.names=TRUE))
  #print(names(as.list(env)))
  #print(evals("pl.hist",env=env))
  report$add(pl.hist,
             caption=paste("Abundance profile of samples maximally different with regard to ",
                           sample.group.name,"(only features that were used for selection are shown)")
  )
  
  m_a.sel$count = m_a$count[mask.sel,]
  
  pl.hist = plot.abund.meta(m_a=m_a.sel,
                            id.vars=c(sample.group.name),
                            geom="bar",
                            id.var.dodge=sample.group.name
  )$plot
  #env=as.environment(as.list(environment(), all.names=TRUE))
  #print(names(as.list(env)))
  #print(evals("pl.hist",env=env))
  report$add(pl.hist,
             caption=paste("Abundance profile of samples maximally different with regard to ",
                           sample.group.name,"(the most abundant features are shown, even those not used for selection)")
  )
  report$pop.section()
}

proc.t1d.som <- function() {
  library(kononen)
  m.som = som(m, grid = somgrid(5, 5, "hexagonal"))
  m.pred = predict(m.som,newdata=m,trainX=m,trainY=m_a$attr$T1D)
  table(m_a$attr$T1D,m.pred$prediction)
  #plot assignment of original data
  plot(m.som,type="mapping",col=as.numeric(m_a$attr$T1D=="T1D")+1,pch=1)
  plot(m.som,type="mapping",col=as.numeric(m_a$attr$Batch),pch=2)
}


## Build mixed effects model for count data of a single feature
## as a function of attr metadata var (e.g. patient-control classification)
## (Binomial regression model. Add a sample specific random effect
## to account for overdispersion.)
## Value: p value of testing against a null hypothesis that 
## model coefficients are zero (ignoring intercept; two-sided)
test.counts.glmer.col <- function(taxa.count,
                                  attr,
                                  taxa.count_rowsum,
                                  formula_rhs,
                                  linfct) {
  
  
  ## T1D groups samples into patient and control groups (fixed effect).
  ## FamilyID groups patients and their control siblings into 
  ## matched groups (random effect).
  ## We add the unique ID of each observation (SampleID) as a random effect, in order 
  ## to account for overdispersion.
  ## This approach has been used in a T1D study [10.1371/journal.pone.0025792]:
  ## (http://www.plosone.org/article/info%3Adoi%2F10.1371%2Fjournal.pone.0025792)
  ## Specific implementations are described in
  ## (http://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html)
  ## and
  ## (http://depts.washington.edu/cshrb/wordpress/wp-content/uploads/2013/04/Tutorials-Tutorial-on-Count-Regression-R-code.txt)
  ## 
  ##ANOVA comparison for models with and without the overdispersion term in the
  ## case of Streptococcus is below:
  ##> anova(dat.glm.0,dat.glm)
  ##Data: dat
  ##Models:
  ##  dat.glm.0: y ~ T1D + (1 | FamilyID)
  ##dat.glm: y ~ T1D + (1 | FamilyID) + (1 | SampleID)
  ##Df    AIC    BIC   logLik deviance  Chisq Chi Df Pr(>Chisq)    
  ##dat.glm.0  3 5837.3 5846.1 -2915.66   5831.3                             
  ##dat.glm    4 1287.0 1298.7  -639.49   1279.0 4552.4      1  < 2.2e-16 ***
  ##  ---
  ##  Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1  
  response = cbind(taxa.count,taxa.count_rowsum-taxa.count)
  test_glmer_inner <- function() {
    
    ##default optimizer was not converging often, 
    ##the developer recommended using bobyga:
    ##http://stackoverflow.com/questions/21344555/convergence-error-for-development-version-of-lme4
    
    dat.glmer = glmer(as.formula(paste("response",formula_rhs,sep="~")),
                      data=attr,
                      family="binomial",
                      control=glmerControl(optimizer="bobyqa")
                      ## need to load package "optimx" for the optimizer below
                      #control=glmerControl(optimizer="optimx",
                      #                    optCtrl=list(method="L-BFGS-B") #"nlminb"
                      #)
    )
    ## For some reason, T1D name soemtimes becomes T1DT1D in the output - it looks like
    ## there is a bug in glmer, and glmer corrupts something in its internal state,
    ## and generates dummy indicator variables for the factors with names build from
    ## contrasts values rather than the original factor levels, when the whole
    ## process is executed again in the same R session. Run the whole analysis each time
    ## in a fresh R session.
    
    #p_val = summary(dat.glmer)$coefficients[,"Pr(>|z|)"]["T1DT1D"]
    ## If we have multiple fixed effects factors, then the right way to obtain
    ##p-vals is below, using multiple testing correction
    #library(multcomp)
    dat.glht = glht(dat.glmer,linfct=linfct)
    dat.glht.summ = summary(dat.glht)
    p_val = as.numeric(dat.glht.summ$test$pvalues)
    p_val
  }
  
  #warn_saved = options()$warn
  #options(warn=2)
  ok=T
  #somehow if both warning and error arguments are used
  #in a single call to tryCatch, then error handler is not
  #called
  p_val = tryCatch(tryCatch(test_glmer_inner(),
                            
                            warning=function(w) {
                              warning(paste("Warning caught in glmer; results converted to NA: ",w,"\n"))
                              ok<<-F
                            }
  ),                   error=function(e) {
    warning(paste("Error caught in gmler; results converted to NA: ",e,"\n"))
    ok<<-F
  }
  )
  
  #options(warn=warn_saved)
  
  if (!ok) {
    
    print("Warnings or errors in glmer, returning NA")
    
    #cat(xtabs(~(response[,"taxa.count"]>0)+attr$Batch+attr$T1D))
    
    p_val = NA
    
  }
  #print(paste("p_val=",p_val))
  p_val
}

## Test all features pairwise using Poisson regression model with 
## a subject specific random effect
test.counts.glmer <- function(m_a,
                              formula_rhs,
                              linfct,
                              alpha=0.05,
                              p_adjust_method="BH") {
  ##http://thebiobucket.blogspot.com/2011/06/glmm-with-custom-multiple-comparisons.html
  library(lme4)
  library(multcomp)
  count = m_a$count
  count_rowsum = rowSums(count)
  #attr = m_a$attr[,c("T1D","FamilyID","SampleID","Batch")]
  attr = m_a$attr
  # Not working properly with .parallel=T in the presence
  # of warnings or errors - getting all NAs no matter what
  # I try for error handling in tryCatch
  ##TODO: try foreach
  p_vals = plyr::aaply(count,
                       2,
                       test.counts.glmer.col,
                       attr,
                       count_rowsum,
                       formula_rhs=formula_rhs,
                       linfct=linfct,
                       .inform=F,
                       .parallel=F,
                       .paropts=list(.packages=c("lme4","multcomp")))
  names(p_vals) = colnames(count)
  p_vals = p_vals[order(p_vals,na.last=T)]
  p_vals_signif = p_vals[!is.na(p_vals) & p_vals<=alpha]
  p_vals_adj = p.adjust(p_vals,method=p_adjust_method)
  p_vals_adj_signif = p_vals_adj[!is.na(p_vals_adj) & p_vals_adj<=alpha]
  ret = llist(p_vals,
              p_vals_adj,
              p_vals_signif,
              p_vals_adj_signif,
              formula_rhs,
              linfct,
              alpha,
              p_adjust_method)
  class(ret) = "taxa_count_glmer"
  
  return(ret)
}

print.taxa_count_glmer <- function(x,...) {
  cat ("Test results from fitting mixed effects binomial model for each taxa\n")
  cat ("Right-hand side of formula:\n")
  cat (x$formula_rhs,"\n")
  cat ("Specification of the linear hypotheses (see glht):\n")
  cat (x$linfct,"\n")
  cat ("p-values that pass significance cutoff before multiple testing correction:\n")
  print(x$p_vals_signif)
  cat ("p-values that pass significance cutoff after BH multiple testing correction:\n")
  print(x$p_vals_adj_signif)
}

report.counts.glmer <- function(report,x,...) {
  report$add.p("Test results from fitting mixed effects binomial model for each feature")
  report$add.p("Right-hand side of formula:")
  report$add.p(x$formula_rhs)
  report$add.p("Specification of the linear hypotheses (see glht):")
  report$add.p(x$linfct)
  
  ## our default table style of "grid" somehow screws up row representation, so we use "rmarkdown" here
  report$add.vector(x$p_vals_signif,"p_value",
                    caption="p-values that pass significance cutoff before multiple testing correction.",
                    style="rmarkdown")
  report$add.vector(x$p_vals_adj_signif,"p_value",
                    caption="p-values that pass significance cutoff after BH multiple testing correction.",
                    style="rmarkdown")
  report$add.vector(x$p_vals_adj[1:min(20,length(x$p_vals_adj))],"p_value",
                    caption="Top 20 p-values after BH multiple testing correction.",
                    style="rmarkdown")
  failed.names = names(x$p_vals[is.na(x$p_vals)])
  if(length(failed.names)>0) {
    report$add(list(failed.names),
               caption="features for which the model could not be built for any reason")
  }
}

## Wilcoxon rank-sum effect size r taken from Andy Fields
## two-sided
wilcox.eff.size.r <-function(pval, n){
  z<- qnorm(pval/2)
  r<- z/ sqrt(n)
  return(abs(r))
}

## x is the sample matrix, with samples in COLUMNS
##Note: Do not pass here statistic from coin::wilcox_test, it computes some other statistic, even when using statistic(wt,"linear"), we
##are getting values that are larger than the product of group sizes

wilcox.eff.size <- function(x,stat,pval=NULL,group=NULL,type="unpaired") {
  
  if(type == "unpaired") {
    nsamp = ncol(x)
    group = factor(group)
    stopifnot(nsamp==length(group))
    taby <- table(group)
    stopifnot(nlevels(group) == 2)
    npairs = prod(taby)
    common.lang = stat/npairs
    ##statistics from RankingWilcoxon is for the group with most samples,
    ##but we want in the order of the levels
    if(taby[2]>taby[1]) {
      common.lang = 1 - common.lang
    }
  }
  else if(type == "onesample") {
    ## assumes zero fudge was done when computing p-values
    nsamp = rowSums(x!=0)
    common.lang = (rowSums(x > 0) + 0.5*rowSums(x==0))/
      ncol(x)
  }
  else {
    stop("type parameter must be one of c('unpaired','onesample')")
  }
  rbs = 2*common.lang - 1
  r = NULL
  if(!is.null(pval)) {
    r = wilcox.eff.size.r(pval,nsamp)
  }
  return (data.frame(common.lang.eff.size=common.lang,
                     rank.biserial.corr.eff.size=rbs,
                     r.eff.size=r)
  )
}

## GeneSelector Wilcoxon signed rank implementation does not work correctly with zeros
## (possibly because it was designed for microarrays where ties never happen?).
## Here we: 
## - fix treatment of zeros in calculating statistic by using "zero fudge"
## - call exactRankTests::wilcox.exact by default when p-values are requested.
## wilcox.exact implements exact p-value computation (by default, at n <50)
## even when ties are present.
## Our statistic now is computed in the same way as in wilcox.exact.
## Call RankingWilcoxonAT(...,pvalues=F) in GeneSelector ranking replications,
## and call RankingWilcoxonAT(...,pvalues=T) to get p-values just once on the
## full dataset (~10x slower than statistic only).

RankingWilcoxonAT <- function (x, y, type = c("unpaired", "paired", "onesample"), 
                               pvalues = FALSE, gene.names = NULL, 
                               pvalues.impl="wilcox.exact",...) 
{
  mode(x) <- "numeric"
  if (length(y) != ncol(x)) 
    stop("Length of y is not equal to the number of columns of the expression matrix \n.")
  type <- match.arg(type)
  if (!is.element(type, eval(formals(RankingWilcoxonAT)$type))) 
    stop("Argument 'type' must be one of 'unpaired', 'paired' or 'onesample'. \n")
  y <- factor(y)
  if (type == "unpaired") {
    if (nlevels(y) != 2) 
      stop("Type has been chosen 'unpaired', but y has not exactly two levels ! \n")
    taby <- table(y)
    levy <- names(taby)[which.max(taby)]
    ind <- (y == levy)
    Rx <- apply(x, 1, rank)
    r1 <- colSums(Rx[ind, , drop = FALSE]) - sum(1:max(taby))
    e1 <- taby[1] * taby[2]/2
    if (pvalues) {
      if(pvalues.impl=="genesel") {
        maxr <- sum((min(taby) + 1):length(y)) - sum(1:sum(ind))
        pvals <- 2 * pwilcox(ifelse(r1 < e1, maxr - r1, r1), 
                             taby[1], taby[2], lower.tail = FALSE)
      }
      else {
        require("exactRankTests")
        pvals = apply(x,1,function(z) {
          wilcox.exact(x=z[ind],y=z[!ind],paired=F)$p.value
        })
      }
      
    }
    else pvals <- rep(NA, nrow(x))
  }
  if (type == "paired") {
    tab <- table(y)
    if (length(tab) != 2) 
      stop("Type has been chosen 'paired', but y has not exactly two levels. \n")
    xx1 <- x[, 1:tab[1]]
    xx2 <- x[, -c(1:tab[1])]
    if (tab[1] != tab[2] || length(unique(y[1:tab[1]])) != 
        1 | length(unique(y[-c(1:tab[1])])) != 1) 
      stop("Incorrect coding for type='paired'. \n")
    diffxx <- xx2 - xx1
    r1 <- apply(diffxx, 1, function(z) {
      ##zero fudge
      z = z[z!=0]
      zz <- rank(abs(z))
      sum(zz[z > 0])
    })
    ##zero fudge
    ly <- rowSums(diffxx!=0)
    e1 <- (ly/2) * (ly/2 + 1)/4
    if (pvalues) {
      if(pvalues.impl=="genesel") {
        maxr = (1+ly/2)*ly/4
        pvals <- 2 * psignrank(ifelse(r1 < e1, maxr - r1, 
                                      r1), n = ly/2, lower.tail = FALSE)
      }
      else {
        require("exactRankTests")
        pvals = apply(diffxx,1,function(z) {
          wilcox.exact(x=z,y=NULL,paired=F)$p.value 
        })
      }
    }
    else pvals <- rep(NA, nrow(x))
  }
  if (type == "onesample") {
    if (length(unique(y)) != 1) 
      warning("Type has been chosen 'onesample', but y has more than one level. \n")
    
    r1 <- apply(x, 1, function(z) {
      ##zero fudge
      z = z[z!=0]
      zz <- rank(abs(z))
      sum(zz[z > 0])
    })
    ##zero fudge
    ly <- rowSums(x!=0)
    e1 <- (ly) * (ly + 1)/4
    maxr = (1+ly)*ly/2
    if (pvalues) {
      if(pvalues.impl=="genesel") {
        pvals <- 2 * psignrank(ifelse(r1 < e1, maxr - r1, 
                                      r1), n = ly, lower.tail = FALSE)
      }
      else {
        require("exactRankTests")
        pvals = apply(x,1,function(z) {
          wilcox.exact(x=z,y=NULL,paired=F)$p.value
        })
      }
    }
    else pvals <- rep(NA, nrow(x))
  }
  statistic <- r1
  ranking <- rank(-abs(r1 - e1), ties.method = "first")
  if (!is.null(gene.names)) 
    names(pvals) <- names(statistic) <- gene.names
  else {
    if (!is.null(rownames(x))) 
      names(pvals) <- names(statistic) <- rownames(x)
  }
  new("GeneRanking", x = x, y = as.factor(y), statistic = statistic, 
      ranking = ranking, pval = pvals, type = type, method = "WilcoxonAT")
}

## We have to patch this in order to add RankingWilcoxonAT to hard-wired
## list of methods
RepeatRankingAT <- function (R, P, scheme = c("subsampling", "labelexchange"), iter = 10, 
                             varlist = list(genewise = FALSE, factor = 1/5), ...) 
{
  scheme <- match.arg(scheme)
  if (!is.element(scheme, c("subsampling", "labelexchange"))) 
    stop("'scheme' must be  either 'subsampling' or 'labelexchange'")
  x <- R@x
  y <- R@y
  Pm <- P@foldmatrix
  type <- R@type
  iter <- ncol(Pm)
  rankm <- pvalm <- statisticm <- matrix(nrow = nrow(x), ncol = iter)
  rankfun <- switch(R@method, ordinaryT = RankingTstat, WelchT = RankingWelchT, 
                    BaldiLongT = RankingBaldiLong, Bstat = RankingBstat, 
                    Ebam = RankingEbam, Foldchange = RankingFC, FoxDimmicT = RankingFoxDimmic, 
                    Limma = RankingLimma, Permutation = RankingPermutation, 
                    Sam = RankingSam, ShrinkageT = RankingShrinkageT, SoftthresholdT = RankingSoftthresholdT, 
                    WilcEbam = RankingWilcEbam, Wilcoxon = RankingWilcoxon, 
                    WilcoxonAT = RankingWilcoxonAT)
  if (scheme == "subsampling") {
    for (i in 1:iter) {
      currx <- x[, Pm[, i]]
      curry <- y[Pm[, i]]
      repet <- rankfun(currx, curry, type, ...)
      rankm[, i] <- repet@ranking
      pvalm[, i] <- repet@pval
      statisticm[, i] <- repet@statistic
    }
  }
  if (scheme == "labelexchange") {
    ly <- levels(y)
    nly <- nlevels(y)
    if (nly != 2) 
      stop("scheme 'labelexchange' not allowed if y has only one level \n")
    for (i in 1:iter) {
      curry <- y
      curry[!Pm[, i]] <- ifelse(y[!Pm[, i]] == ly[1], ly[2], 
                                ly[1])
      repet <- rankfun(x, curry, type, ...)
      rankm[, i] <- repet@ranking
      pvalm[, i] <- repet@pval
      statisticm[, i] <- repet@statistic
    }
  }
  colnames(rankm) <- colnames(pvalm) <- colnames(statisticm) <- paste("iter", 
                                                                      1:iter, sep = ".")
  new("RepeatedRanking", original = R, rankings = rankm, pvals = pvalm, 
      statistics = statisticm, scheme = scheme)
}


new_genesel <- function(...) {
  x = new_mgsatres(...)
  class(x) <- append(class(x),"genesel",0)
  return(x)
}

## It is assumed that the count matrix x is transformed/normalized already
## e.g. with decostand(ihs(count),method="standardize",MARGIN=2),
## otherwise set tran_norm to TRUE and the command above will be
## applied
genesel.stability <- function(m_a,
                              group.attr,
                              block.attr=NULL,
                              type="unpaired",
                              replicates=400,
                              samp.fold.ratio=0.5,
                              maxrank=20,
                              comp.log.fold.change=F,
                              ret.data.contrasts=T
) {
  library(GeneSelector)
  
  type.orig = type
  if(type=="paired") {
    type="onesample"
    s.c = sample.contrasts(m_a, group.attr = group.attr, block.attr = block.attr)
    levels.last.first = names(s.c$contrasts)
    m_a.c = s.c$m_a.contr
    x = m_a.c$count
    y.relev = rep(1,nrow(x))
  }
  else {
    x = m_a$count
    y = factor(m_a$attr[,group.attr])
    stopifnot(length(levels(y))==2)
    ##make last level to be first, so that effect sizes and statistics
    ##are computed for last over first
    y.relev = relevel(y,levels(y)[length(levels(y))])
    levels.last.first = levels(y.relev)
  }
  
  x = t(x)
  
  n_feat = nrow(x)
  n_samp = ncol(x)
  
  ranking_method = RankingWilcoxonAT
  
  rnk = ranking_method(x,y.relev,type=type,pvalues=T)
  
  rnk.vals = toplist(rnk,n_feat)
  rnk.vals$ranking = seq(nrow(rnk.vals))
  rnk.vals[rnk.vals$index,] = rnk.vals
  rnk.vals$pval.adjusted = p.adjust(rnk.vals$pval,method="BH")
  rnk.vals = cbind(data.frame(name=rownames(x)),rnk.vals)
  
  rnk.vals = cbind(rnk.vals,
                   wilcox.eff.size(x=x,
                                   stat=rnk.vals$statistic,
                                   pval=rnk.vals$pval,
                                   group=y.relev,
                                   type=type))
  m_a.lfc.paired = NULL
  if(comp.log.fold.change) {
    rnk.vals = cbind(rnk.vals,
                     l2fc.group.mean=t(group.log.fold.change(m_a$count,m_a$attr[,group.attr],base=2))
    )
    group.mean = t(count.summary(m_a$count,mean,m_a$attr[,group.attr],
                                 format="matrix",group.prefix="mean"))
    rnk.vals = cbind(rnk.vals,group.mean)
    if(type.orig=="paired") {
      m_a.g = s.c$m_a.groups
      m_a.lfc.paired = contrasts.groups.log.fold.change(m_a.g)
      
      rnk.vals = cbind(rnk.vals,
                       l2fc.paired.median = plyr::aaply(
                         m_a.lfc.paired$count,
                         2,
                         median
                       )
      )
      group.median = t(count.summary(m_a.g$count,function(x) median(x,na.rm = T),m_a.g$attr[,group.attr],
                                     format="matrix",group.prefix="median.paired"))
      rnk.vals = cbind(rnk.vals,group.median)
    }
  }
  ##make.global(rnk.vals)
  ##make.global(y.relev)
  ##make.global(type)
  
  if(replicates > 0) {
    #minclassize = min(table(y.relev))*samp.fold.ratio*0.9
    ##TODO: expose minclassize and balanced to user. But currently neither works.
    #contrary to the help page, does not work when y is a factor - need to convert
    fold_matr = GenerateFoldMatrix(y = as.numeric(y.relev), 
                                   k=trunc(n_samp*(1-samp.fold.ratio)),
                                   replicates=replicates,
                                   type=type)
    
    rep_rnk = RepeatRankingAT(rnk,fold_matr,scheme="subsampling",pvalues=F)
    ##make.global(rep_rnk)
    #toplist(rep_rnk,show=F)
    #stab_ovr = GetStabilityOverlap(rep_rnk, scheme = "original", decay = "linear")
    ### for a short summary
    #summary(stab_ovr, measure = "intersection", display = "all", position = 10)
    #summary(stab_ovr, measure = "overlapscore", display = "all", position = 10)
    ### for a graphical display
    #plot(stab_ovr)
    #aggr_rnk = AggregateSimple(rep_rnk, measure="mode")
    aggr_rnk = AggregateMC(rep_rnk, maxrank=n_feat)
    #aggr_rnk = AggregateSVD(rep_rnk)
    ##make.global(aggr_rnk)
    #toplist(aggr_rnk)
    #gsel = GeneSelector(list(aggr_rnk), threshold = "BH", maxpval=0.05)
    gsel = GeneSelector(list(aggr_rnk), threshold = "user", maxrank=n_feat)
    #show(gsel)
    #str(gsel)
    #toplist(gsel)
    selected = SelectedGenes(gsel)
    #pvals are always NA somehow, and we do not need the 'index' in
    #the output, so not adding 'selected' to the output
    index = selected$index
  }
  else {
    gsel = NULL
    #index = order(rnk.vals$ranking)
    index = order(rnk.vals$pval.adjusted,rnk.vals$pval)
  }
  rnk.vals$index = NULL
  rnk.vals$ranking = NULL
  rnk.vals = rnk.vals[index,]
  ret = new_genesel(stab_feat=rnk.vals,
                    gsel=gsel,
                    levels.last.first=levels.last.first,
                    n.feat = n_feat,
                    n.samp = n_samp)
  if(type!="unpaired") {
    if(ret.data.contrasts) {
      ret$m_a.contrasts = m_a.c
      ret$m_a.lfc.paired = m_a.lfc.paired
    }
    ret$contrasts = s.c$contrasts
  }
  return(ret)
}

feat.sel.samr <- function(m_a.abs) {
  library(samr) #Tibshirani's package for feature selection in microarrays and RNASeq
  ##TODO: study assumptions of this method on sequence counts. The help page mentions only RNASeq.
  ##Something is probably not right because it only reports "genes down" and empty for "genes up" for
  ##T1D genus count data
  samfit = SAMseq(t(m_a.abs$count),m_a.abs$attr$T1D,resp.type="Two class unpaired",geneid=colnames(m_a.abs$count))
  print(samfit)
  plot(samfit)
}

make.selection.mask.matrix <- function(m.source,val=T) {
  matrix(val,nrow(m.source),ncol(m.source))
}

make.matrix.for.subscript.testing <- function(nrow=3,ncol=4) {
  m = matrix("",nrow,ncol)
  for(irow in 1:nrow) for(icol in 1:ncol) m[irow,icol] = sprintf("c(%s,%s)",irow,icol)
  m
}

make.mask.for.selected.subscript.testing <- function(m.source,m.selected) {
  mask = make.selection.mask.matrix(m.source,F)
  m.ind.sel=foreach(x=as.character(m.selected),.combine=rbind) %do% {
    eval(parse(text=x))
  }
  print(m.ind.sel)
  mask[m.ind.sel] = T
  mask
}

rank.biserial.corr <- function(x,y) {
  rx = rank(c(x,y),na.last=NA,ties.method="average")
  n.x = length(x)
  n.y = length(y)
  comm.lang = (sum(rx[1:n.x]) - n.x*(n.x+1)/2)/(n.x*n.y)
  return (2*comm.lang-1)
}

test.dist.matr.within.between <- function(m_a,
                                          group.attr,
                                          block.attr,
                                          n.perm=4000,
                                          dist.metr="bray",
                                          col.trans="range",
                                          data.descr="proportions of counts",
                                          norm.count.task=NULL) {
  require(vegan)
  require(permute)
  
  report$add.header(sprintf('Comparison and test of significant difference for profile
dissimilarities within and between blocks defined
by attribute %s across groups defined by attribute %s',block.attr,group.attr))
  
  s.c = sample.contrasts(m_a, group.attr = group.attr, block.attr = block.attr)
  
  m_a = s.c$m_a.groups
  
  if(!is.null(norm.count.task)) {
    m_a <- norm.count.report(m_a,
                             descr="Dist.Matr.Within.Between",
                             norm.count.task=norm.count.task)
  }
  
  ##Negative values break bray-curtis and jaccard distances; we standardize to "range" to reduce
  ##the influence of very abundant species:
  
  if(!is.null(col.trans) && col.trans != "ident") {
    m_a$count = decostand(m_a$count,method=col.trans,MARGIN=2)
    col.trans.descr = sprintf(" Profile columns are normalized with %s method of decostand function.",col.trans)
  }
  else {
    col.trans.descr = ""
  }
  
  dist.metr.descr = sprintf(" Dissimilarity index is %s.",dist.metr)
  
  ##check that there is strictly one observation in each cell of (block,group)
  group.attr.lev = levels(factor(m_a$attr[,group.attr]))
  block.attr.lev = levels(factor(m_a$attr[,block.attr]))
  stopifnot(length(group.attr.lev)==2)
  fam.counts = as.matrix(with(m_a$attr,xtabs(as.formula(sprintf("~%s+%s",block.attr,group.attr)))))
  ##sort by keys
  m_a = subset.m_a(m_a,subset=order(m_a$attr[,block.attr],m_a$attr[,group.attr]))
  dd = as.matrix(vegdist(m_a$count,dist.metr))
  dd = dd[m_a$attr[,group.attr]==group.attr.lev[1],
          m_a$attr[,group.attr]==group.attr.lev[2]]
  
  block = m_a$attr[,block.attr,drop=F]
  
  dd.block.col = block[colnames(dd),]
  dd.block.row = block[rownames(dd),]
  
  mask.within = make.selection.mask.matrix(dd,F)
  
  for(irow in 1:nrow(dd)) {
    mask.within[irow,] = (dd.block.col == dd.block.row[irow])
  }
  
  st.obs = rank.biserial.corr(dd[!mask.within],dd[mask.within])
  
  n.col = ncol(dd)
  ## how() is masked by kernlab
  ## permute blocks ("plots"), otherwise keep the order unchanged
  ## All permute machinery refuses to work with unbalanced datasets
  ctrl = permute::how(
    within=Within(type="none"),
    plots=Plots(type="free",strata=dd.block.col)
  )
  ## default check=T gives some cryptic error message despite
  ## reasonably looking permutations being generated with check=F
  perm = permute::shuffleSet(n=n.col,
                             nset=n.perm-1,
                             control=ctrl,
                             check = F
  )
  perm = rbind(perm,1:n.col)
  
  
  st.perm = foreach(i.iter=1:n.perm,.combine=c,.export=c("rank.biserial.corr")) %dopar% {
    i.col = perm[i.iter,]
    mask.within.i = mask.within[,i.col]
    rank.biserial.corr(dd[!mask.within.i],dd[mask.within.i])
  }
  
  p.val = mean(st.perm>=st.obs)
  
  report$add.descr(paste(
    sprintf('%s%s The matrix of 
  dissimilarities D is formed where rows correspond to observations with level %s
  of %s, and columns - to level %s. The elements of 
  this matrix corresponding to rows and columns with the same 
  level of %s are called \"within\" block dissimilarities, while
  the elements drawn from combinations of rows and columns
  where %s are not equal are called \"between\" blocks dissimilarities.',
            dist.metr.descr,
            col.trans.descr,
            group.attr.lev[1],group.attr,group.attr.lev[2],
            block.attr,block.attr),
    'The null hypothesis is that the observed difference of \"between\" and \"within\" 
  block dissimilarities is consistent with what could be expected 
  if the block structure was assigned to the observations at random.',
    'The alternative hypothesis is that the \"between\"/\"within\" 
  difference is larger than would have been expected from a random block assignment.',
    sprintf('We simulate %s matrices in which both \"between\" and \"within\" 
  dissimilarities come from the null distribution 
  by permuting the %s labels of the columns
  of matrix D.',n.perm,block.attr), 
    'The rank biserial correlation (Grissom, R. J., and J. J. Kim. 
  \"Effect Sizes for Research: Univariate 
  and Multivariate Applications, 2nd Edn New York.\" NY: Taylor and Francis (2012)) is
  computed between the observed \"between\" and \"within\" dissimilarities both in the observed and
  simulated samples. Positive values of this correlation statistic would indicate 
  \"between\" dissimilarities being stochastically larger than \"within\" dissimilarities.
  The p-value is estimated as the fraction of the simulated statistic values that are as high or higher 
  than the observed value.',
    sprintf('The estimated p-value was %f and the observed value of the statistic was %f.',
            p.val,st.obs)))
  dd.pl = rbind(
    data.frame(x=as.numeric(dd[!mask.within]),
               group="Between"),
    data.frame(x=as.numeric(dd[mask.within]),
               group="Within")
  )
  
  g = show.distr.group(dd.pl$x,dd.pl$group)
  g = g + 
    #geom_vline(xintercept=st.obs,color="blue",size=rel(1.5),linetype="dashed") +
    #geom_vline(xintercept=st.ci[4],color="blue",size=rel(1),linetype="dotted") +
    #geom_vline(xintercept=st.ci[5],color="blue",size=rel(1),linetype="dotted") +
    xlab("Dissimilarity") +
    ylab("Empirical Distribution Density")
  report$add(g,caption=sprintf('Emprical distribution density plots of the 
             profile-profile
             dissimilarities observed between and within %s blocks.
             Distances were computed only across levels of %s variable.',
                               block.attr,group.attr))
  res = list(p.val=p.val,statistic=st.obs)
  return(res)
}

## Taken from (http://depts.washington.edu/cshrb/wordpress/wp-content/uploads/2013/04/Tutorials-Tutorial-on-Count-Regression-R-code.txt)
### Small utility function to get (conditional) rate ratios and
### 95% CI from fitted glmer() object -- only appropriate for 
### a Poisson model (though, we might do similar things with a
### binomial outcome)
lmerCI <- function(obj, rnd = 2) {
  cc <- fixef(obj)
  se <- sqrt(diag(vcov(obj)))
  upper <- cc + 1.96*se
  lower <- cc - 1.96*se
  out <- data.frame(RR = round(exp(cc), rnd), 
                    upper = round(exp(upper), rnd), 
                    lower = round(exp(lower), rnd))
  out
}

## Combine samples within blocks as sum of values within each block
## weighed by contrasts coefficients corresponding to levels of
## group variable. Contrasts must be a named vector with elements 
## matching levels of the group variable. Levels that are absent from
## the contrasts argument are weighed as zero.
## Samples within each block are first averaged for each level of the
## group variable. If any of the non-zero contrast levels is not matched
## (missing from data), the entire block is dropped.
## Default value of contrasts will be set to result in (last-first) levels.
## Value: named list with m_a objects for contrasts and balanced observations
## attr values will have just one record picked from e.g. a pair that formed
## the contrast; therefore, if nay attr values are used later, they must be
## already indetical within the grouped original records. No checks are made here 
## about that requirement.
sample.contrasts <- function(m_a,group.attr,block.attr,contrasts=NULL,return.groups=T) {
  #m_m = melt(m_a$count,varnames=c("SampleID","feature"),value.name="cnt")
  #attr_sub = m_a$attr[,c("SampleID",pair.attr)]
  #m_m = join(m_m,attr_sub,by="SampleID",type="inner",match="first")
  attr.names = c(block.attr,group.attr)
  attr = m_a$attr[,attr.names]
  if(is.null(contrasts)) {
    lev = levels(factor(attr[,group.attr]))
    #default contrast is like in DESeq2 - (last - first) levels
    lev = c(lev[length(lev)],lev[1])
    contrasts = c(1,-1)
    names(contrasts) = lev
  }
  contrasts.ret = contrasts
  #ord = order(attr[,block.attr],attr[,group.attr])
  #dat = cbind(attr[ord,],m_a$count[ord,])
  dat = cbind(attr,m_a$count)
  attr.mask = names(dat) %in% attr.names 
  mean.counts = function(x) {
    colMeans(x[,!attr.mask,drop=F])
    #cbind(x[1,attr.mask,drop=F],colMeans(x[,!attr.mask,drop=F]))
  }
  dat = plyr::ddply(dat,c(block.attr,group.attr),mean.counts)
  
  ## convert contrasts to data.frame for joining
  contrasts = data.frame(names(contrasts),contrasts)
  contr.attr = ".contrast"
  names(contrasts) = c(group.attr,contr.attr)
  
  ## this will drop all group levels not in contrasts
  dat = plyr::join(dat,contrasts,by=group.attr,type="inner")
  
  if(return.groups) {
    dat.groups = dat
  }
  
  attr.mask = names(dat) %in% c(attr.names,contr.attr)
  dat = cbind(dat[,block.attr,drop=F],
              dat[,!attr.mask,drop=F]*dat[,contr.attr,drop=T])
  
  attr.mask = names(dat) %in% c(attr.names,contr.attr)
  
  n.contr = nrow(contrasts)
  
  ## if() will drop all blocks where not all contrasts matched
  dat = plyr::ddply(dat,c(block.attr),
                    function(x) {
                      if(nrow(x) == n.contr) {
                        colSums(x[,!attr.mask,drop=F])
                      }
                      else {
                        NULL
                      }
                    }
  )
  
  rownames(dat) = dat[,block.attr]
  
  attr.mask = names(dat) %in% c(attr.names,contr.attr)
  m_a.contr = list(count=as.matrix(dat[,!attr.mask]),
                   attr=plyr::join(dat[,block.attr,drop=F],
                                   m_a$attr,
                                   by=block.attr,
                                   match="first")
  )
  with(m_a.contr, stopifnot(all(rownames(count)==rownames(attr))))
  
  if(return.groups) {
    ##only groups with all requested contrast levels
    dat.groups = dat.groups[dat.groups[,block.attr] %in% rownames(m_a.contr$count),]
    rownames(dat.groups) = plyr::maply(dat.groups[,attr.names],paste,sep=".",.expand=F)
    attr.mask = names(dat.groups) %in% c(attr.names,contr.attr)
    m_a.groups = list(count=as.matrix(dat.groups[,!attr.mask]),
                      attr=plyr::join(dat.groups[,attr.mask,drop=F],
                                      m_a$attr,
                                      by=attr.names,
                                      match="first"))
    with(m_a.groups, stopifnot(all(rownames(count)==rownames(attr))))
  }
  else {
    m_a.groups = NULL
  }
  return (list(m_a.contr=m_a.contr,m_a.groups=m_a.groups,contrasts=contrasts.ret))
}

contrasts.groups.log.fold.change <- function(m_a.g,base=2,
                                             contrasts=c(1,-1),
                                             offset=.Machine$double.eps*100) {
  lfc = log(m_a.g$count[m_a.g$attr$.contrast==contrasts[1],]+offset,base=base) - 
    log(m_a.g$count[m_a.g$attr$.contrast==contrasts[2],]+offset,base=base)
  m_a.lfc = list(count=lfc,
                 attr=m_a.g$attr[m_a.g$attr$.contrast==contrasts[1],])
  with(m_a.lfc, stopifnot(all(rownames(count)==rownames(attr))))
  
  return (m_a.lfc)
}

report.sample.count.summary <- function(m_a,meta.x.vars=c(),group.vars=NULL,
                                        show.sample.totals=T, show.sample.means=T,
                                        sub.report=T) {
  report.section = report$add.header("Summary of total counts per sample",section.action="push",sub=sub.report)
  
  m_a.summ=make.sample.summaries(m_a)
  
  if(show.sample.totals) {
    report$add.table(cbind(m_a.summ$count[,"count.sum",drop=F],m_a.summ$attr),
                     caption="Total counts per sample",
                     show.row.names=F)
  }  
  
  if(show.sample.means) {
    report$add.vector(c(summary(m_a.summ$count[,"count.sum"])),caption="Summary of total counts per sample")
    report$add.table(plyr::ddply(join_count_df(m_a.summ),
                                 group.vars,
                                 plyr::summarise,
                                 Min.Count.Sum=min(count.sum),
                                 Max.Count.Sum=max(count.sum),
                                 Mean.Count.Sum=mean(count.sum),
                                 Median.Count.Sum=median(count.sum),
                                 Q25.Count.Sum=quantile(count.sum,0.25),
                                 Q75.Count.Sum=quantile(count.sum,0.75)
    ),
    caption="Group summaries of total counts per sample")
  }
  
  if(!is.null(group.vars)) {
    for(group.var in group.vars) {
      report$add(kruskal.test(m_a.summ$count[,"count.sum"],m_a.summ$attr[,group.var]),
                 caption=paste("Test for difference of total counts per sample across groups defined by",
                               group.var)
      )
      group.mean = count.summary(m_a$count,mean,m_a$attr[,group.var],format="matrix")
      report$add(friedman.test(t(group.mean)),
                 caption=paste("Test for difference in group means across all features, where groups
               are defined by",group.var))
      report$add.vector(rowMeans(group.mean),
                        caption=paste("Mean values of group means across all features, where groups
               are defined by",group.var))
    }
  }
  
  if(!(is.null(meta.x.vars) | is.null(group.vars))) {
    
    report$add.header("Iterating over meta data variables")
    report$push.section(report.section)
    
    for(x.var in meta.x.vars) {
      
      report$add.header("Iterating over group variables")
      report$push.section(report.section)
      
      for(group.var in group.vars) {
        show.sample.summaries.meta(m_a=m_a.summ,
                                   x.var=x.var,
                                   group.var=group.var)
      }
      report$pop.section()
    }
    
    report$pop.section()
    
  }
  
}


read.t1d.mg <-function(annot.type,level) {
  annot.dir = "../BATCH_01_02_META"
  mgrast.dir = paste(annot.dir,"BATCH_01-02_METAGENOMICS_MGRAST",sep="/")
  if (annot.type == "humann") {
    counts = read.humann.summary(paste(annot.dir,"humann/output/04b-keg-mpt-cop-nul-nve-nve.txt",sep="/"))
    counts = count.filter(counts,col_ignore=c(),min_max_frac=0,min_max=0,min_row_sum=0,other_cnt="other")      
  }
  else if (annot.type %in% c("cog","kegg","subsys")) {
    
    counts = read.mgrast.summary(paste(mgrast.dir,paste(annot.type,"tsv",sep="."),sep="/"),
                                 file_name.id.map=paste(mgrast.dir,"mgrast_to_samp_id.tsv",sep="/"))
    print(level)
    counts = mgrast.to.abund.df(counts,level)
    counts = count.filter(counts,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=0,other_cnt="other")
  }
  ##make.global(counts)
  
  meta = load.meta.t1d("annotation20130819.csv",as.merged=T)
  return (merge.counts.with.meta(counts,meta))
}

proc.t1d.mg <- function() {
  annot.types = c("humann","cog","subsys","kegg")
  #annot.types = c("subsys")
  do.std.plots = T
  do.tests = T
  for (annot.type in annot.types) {
    
    levels = annot.levels(annot.type)
    if(annot.type=="humann") {
      do.diversity = F
    }
    else {
      do.diversity = T
    }
    do.diversity = F
    
    for (level in levels) {
      
      label = paste(annot.type,level,sep=".")
      print (paste("Working on",label))
      
      count.meta = read.t1d.mg(annot.type,level)
      count.meta.data = count.meta$data
      count.meta.attr.names = count.meta$attr.names
      
      if (do.tests) {
        
        res.tests = tryCatchAndWarn({test.counts.t1d(count.meta.data,count.meta.attr.names,
                                                     label=label,
                                                     stability.transform.counts="ihs",
                                                     do.stability=T,
                                                     do.tests=T)})
      }
      if (do.std.plots) {
        tryCatchAndWarn({
          std.plots(count.meta.data,count.meta.attr.names,id.vars.list=
                      list(
                        c("T1D","age.quant"),
                        #c("Family","T1D"),
                        #c("SampleID","Batch"),
                        c("T1D"),
                        c("T1D","Batch")
                      ),
                    label=label,
                    do.diversity=do.diversity,
                    res.tests=res.tests
          )
        })
      }
    }
  }
}



load.meta.mr_oralc <- function(file_name) {
  
  #Sebastian's code from AnUnivariate.r
  #meta =read.csv(file_name, as.is=TRUE, header=TRUE,stringsAsFactors=T)
  meta =read.delim(file_name, header=TRUE,stringsAsFactors=T)
  row.names(meta) = meta$ID
  return (meta)
}


read.mr_oralc <- function(taxa.level=3) {
  #moth.taxa <- read.mothur.taxa.summary("X:/sszpakow/BATCH_03_16S/ab1ec.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary.txt")
  #moth.taxa <- read.mothur.taxa.summary("43aa6.files_x1.sorted.0.03.cons.tax.summary.otu.taxsummary.txt")
  moth.taxa <- read.mothur.taxa.summary("5e2c2.files_x1.sorted.0.03.cons.tax.summary.seq.taxsummary")
  taxa.lev.all = multi.mothur.to.abund.df(moth.taxa,taxa.level)
  taxa.lev = count.filter(taxa.lev.all,col_ignore=c(),min_max_frac=0,min_max=30,min_row_sum=500,other_cnt="other")
  #taxa.lev = taxa.lev.all
  meta = load.meta.mr_oralc("UCFTumorTissueSampleInventory.cleaned.txt")
  return (merge.counts.with.meta(taxa.lev,meta))
}



count.summary <- function(count,fun,group,format="data.frame",group.prefix=NULL) {
  library(data.table)
  .group = group
  .group = data.frame(.group)
  n.col.sel = 1:ncol(count)
  col.sel = colnames(count)
  dt = cbind(data.table(count),.group)
  by.col = colnames(.group)
  group.summ = as.data.frame(dt[,lapply(.SD,fun),by=by.col,.SDcols=col.sel])
  if(format == "data.frame") {
    return (group.summ)
  }
  else if(format == "matrix") {
    x = as.matrix(group.summ[,!names(group.summ) %in% names(.group),drop=F])
    rnames = apply(group.summ[,names(.group),drop=F],1,function(x) paste(x,collapse="."))
    if(!is.null(group.prefix)) {
      rnames = paste(group.prefix,rnames,sep=".")
    }
    
    rownames(x) = rnames
    return (x)
  }
  else stop(paste("Unknown format argument value",format))
}

group.mean.ratio <- function(count,group,row.names.pref="") {
  group.mean = count.summary(count,mean,group)
  stopifnot(dim(group.mean)[1] == 2)
  x = (group.mean[2,-1] / (group.mean[1,-1]+.Machine$double.eps))
  row.names(x) = paste(row.names.pref,paste(group.mean[2,1],group.mean[1,1],sep=".by."),sep=".")
  return (as.matrix(x))
}

## http://stackoverflow.com/questions/32513189/fast-and-elegant-way-to-calculate-fold-change-between-several-groups-for-many-va
fold.change <- function(mat,key,aggr.fun=mean,comb.fun=function(x,y) "/"(x,y),
                        out.format=c("long","wide")){
  library(purrr)
  out.format = out.format[1]
  key = as.data.frame(key)
  key.names = colnames(key)
  mat = as.matrix(mat)
  x <- data.frame(key,mat) %>%  slice_rows(key.names) %>% by_slice(map, aggr.fun)
  
  key.grouped = as.data.frame(x[,1:length(key.names),drop=F])
  x <- as.matrix(x[,-(1:length(key.names)),drop=F])
  
  # calculate changes between all rows
  i <- combn(nrow(x), 2)
  x <- comb.fun(x[i[1,],,drop=F] , x[i[2,],,drop=F])
  key.1 = key.grouped[i[1,],,drop=F]
  key.2 = key.grouped[i[2,],,drop=F]
  colnames(x) = colnames(mat)
  rownames(x) <- paste(
    plyr::maply(key.1,paste,sep=".",.expand=F),
    plyr::maply(key.2,paste,sep=".",.expand=F),
    sep = "-")
  rownames(key.1) = rownames(x)
  rownames(key.2) = rownames(x)
  
  ret = list(key.1=key.1,key.2=key.2,long=x)
  if(out.format == "wide") {
    stopifnot(ncol(key.1)==1)
    stopifnot(ncol(x)==1)
    key.grouped = key.grouped[,1]
    n = length(key.grouped)
    m = matrix(0,nrow=n,ncol=n,dimnames = list(key.grouped,key.grouped))
    for(i in 1:nrow(x)) {
      m[key.1[i,1],key.2[i,1]] = x[i,1]
      m[key.2[i,1],key.1[i,1]] = 1/x[i,1]
    }
    diag(m) = 1
    ret$wide = m
  }
  return (ret)
}

group.log.fold.change <- function(count,group,base=2) {
  row.names.pref = sprintf("l%sfc",base)
  return (log(group.mean.ratio(count=count,group=group,row.names.pref=row.names.pref)
              +.Machine$double.eps,base=base))
}

n.cases.pepe.cont <- function(tpr1,tpr0,r,fpr0,alpha,beta,eps,k=NULL) {
  theta = (qnorm(1-alpha) + qnorm(1-beta))**2
  k.opt = 1/r*sqrt((tpr1*(1-tpr1))/(fpr0*(1-fpr0)))
  if(is.null(k)) {
    k = k.opt
  }
  n.cases = theta * ( tpr1*(1-tpr1) + k*(r**2)*fpr0*(1-fpr0)) / ((tpr1 - tpr0)**2)
  n.controls.k = ceil(n.cases/k)
  n.controls.eps = ceil((qnorm(1-beta)/eps)**2 * fpr0 * (1-fpr0))
  n.controls = max(n.controls.k,n.controls.eps)
  n.total = n.cases + n.controls
  return(list(n.total=n.total,
              n.cases=n.cases,
              k.opt=k.opt,
              n.controls=n.controls,
              n.controls.eps=n.controls.eps,
              n.controls.k=n.controls.k,
              k=k,
              tpr1=tpr1,
              tpr0=tpr0,
              fpr0=fpr0,
              alpha=alpha,
              beta=beta,
              eps=eps))
}

deriv.ser <- function(x,y,x0) {
  library(pspline)
  predict(sm.spline(x, y), x0, 1)
}

marker.ver.power <- function(sm.df,tpr0,fpr0) {
  alpha=0.005
  beta=0.2
  eps=0.03
  k=0.28
  spe0=1-fpr0
  ss0 = sm.df[which(sm.df$specificities>=spe0)[1],]
  tpr1=ss0$sensitivities
  r = as.numeric(-deriv.ser(sm.df$specificities,sm.df$sensitivities,spe0))
  res = list(tpr1=tpr1,tpr0=tpr0,r=r,fpr0=fpr0,alpha=alpha,beta=beta,eps=eps,k=k)
  if(tpr1>tpr0) {
    res = n.cases.pepe.cont(tpr1=tpr1,tpr0=tpr0,r=r,fpr0=fpr0,alpha=alpha,beta=beta,eps=eps,k=k)
  }
  else {
    res$k = NULL #this actually deletes k, but res$k = k does not even if k is NULL
  }
  return(res)
}

show.partial.auc.roc <- function(response,predictor,predictor.descr,tpr0=0.25,fpr0=0.10) {
  spe.range=c(1,1-fpr0)
  se.range=NULL
  require(pROC)
  ro = pROC::roc(response, predictor, 
                 percent=F,
                 auc=F)
  smooth.bi = smooth(ro,method="bi")
  sm.df = with(smooth.bi,data.frame(specificities,sensitivities))
  pwr = marker.ver.power(sm.df,tpr0=tpr0,fpr0=fpr0)
  
  report$add({
    plot.roc(ro, 
             partial.auc=spe.range, 
             partial.auc.correct=TRUE, 
             # define a partial AUC (pAUC)  
             print.auc=TRUE, 
             #display pAUC value on the plot with following options:  
             print.auc.pattern=sprintf("Normalized partial AUC (%s-%s SP): %%.2f\n\
FPR cutoff: %.2f\n\
Smoothed TPR(FPR cutoff): %.2f",
                                       spe.range[1],spe.range[2],
                                       fpr0,
                                       pwr$tpr1
             ),
             print.auc.col="#1c61b6",  
             auc.polygon=TRUE, 
             auc.polygon.col="#1c61b6", 
             # show pAUC as a polygon  
             max.auc.polygon=TRUE, 
             max.auc.polygon.col="#1c61b622", 
             # also show the 100% polygon  
             main=NULL,
             grid.h=pwr$tpr1
    )
    #ro.si.se = ci.se(ro,specificities=spe.range[2])
    #plot(ro.si.se,col="green")
    plot.roc(smooth.bi,add=T,col="red")
    #plot.roc(smooth(ro,method="de",bw="SJ"),add=T,col="yellow")
    if(!is.null(se.range)) {
      plot(ro,
           add=TRUE, type="n", 
           # add to plot, but don't re-add the ROC itself (useless)  
           partial.auc=se.range, 
           partial.auc.correct=TRUE,  
           partial.auc.focus="se", 
           # focus pAUC on the sensitivity  
           print.auc=TRUE, 
           print.auc.pattern=sprintf("Corrected pAUC (%s-%s SE):\n%%.1f%%",
                                     se.range[1],se.range[2]),
           print.auc.col="#008600",  
           print.auc.y=40, 
           # do not print auc over the previous one  
           auc.polygon=TRUE, 
           auc.polygon.col="#008600",  
           max.auc.polygon=TRUE, 
           max.auc.polygon.col="#00860022"
      )   
    }
  },caption=sprintf("ROC curve for %s",predictor.descr))
  return (list(power=pwr,roc=roc,smooth.bi=smooth.bi))
}

counts.distro.report <- function(m_a,group.attr,descr) {
  report.section = report$add.header(sprintf("Empirical distributions of individual features %s",
                                             descr),
                                     section.action="push", sub=T)
  for(feat.name in colnames(m_a$count)) {
    report$add.header(paste("Feature name",feat.name))
    report$push.section(report.section)
    x = m_a$count[,feat.name]
    g = factor(m_a$attr[,group.attr])
    if(T) {
      report$add(show.distr.group(x,g),
                 caption=sprintf("Empirical distribution of %s grouped by %s",
                                 feat.name,group.attr))
      for(lev in levels(g)) {
        xl = x[g==lev]
        lev.descr = sprintf("level %s of %s",lev,group.attr)
        report$add.printed(summary(xl),
                           caption=sprintf("Summary of %s for %s",feat.name,lev.descr))
        report$add.printed(format(descdist(xl,graph=F)),
                           caption=sprintf("Additional descriptive parameters of %s for %s",
                                           feat.name,lev.descr))
        for(discrete in c(F,T)) {
          report$add(descdist(xl,boot=1000,discrete=discrete),
                     caption=sprintf("Skewness-kurtosis of %s for %s",
                                     feat.name,lev.descr))
        }
      }
    }
    #report$add(roc(g,x,plot=T,smooth=F,ci=F,print.thres=T,grid=c(0.1,0.2)),
    #           caption=sprintf("ROC of %s for predicting %s",feat.name,group.attr))
    tryCatchAndWarn({ 
      
      pwr.res = show.partial.auc.roc(g,x,predictor.descr=feat.name)
      
      report$add.table(as.data.frame(pwr.res$power),
                       caption=sprintf("Power analysis for simple predictor based on %s level",feat.name)
      )
    })
    report$pop.section()
  }
}

## power analysis of biomarker verification study
verification.power <- function(m_a,
                               group.attr,
                               id.markers=NULL) {
  report.section = report$add.header("Power analysis of a verification study",
                                     section.action="push", sub=F)
  
  if(!is.null(id.markers)) {
    m_a = subset.m_a(m_a,select.count=id.markers)
  }
  report$add.header("Empirical distributions of individual features")
  report$push.section(report.section)  
  for(norm.method in c("ident")) {
    m_a.norm = norm.count.report(m_a,norm.count.task=list(method=norm.method))
    counts.distro.report(m_a.norm,group.attr=group.attr,
                         descr=sprintf("transformation method %s",norm.method))
  }
  report$pop.section()
}

test.run.verification.power <- function() {
  report <<- PandocAT$new()
  verification.power(m_a,group.attr="Group",id.markers= c("P17050","Q9BTY2","P04066"))
  report$save()
}

power.pieper.t1d <- function(
  n = 50,
  alpha.sim = 0.05,
  alpha.orig = 0.05,
  min.mean = 100,
  mult.adj = "fdrtool",
  data.file="aim3/September 28 analysis_T1D.txt",
  R = 4000,
  feat.sig.names=NULL,
  targeted=F
) {
  if(T) {
    taxa.meta = read.pieper.t1d(file.name=data.file)
    aggr_var = "SubjectID"
    taxa.meta = aggregate.by.meta.data(taxa.meta$data,
                                       aggr_var,
                                       taxa.meta$attr.names)
    
    
  }
  taxa.meta.attr.names = taxa.meta$attr.names    
  dim.data.orig = dim(count_matr_from_df(taxa.meta$data,taxa.meta.attr.names))
  taxa.meta.data.raw = count.filter(taxa.meta$data,col_ignore=taxa.meta.attr.names,
                                    min_max_frac=0,min_max=0,min_mean=min.mean,min_row_sum=0,other_cnt=NULL)
  print(dim(taxa.meta.data.raw))
  
  taxa.meta.data = norm.meta.data(taxa.meta.data.raw,col_ignore=taxa.meta.attr.names,norm.func=ihs)
  
  
  #feature.names = get.feature.names(taxa.meta.data,taxa.meta.attr.names)
  #print(feature.names)
  #pvals = wilcox.test.multi(data=taxa.meta.data,resp.vars=feature.names,group.var="group",subset=NULL)
  m = count_matr_from_df(taxa.meta.data,taxa.meta.attr.names)
  
  if(!is.null(feat.sig.names)) {
    ind.sig = which(colnames(m) %in% feat.sig.names)
  }
  else {
    #ind.sig = which(pvals.adj<=alpha.orig)
    ind.sig = which(eff.raw <= 0.66 | eff.raw >= 1.5)
  }
  
  if(targeted) {
    m = m[,ind.sig]
    ind.sig=1:length(ind.sig)
  }
  
  print(length(ind.sig))
  #return(0)  
  
  dim.data.filt = dim(m)
  #mtp.res = MTP(X=t(m),Y=taxa.meta.data$group,get.adjp=F)
  
  group = taxa.meta.data$group
  m.raw = count_matr_from_df(taxa.meta.data.raw,taxa.meta.attr.names)
  eff.raw = group.mean.ratio(m.raw,group)
  
  pvals = wilcox.test.multi.fast(m,group)  
  ##make.global(pvals)
  
  if(mult.adj=="fdrtool") {
    pvals.adj = fdrtool(pvals,statistic="pvalue",plot=F,verbose=F)$qval
  }
  else {
    pvals.adj = p.adjust(pvals,method=mult.adj)
  }
  
  ##make.global(pvals.adj)
  
  group.mean = count_matr_from_df(count.summary(m,mean,group),c(".group"))
  group.sd = count_matr_from_df(count.summary(m,sd,group),c(".group"))
  group.mean.sig = group.mean[,ind.sig]
  #print(group.mean.sig)
  group.sd.sig = group.sd[,ind.sig]
  #print(group.sd.sig)
  #group.n = c(table(group))[rownames(group.mean.sig)]
  group.n = count(group)
  cohens.d = foreach(i.var=seq(ncol(group.mean.sig)),.combine=c) %do% {
    cohens.d.from.mom(mean.gr=group.mean.sig[,i.var],var.gr=group.sd.sig[,i.var]**2,n.gr=group.n)
  }
  #print(cohens.d)
  pvals.boot = boot(data=m,
                    statistic=booted.wilcox.test.multi.fast,
                    R=R,
                    n=n,
                    strata=group,
                    group=group,
                    test.method=wilcox.test.multi.fast)$t
  ## some numerical instability creates pvals that are slightly larger than 1
  pvals.boot[pvals.boot>1] = 1  
  
  if(mult.adj=="fdrtool") {
    pvals.boot.adj = plyr::aaply(pvals.boot,1,function(x) fdrtool(x,statistic="pvalue",plot=F,verbose=F)$qval)
  }
  else {
    pvals.boot.adj = plyr::aaply(pvals.boot,1,p.adjust,method=mult.adj)
  }
  
  power.sig = colMeans(pvals.boot.adj[,ind.sig] <= alpha.sim)
  names(power.sig) = names(pvals.adj[ind.sig])
  
  #print(power.sig)
  #print(mean(power.sig))
  return(list(cohens.d=cohens.d,
              mean.cohens.d=mean(cohens.d),
              group.mean.sig=group.mean.sig,
              group.sd.sig=group.sd.sig,
              group.n=group.n,
              power.sig=power.sig,
              mean.power.sig=mean(power.sig),
              pvals.adj.sig=pvals.adj[ind.sig],
              dim.data.orig=dim.data.orig,
              dim.data.filt=dim.data.filt,
              alpha.orig=alpha.orig,
              alpha.sim=alpha.sim))
}

pediatric.cancer.2013.aim2 <- function() {
  doc = "(http://www.ncbi.nlm.nih.gov/pubmed/17414136)
  
  Fecal specimens were collected from 148 consecutive pediatric patients 
  (79 with Crohn disease, 62 with ulcerative colitis, and 7 with irritable 
  bowel syndrome) and 22 healthy control individuals.
  
  Lactoferrin levels were significantly higher in patients with ulcerative 
  colitis (1880 +/- 565 microg/mL) (mean +/- SE) or Crohn disease 
  (1701 +/- 382 microg/mL) than in healthy control individuals under 21 
  years of age (1.17 +/- 0.47 microg/mL, P < 0.001). 
  "
  
  mean.gr = c(1880,1.17)
  se.gr = c(565,0.47)
  n.gr = c(79,22)
  
  var.gr = se.gr**2*n.gr
  sd.gr = var.gr**0.5
  
  mom.f.gr = foreach(i.gr=1:2,.combine=rbind) %do% {
    mom.f(var.x=var.gr[i.gr],mean.x=mean.gr[i.gr],f=ihs,f.d1=ihs.d1,f.d2=ihs.d2)
  }
  cohens.d = cohens.d.from.mom(mean.gr=mom.f.gr[,"mean.f"],var.gr=mom.f.gr[,"var.f"],n.gr=n.gr)
  return(list(mom.f.gr=mom.f.gr,mean.gr=mean.gr,
              se.gr=se.gr,var.gr=var.gr,sd.gr=sd.gr,
              cohens.d=cohens.d))
}

power.pediatric.cancer.2013<-function() {
  power.res = power.pieper.t1d(
    n = 50,
    data.file="aim3/September 28 analysis_T1D.txt"
  )
  print("Power results Aim 3:")
  print(power.res)
  mom = pediatric.cancer.2013.aim2()
  print("Power results Aim 2:")
  print(paste("Cohen's d Lactoferrin:",mom$cohens.d))
  print(paste("Mean Cohen's d for significant proteins from Aim 3:",power.res$mean.cohens.d))
}

power.pieper.prostate.cancer.2014<-function() {
  power.res = power.pieper.t1d(
    n = 170,
    ## this is under the proposal directory
    data.file = "data/T1D_proteome/Original Collapesed APEX (All information).AT.tsv",
    min.mean = 1200,
    alpha.sim = 0.05,
    alpha.orig = 0.05,
    R = 400
  )
  print("Power results:")
  print(power.res)
}

power.proteom.bladder.cancer.probiotic <-function() {
  power.res = power.pieper.t1d(
    n = 170,
    ## this is under the proposal directory
    data.file = "data/T1D_proteome/Original Collapesed APEX (All information).AT.tsv",
    min.mean = 1200,
    alpha.sim = 0.05,
    alpha.orig = 0.05,
    R = 400
  )
  print("Power results:")
  print(power.res)
}

power.madupu.kidney_diabetes<-function() {
  prot.ids = as.data.frame(
    matrix(
      c(
        "MAN2B1","O00754",
        "GP5","P40197",
        "FUCA2","Q9BTY2",
        "ATP1A1","P05023",
        "CDH5","P33151",
        "ACE2","Q9BYF1"
      ),
      ncol=2,
      byrow=T
    )
  )
  names(prot.ids) = c("name","uniref100")
  
  power.res = power.pieper.t1d(
    n = 300,
    ## this is under the proposal directory
    data.file = "data/T1D_proteome/Original Collapesed APEX (All information).AT.tsv",
    min.mean = 0,
    alpha.sim = 0.05,
    alpha.orig = 0.05,
    R = 400,
    mult.adj="bonferroni",
    feat.sig.names=prot.ids$uniref100,
    targeted=T
  )
  print("Power results:")
  print(power.res)  
}


wilcox.power <- function(m_a,
                         group.attr="Group",
                         alpha=0.05,
                         power.target=0.8,
                         n = 300,
                         R = 100,
                         mult.adj="BH",
                         id.markers)
{
  
  orig.res = genesel.stability(m_a,
                               group.attr=group.attr,
                               block.attr=NULL,
                               type="unpaired",
                               replicates=0,
                               samp.fold.ratio=0.5,
                               maxrank=20,
                               comp.log.fold.change=T,
                               ret.data.contrasts=T
  )
  
  m = m_a$count
  group = m_a$attr[,group.attr]
  
  ind.ob = seq(nrow(m))
  plus = function(x,y) (x + y)
  fin = function(x) (x/R)
  power = foreach(i.rep = seq(R),
                  .combine=plus,
                  .final=fin,
                  .packages=c("GeneSelector","fdrtool"),
                  .export=c("wilcox.test.multi.fast","RankingWilcoxonAT")) %dopar% {
                    ##replace=TRUE to select more samples than in original dataset
                    ##TODO: apply strata to keep group count ratio
                    ind.n = sample(ind.ob, n, replace = TRUE, prob = NULL)
                    d <- m[ind.n,,drop=F]
                    g <- group[ind.n]
                    #print(summary(g))
                    pvals = wilcox.test.multi.fast(d,g,pvalues.impl="wilcox.exact")    
                    if(mult.adj=="fdrtool") {
                      pvals.adj = fdrtool(pvals,statistic="pvalue",plot=F,verbose=F)$qval
                    }
                    else {
                      pvals.adj = p.adjust(pvals,method=mult.adj)
                    }
                    as.numeric(pvals.adj <= alpha)
                  }
  #colnames(pvals.boot.adj) = colnames(m)
  #power = colMeans(pvals.boot.adj <= alpha)
  power = data.frame(power=power)
  rownames(power) = colnames(m)
  #power$name = rownames(power)
  res = merge(orig.res$stab_feat,power,by="row.names")
  stopifnot(nrow(res)==nrow(power))
  rownames(res) = res$name
  res$selection = F
  res$selection[rownames(res) %in% id.markers] = T
  report$add.table(res)
  
  require(splines)
  pl = ggplot(res, aes(abs(r.eff.size),power)) + geom_point(aes(color=selection))
  if(max(res$power)>=0.95) {
    pl = pl + geom_smooth(method = "gam", se=T, color="black",family="betar")
  }
  else {
    pl = pl + geom_smooth(method="lm",formula = y ~ ns(x,df=100), se=T, color="black")
  }
  pl = pl + 
    geom_hline(y=power.target,color="blue") +
    xlab("Abs value of Cohen's r")
  #geom_smooth(method = "gam", se=T, color="black",family="betar")
  report$add(pl,
             caption=sprintf("Power of Wilcoxon rank-sum test as a function of absolute 
                             value of standardized r effect size at %s significance cutoff
                             and sample size %s (total in both groups). 
                             Each point represents one feature. 
                             Externally selected features are colored.
                             Blue line represents target power value of %s", alpha,n,power.target)
  )
  return (res)
}
