#!/usr/bin/env Rscript

load_sample_manifest <- function(sample_sheet_file,
                                 study_design_file=NULL,
                                 sample_sheet_root_path=NULL,
                                 normalize_paths=T) {
  samp_manifest = data.table::fread(sample_sheet_file,header = T,
                                    colClasses=c(SampleID="character"))
  if(!is.null(sample_sheet_root_path)) {
    samp_manifest[,file1 := file.path(sample_sheet_root_path,file1)]
    samp_manifest[,file2 := file.path(sample_sheet_root_path,file2)]
  }
  if( normalize_paths ) {
    samp_manifest[,file1 := normalizePath(file1)]
    samp_manifest[,file2 := normalizePath(file2)]
  }
  meta = samp_manifest
  if(!is.null(study_design_file)) {
    study_des = data.table::fread(study_design_file,header = T,
                                    colClasses=c(SampleID="character"))
    meta = meta[study_des,on="SampleID",nomatch=0]
    if(!(nrow(meta)==nrow(study_des) && all(meta$SampleID==study_des$SampleID))) {
        miss_samp_id = setdiff(study_des$SampleID,meta$SampleID)
        if(length(miss_samp_id)>0) {
            warning(sprintf("The following SampleID values from study design file 
                 are missing from sample sheet file: %s",
                 paste(miss_samp_id,collapse=",")))
        }
        else {
            stop("Detected one-to-many match between SampleID values in
                 study design and sample sheet files")
        }
    }
  }
  return (meta)
}

validate_sample_manifest <- function(meta) {
  ##TODO: consider NAs
  missing_files = meta[!(file.exists(file1) & file.exists(file2))]
  return (list(missing_files=missing_files))
}

report_seq_quality <- function(meta,descr="",threads=0,n_reads=100000) {
  #report.section = report$get.section()
  #report.section = report$add.header("Quality profiles",section.action = "push", report.section=report.section, sub = T)
  report.section = report$add.header(sprintf("Quality profiles. %s",descr),section.action = "push", sub = T)
  plots = foreach(i_samp = seq(nrow(meta)),
          .packages = c("data.table","dada2","ggplot2")) %dopar% {
    rec = meta[i_samp]
    res = list(file1=rec$file1,file2=rec$file2,SampleID=rec$SampleID)
    res$plot1 = dada2::plotQualityProfile(rec$file1,n=n_reads)
    if(!is.na(rec$file2)) {
      res$plot2 = dada2::plotQualityProfile(rec$file2,n=n_reads)
    }
    res
  }
  for(rec in plots) {
    report$add.header(sprintf("Quality profile for sample %s",rec$SampleID),report.section=report.section, section.action = "push")
    report$add(rec$plot1,
               caption = sprintf("Quality profile for sequence file %s", rec$file1))
    if(!is.na(rec$file2)) {
      report$add(rec$plot2,
                 caption = sprintf("Quality profile for sequence file %s", rec$file2))
    }
    report$pop.section()
  }
  report$pop.section()
}

#maxN=0, maxEE=2, truncQ=2, trimLeft=c(10, 15), truncLen=c(240,160), 
#for 300x2 kit, truncLen=c(240,160) results in very low number of output sequences
#TODO: now DADA2 author recommends trimLeft=0 "assuming the primers have been removed", 
#but says trimLeft=10 is still OK. Explore.
#https://github.com/benjjneb/dada2/issues/226#issuecomment-298917705
dada_trim_filter <- function(meta,seq_dir_filtered=NULL,
                             use_existing_files=F,
                             maxN=0, maxEE=2, truncQ=2, trimLeft=c(10, 10), 
                             truncLen=c(240,160), 
                             rm.phix=TRUE, compress=TRUE, verbose=TRUE) {
  library(foreach)
  if(is.null(seq_dir_filtered)) {
    seq_dir_filtered = "dada_filtered_seq"
  }
  dir.create(seq_dir_filtered,showWarnings = F,recursive = T)
  meta_filt = data.table::copy(meta)
  ##TODO: make it so that assert below is not necessary
  meta_filt[,`:=`(
    file1 = file.path(seq_dir_filtered,sprintf("%s_%s",SampleID,basename(file1))),
    file2 = file.path(seq_dir_filtered,sprintf("%s_%s",SampleID,basename(file2)))
  )]
  stopifnot(anyDuplicated(meta_filt$file1)==0)
  stopifnot(anyDuplicated(meta_filt$file2)==0)
  foreach(i_rec = seq(nrow(meta)),
          .packages = "data.table") %dopar% {
            rec = meta[i_rec]
            rec_filt = meta_filt[i_rec]
            files = c(rec$file1,rec$file2)
            filt_files = c(rec_filt$file1,rec_filt$file2)
            if(!(all(file.exists(filt_files)) && use_existing_files)){
              dada2::fastqPairedFilter(files, filt_files, 
                                       maxN=maxN, maxEE=maxEE, truncQ=truncQ, 
                                       trimLeft=trimLeft, 
                                       truncLen=truncLen, rm.phix=rm.phix,
                                       compress=compress, verbose=verbose)
            }
          }
  return (meta_filt)
}

dada_dereplicate <- function(meta) {
  paths = list(file1=meta$file1,file2=meta$file2)
  res = foreach(path=paths) %do% {
    derep_dat = foreach(fpath=path) %dopar% { 
      dada2::derepFastq(fpath, verbose=F)
    }
    names(derep_dat) <- meta$SampleID
    derep_dat
  }
  names(res) = names(paths)
  return (res)
}

assert_same_length_lists_ <- function(x){
  n = length(x[[1]])
  stopifnot(all(sapply(x,length)==n))
  return (n)
}

dada_estimate_error_rate <- function(dat,n_subset=20,threads=0) {
  if(threads<1) threads = TRUE
  n_samp = assert_same_length_lists_(dat)
  ind_sel = sample(n_samp,min(n_samp,n_subset))
  res = list()
  for(i_dat in seq_along(dat)) {
    res[[i_dat]] <- dada2::dada(dat[[i_dat]][ind_sel], err=NULL, 
                                selfConsist = TRUE, multithread=threads)[[1]]$err_out
  }
  names(res) = names(dat)
  return (res)
}

dada_infer <- function(dat,err_dat,pool=F,threads=TRUE) {
  if(threads<1) threads = TRUE
  res = list()
  for(i_dat in seq_along(dat)) {
    res[[i_dat]] <- dada2::dada(dat[[i_dat]], err=err_dat[[i_dat]], 
                                pool=pool, selfConsist=FALSE, multithread=threads)
  }
  names(res) = names(dat)
  return (res)
}

format_timer <- function(msg,ptm) {
    time_used = proc.time() - ptm
    time_used = paste(names(time_used),time_used,sep="=",collapse = ", ")
    sprintf("%s. Time used so far: %s",msg,time_used)
}

report_timer <- function(msg,ptm) {
    msg = format_timer(msg,ptm)
    report$add.descr(msg)
    message(msg)
}

dada_tab_to_mothur <- function(tab_file,otu_shared_file,cons_taxonomy_file,label="0.03") {
  library(data.table)
  tab = data.table::fread(tab_file,header = T)
  ind_seq = which(colnames(tab)=="sequence")
  taxa = tab[,(ind_seq+1):ncol(tab),with=F]
  otu = t(tab[,1:(ind_seq-1),with=F])
  colnames(otu) = taxa$OTU
  otu_sums = colSums(otu)
  otu_sums = data.table(OTU=names(otu_sums),Size=otu_sums)
  otu = data.table(label=label,Group=rownames(otu),numOtus=ncol(otu),otu)
  write.table(otu,otu_shared_file,sep="\t",row.names = F,col.names = T)
  taxa_lin = taxa[,!"OTU",with=F]
  taxa_lin[which(is.na(taxa_lin),arr.ind = T)] = "unclassified"
  taxa_lin = apply(taxa_lin,1,function(row) paste0(paste(sprintf("%s(100)",row),collapse=";"),";"))
  taxa = data.table(taxa[,.(OTU)],Taxonomy=taxa_lin)
  taxa = otu_sums[taxa,on="OTU"]
  write.table(taxa,cons_taxonomy_file,sep="\t",row.names = F,col.names = T)
}

dada_tab_concatenate <- function(tab_files,tab_file,sample_map_file,otu_map_file) {
  library(data.table)
  if(is.null(names(tab_files))) {
    names(tab_files) = sprintf("D%02i",seq_along(tab_files))
  }
  taxa = list()
  count = list()
  for(pref_t in names(tab_files)) {
    tab = data.table::fread(tab_files[[pref_t]],header = T)
    ind_seq = which(colnames(tab)=="sequence")
    tx = tab[,ind_seq:ncol(tab),with=F]
    tx[,pref:=pref_t]
    taxa = rbind(taxa,tx,use.names=T,fill=T)
    cnt = tab[,c(1:(ind_seq-1),ind_seq+1),with=F]
    cnt = list(cnt)
    names(cnt) = pref_t
    count = c(count,cnt)
  }
  uni_otu = taxa[,.(n_sets=.N,sets=paste0(pref,collapse = "")),by=sequence]
  uni_otu[,OTU_uni:=paste0("OTU",index_as_left_padded_str(.I))]
  taxa = taxa[uni_otu,on="sequence"]
  count_tall = NULL
  uni_sample_id = NULL
  for(pref_t in names(count)) {
    cnt = count[[pref_t]][taxa[pref==pref_t,.(OTU,OTU_uni)],on="OTU"]
    cnt[,OTU:=OTU_uni]
    cnt[,OTU_uni:=NULL]
    sample_ids = data.table(pref=pref_t,SampleID=colnames(cnt))[SampleID!="OTU"]
    sample_ids[,SampleID_uni:=paste0(pref_t,SampleID)]
    setnames(cnt,sample_ids$SampleID,sample_ids$SampleID_uni)
    cnt_tall = melt(cnt,id.vars = "OTU",value.name = "count",variable.name = "SampleID")[count!=0]
    count_tall = rbind(count_tall,cnt_tall)
    uni_sample_id = rbind(uni_sample_id,sample_ids)
  }
  extra_fields = c("pref", "n_sets", "sets", "OTU_uni")
  uni_taxa = taxa[,c("OTU",extra_fields),with=F]
  taxa[,OTU:=OTU_uni]
  taxa = taxa[,-extra_fields,with=F][!duplicated(OTU)]
  count = dcast(count_tall,OTU~SampleID,value.var = "count",fill = 0)
  ##OTU is the first field after the above - relying on that downstream
  stopifnot(which(colnames(count)=="OTU")==1)
  ##sanity check
  stopifnot(!any(colSums(count[,-"OTU",with=F])==0))
  count = count[taxa,on="OTU"]
  ind_seq = which(colnames(count)=="sequence")
  setcolorder(count,c(2:ind_seq,1,(ind_seq+1):ncol(count)))
  ind_seq = which(colnames(count)=="sequence")
  ##order rows by total count just for convenience
  count$.order_col=rowSums(count[,1:(ind_seq-1),with=F])
  setorder(count,-.order_col)
  count[,.order_col:=NULL]
  ##final sanity checks
  stopifnot(all(sort(count$OTU)==sort(uni_otu$OTU_uni)))
  write.csv(count,tab_file,row.names = F)
  write.table(uni_sample_id,sample_map_file,sep="\t",row.names = F,col.names = T)
  write.table(uni_taxa,otu_map_file,sep="\t",row.names = F,col.names = T)
}

run_dada_pipeline <- function(load_sample_manifest_task,
                              do_report_seq_quality=F,
                              pool=F,
                              rdp_train_file,
                              threads=0) {
  library(data.table)
  ptm = proc.time()
  report_timer("Starting DADA pipeline",ptm)
  report.section = report$get.section()
  multithread = TRUE
  if(threads>0) {
      multithread = threads
  }
  meta = do.call(load_sample_manifest,load_sample_manifest_task)
  val_res = validate_sample_manifest(meta)
  file.remove("missing_files.csv")
  if(nrow(val_res$missing_files)>0) {
    write.csv(val_res$missing_files,"missing_files.csv",row.names = F)
    stop("Missing files in the manifest, see missing_files.csv")
    #meta = meta[!(SampleID %in% val_res$missing_files$SampleID)]
  }
  report_timer(sprintf("Loaded manifests, processing %s samples",nrow(meta)),ptm)
  if( do_report_seq_quality ) {
    message("QA of input sequences")
    report_seq_quality(meta=meta,descr = "Input sequences")
  }
  
  report_timer("Trimming and filtering input sequences",ptm)
  ## this discards the information about the original sample paths
  meta = dada_trim_filter(meta=meta,use_existing_files = F)
  #meta = dada_trim_filter(meta=meta,use_existing_files = F,truncLen=150,trimLeft=0)
  
  if(do_report_seq_quality) {
    report_timer("QA of trimmed sequences",ptm)
    report_seq_quality(meta=meta,descr = "Trimmed and filtered sequences")
  }
  
  report_timer("Dereplicating",ptm)
  derep_dat = dada_dereplicate(meta = meta)
  
  report_timer("Estimating error rate",ptm)
  err_dat = dada_estimate_error_rate(derep_dat,n_subset = 20,threads=threads)
  
  report_timer("Inferring",ptm)
  dada_dat = dada_infer(derep_dat,err_dat = err_dat,pool=pool,threads=threads)
  
  report_timer("Merging pairs",ptm)
  mergers = dada2::mergePairs(dada_dat[[1]], derep_dat[[1]], 
                              dada_dat[[2]], derep_dat[[2]])
  
  report_timer("Making sequence table",ptm)
  seqtab = dada2::makeSequenceTable(mergers)
  
  report_timer("Removing chimeras",ptm)
  seqtab_nochim = dada2::removeBimeraDenovo(seqtab, verbose=TRUE,
                                            multithread=multithread)
  
  report_timer("Assigning taxonomy",ptm)
  taxa = dada2::assignTaxonomy(seqtab_nochim, rdp_train_file)
  colnames(taxa) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus")
  
  report_timer("Exporting results",ptm)
  seqtab_nochim_t = t(seqtab_nochim)
  tab <- cbind(data.table(seqtab_nochim_t),
               sequence=rownames(seqtab_nochim_t),
               OTU=sprintf("OTU%04i",seq(nrow(seqtab_nochim_t))))
  stopifnot(all(tab$sequence==rownames(taxa)))
  tab <- cbind(tab,taxa)
  write.csv(tab,"tab.csv",row.names = F)
  #grp_abund = tab[S1+S2+S3>0,.(S1=sum(S1),S2=sum(S2),S3=sum(S3)),by=.(Family,Genus)][order(Family,Genus)]
  #grp_abund[,`:=`(S1_r=S1/sum(S1),S2_r=S2/sum(S2),S3_r=S3/sum(S3))]

  dada_tab_to_mothur(tab_file="tab.csv",otu_shared_file="otu.shared",cons_taxonomy_file="cons.taxonomy")
  
}

dada_install_packages <- function() {
    install.packages("optparse")
    ##This installs the current master from GitHub,
    ##and follows the instructions from
    ##http://benjjneb.github.io/dada2/dada-installation.html
    source("http://bioconductor.org/biocLite.R")
    biocLite(suppressUpdates = FALSE)
    biocLite("ShortRead", suppressUpdates = FALSE)
    library("devtools")
    ##If you are on RHEL 6.5 and using GCC from Conda,
    ##you will need to also install binutils from conda-forge,
    ##otherwise the older assembler from RHEL will error out
    ##on unfamiliar assembler instructions generated by the newer
    ##GCC.
    devtools::install_github("benjjneb/dada2")
}

dada_load_packages <- function() {
  library(dada2); packageVersion("dada2")
  library(ShortRead); packageVersion("ShortRead")
  library(ggplot2); packageVersion("ggplot2")
  library(data.table); packageVersion("data.table")
}


this_script_file <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}

dada_pipeline_methods_init <- function() {
  
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
  library(methods) #RScript does not load methods
  ## loads MGSAT code
  source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"))
  source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"))
}

dada_main <- function() {

dada_pipeline_methods_init()

conf = jsonlite::fromJSON(file(file.path(dirname(this_script_file()),"dada.json"),"r"))

## seems to be a bug in optparse - if using add_option, options with hyphens
## or underscores cause errors in parse_args; work fine with make_option.
option_list = list(
optparse::make_option("--study_design_file", help="Study design table"),
optparse::make_option("--sample_sheet_file", help="Sample sheet table"),
optparse::make_option("--sample_sheet_root_path", 
                     help="Path to prepend to file names in sample sheet"),
optparse::make_option("--database_path", default=conf$database,
                     help="Path to reference database directory"),
optparse::make_option("--report_seq_quality", default=FALSE,
                     action="store_true",
                     help="Report QC metrics for each sample"),
optparse::make_option("--pool", default=FALSE,
                      action="store_true",
                      help="Pool samples for error correction. Note: this has major effect on speed and, in some cases, results."),
optparse::make_option("--threads", type="integer", default=0,
                     help="Max number of threads to use [0=total cores on machine]")
)

parser = optparse::OptionParser(option_list=option_list)
args = optparse::parse_args(parser,positional_arguments = TRUE)

ops = args$options

## number of cores to use on multicore machines

if (ops$threads > 0) {
    num.cores = ops$threads
}
else {
    ## By default, use all cores
    num.cores = parallel::detectCores()
}

options(mc.cores=num.cores)
options(boot.ncpus=num.cores)
## parallel backend
options(boot.parallel="parallel")

library("BiocParallel")
#register(SnowParam(num.cores))

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

evalsOptions("graph.output","svg")

report <<- PandocAT$new(author="tovchigrechkoa@medimmune.com",
                       title="Using DADA2 to annotate 16S sequences",
                       incremental.save=F)


cl = start.cluster.project()

rdp_train_file = file.path(ops$database_path,"rdp_train_set_14.fa.gz")

load_sample_manifest_task = list(sample_sheet_file = ops$sample_sheet_file, 
                                 study_design_file = ops$study_design_file,
                                 sample_sheet_root_path = ops$sample_sheet_root_path)

run_dada_pipeline(load_sample_manifest_task=load_sample_manifest_task,
                  do_report_seq_quality=ops$report_seq_quality,
                  rdp_train_file=rdp_train_file,
                  threads=ops$threads)

stop.cluster.project(cl)

report$save()

}

