#!/usr/bin/env Rscript

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

SRC_DIR = dirname(this_script_file())

source(file.path(SRC_DIR,"dada_pipeline_methods.R"),local=T)

dada_main()
