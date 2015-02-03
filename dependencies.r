
vanilla_packages = c(        
  #for quantcut; mask as little as possible (masks permute, 
  "gtools", #pos = "package:base", 
  "reshape2",
  "plyr",
  # for str_split with max splits
  "stringr", 
  "ggplot2",
  "vegan", 
  #"BiodiversityR", 
  "LiblineaR", 
  ## glmnet is exported to snow cluster by c060, but it
  ## forgets to load it first, so we do it here
  "glmnet",
  "stabs",
  #"c060", 
  "geoR", 
  "foreach", 
  "iterators", 
  "doSNOW", 
  "fdrtool", 
  #"elasticnet", 
  #"BioMark", 
  "HMP", 
  "BatchJobs", 
  "boot",
  "vegan",
  #"knitr",
  "lattice",
  "date",
  "timeDate",
  #for llist
  "Hmisc",
  "kernlab",
  "lme4",
  "pander",
  ##fitting to parametric distributions and goodness of fit tests
  "fitdistrplus",
  ##SQL-like aggregate tables, used for metadata summaries
  "doBy",
  ##correct treatment of zeros in signed-rank Wilcoxon test
  ##with an interface compatible with base::wilcox.test while
  ##being relatively fast
  "exactRankTests",
  ##ROC curve power analysis for biomarker verification
  "pROC"
)

bio_packages = c(
  "multtest",
  "GeneSelector",
  "RColorBrewer",
  "Heatplus",
  "DESeq2"
)

install_required_packages <- function() {
  install.packages(vanilla_packages)
  source("http://bioconductor.org/biocLite.R")
  biocLite(bio_packages)
}

packages = c(vanilla_packages,bio_packages)

load_required_packages <- function() {
  for (package in packages) {
    suppressMessages(library(package,character.only=T))
  }  
}
