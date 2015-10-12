
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
  "data.table",
  "dplyr",
  ##ROC curve power analysis for biomarker verification
  "pROC",
  ##interpolation and numerical derivatives (including for ROC analysis)
  "pspline",
  ##permutation test with multiple testing correction for dependent tests
  ##to be used for testing series of Hill numbers (use our modified function
  ##mcp.wy, but it depends on one other method from simboot for now)
  "simboot",
  ## automatic selection of the number of clusters for pam
  "fpc",
  ## descriptive statistics and multinomial confidence intervals
  "DescTools",
  ## to install development versions of packages from GitHub
  "devtools",
  ## Matrix of ggplots and other extensions
  ##GGally
  ## binnedplot to diagnose logistic regression models
  "arm",
  ## mixed() to get p-values for lmer
  "afex",
  ## htmlwidgets-based packages
  "htmlwidgets",
  "d3heatmap",
  "threejs"
)

vanilla_packages_github = c(
  "zdk123/SpiecEasi",
  "juba/scatterD3",
  "bokeh/rbokeh",
  "bwlewis/doRedis"
)

bio_packages = c(
  "multtest",
  "GeneSelector",
  "RColorBrewer",
  "Heatplus",
  "DESeq2",
  "ComplexHeatmap",
  "phyloseq",
  ##third-party R implementation of Holmes 2012 model-based 
  ##and classification clustering algorithm.
  ##Needs GSL installed on the system in order to build
  "DirichletMultinomial"
)

install_required_packages <- function() {
  install.packages(vanilla_packages)
  source("http://bioconductor.org/biocLite.R")
  biocLite(bio_packages)
  library(devtools) ## for install_github
  for(pkg in vanilla_packages_github) {
    install_github(pkg)
  }
}

packages = c(vanilla_packages,bio_packages)

load_required_packages <- function() {
  for (package in packages) {
    suppressMessages(library(package,character.only=T))
  }  
}
