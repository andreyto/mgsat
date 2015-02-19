
## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"
source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)
## loads dependency packages (which already must be installed)
load_required_packages()

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

library(fitdistrplus)

m_a = read.table.m_a("prot.example")

## leave with try.debug=F for production runs
set_trace_options(try.debug=T)

## set incremental.save=T only for debugging or demonstration runs - it forces 
## report generation after adding every header section, thus slowing down
## a long run. But then incremental.save=T, you can open HTML report file in
## a Web browser and refresh it periodically to see it grow.
report <- PandocAT$new(author="noone@email.com",
                       title="Example of the power analysis for protein biomarker verification study",
                       incremental.save=F)

id.markers = c("M1","M2","M3","M4","M5","M6","M7","M8","M9",
               "M10","M11","M12","M13","M14","M15","M16",
               "M17","M18","M19")

verification.power(m_a=m_a,
                   group.attr="Group",
                   id.markers=id.markers)

report$save()

## Tasks
## 1. Use show.distr.group function to plot empirical 
##    density plots for three proteins from the list above
