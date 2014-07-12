
options(mc.cores=4)
options(boot.ncpus=4)
options(boot.parallel="snow")
#cl<-makeCluster(getOption("mc.cores", 2L)) #number of CPU cores
#registerDoSNOW(cl)  

#power.nistal()
#proc.choc()
#meta = load.meta.t1d("annotation20130819.csv")
#mgrast.dir = "../BATCH_01_02_META/BATCH_01-02_METAGENOMICS_MGRAST"
#mgr = read.mgrast.summary(paste(mgrast.dir,"cog.tsv",sep="/"),file_name.id.map=paste(mgrast.dir,"mgrast_to_samp_id.tsv",sep="/"))
#mgr.cnt = mgrast.to.abund.df(mgr,"level.2")
#taxa.meta = read.t1d(3)
#taxa.meta.data = taxa.meta$data
#taxa.meta.attr.names = taxa.meta$attr.names

#sink("analysis.log",split=T)

# panderOptions("table.style","rmarkdown")
# panderOptions("table.split.table",180)
# panderOptions("table.alignment.default","left")
# panderOptions("evals.messages",F)
# 
# evalsOptions("cache",F)
# evalsOptions("cache.mode","environment")
# #evalsOptions("output",c("result"))
# evalsOptions("output",c("all"))
# evalsOptions("graph.unify",F)
# evalsOptions("res",100)
# evalsOptions("hi.res",T)

MGSAT_SRC = "~/work/mgsat"

source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

set_trace_options(try.debug=F)

setwd("~/work/CHOC/stage_3")


report <- PandocAT$new("atovtchi@jcvi.org","CHOC pre-B ALL pilot 16S Analysis")

proc.choc()

#print(report$body)

report$save("report")
Pandoc.convert("report.md",format="html",open=F)
Pandoc.convert("report.md",format="docx",open=F)
#Pandoc.convert("report.md",format="pdf",open=F)



#print(plot.stability.selection.c060.at(s),rank="mean")
#spca.res = spca(count,K=6,type="predictor",sparse="varnum",trace=T,para=c(7,4,4,1,1,1))
##This has selected a single variable (Gordonibacter) from 16S level 3 with normalization "ident".
##Note that by default this function does not scale the predictor variables.
#hc.res <- get.biom(X = m_a$count, Y = m_a$attr$T1D, fmethod = c("studentt", "pls", "vip"), type = "HC")
#taxa.meta$data[taxa.meta$data$Batch==3,c("Gordonibacter_0.1.1.1.3.1.6","T1D")]
#taxa.meta$data[taxa.meta$data$Batch==3,c("Akkermansia_0.1.8.2.1.1.1","T1D")]

#taxa.meta = read.mr_oralc(3)
#proc.mr_oralc()
#power.pediatric.cancer.2013()
#stopCluster(cl)

#sink(NULL)
