
## location of MGSAT code
MGSAT_SRC = "~/work/mgsat"
source(paste(MGSAT_SRC,"dependencies.r",sep="/"),local=T)
## loads dependency packages (which already must be installed)
load_required_packages()

## loads MGSAT code
source(paste(MGSAT_SRC,"report_pandoc.r",sep="/"),local=T)
source(paste(MGSAT_SRC,"power_and_tests.r",sep="/"),local=T)

library(fitdistrplus)

pool.sd <- function(n.1,sd.1,n.2,sd.2) {
  sqrt(((n.1-1)*sd.1**2+(n.2-1)*sd.2**2)/(n.1+n.2))
}

diversity.power.example <- function() {
  m_a = read.table.m_a("divrich.counts")
  m_a = norm.count.m_a(m_a,method="norm.ihs")
  
  x = m_a$count[,"N1"]
  group = m_a$attr[,"Sample.type"]
  print(show.distr.group(x = x, group = group))
  x.1 = x[group=="patient"]
  x.2 = x[group=="control"]
  descdist(x.1,boot=1000)
  plotdist(x.1, "norm", para=list(mean=mean(x.1), sd=sd(x.1)))
  
  print(t.test(x~group,data=data.frame(x=x,group=group)))
  print(wilcox.exact(x~group,data=data.frame(x=x,group=group)))
  
  sd.1 = sd(x.1)
  sd.2 = sd(x.2)
  
  m.1 = mean(x.1)
  m.2 = mean(x.2)
  
  ## find pooled ("averaged" between groups) variance
  sd.p = pool.sd(n.1=table(group)["patient"],sd.1=sd.1,
                 n.2=table(group)["control"],sd.2=sd.2)
  
  ## power for 20 samples per group
  print(power.t.test(20,delta=m.2-m.1,sd = sd.p))
  
  ## compute power for a range of sample sizes and plot
  n = 10:100
  power = sapply(n,function(n) power.t.test(n,delta=m.2-m.1,sd = sd.p)$power)
  plot(power~n)
  
}

## Taxa presence/absence (rates of occurrence) power analysis

occurrence.power.example <- function() {
  m_a = read.table.m_a("seq.raw")
  x = m_a$count[,"Akkermansia"]
  group = m_a$attr[,"DietStatus"]
  
  x = factor(ifelse(x>0,"present","absent"))
  
  t = table(x,group)
  
  print(t)
  
  print(chisq.test(t))
  
  ## Monte Carlo instead of chi-squared approximation
  print(chisq.test(t,simulate.p.value = TRUE, B = 10000))
  
  ##install.packages("lsr")
  ##install.packages("pwr")
  library(lsr)
  library(pwr)
  cramers.v = cramersV(t,simulate.p.value = TRUE, B = 10000)
  
  print(paste("cramers.v",cramers.v))
  
  print(pwr.chisq.test(w=cramers.v,N=40,df=1,sig.level=0.05))
  
  N = 20:100
  power = pwr.chisq.test(w=cramers.v,N=N,df=1,sig.level=0.05)$power
  
  plot(power~N)
}
