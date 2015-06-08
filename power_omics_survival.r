library(powerSurvEpi)
library(Hmisc)


cpower.example <- function(
  morts = seq(40,65,length=50),
  mort.n = 50,
  red = c(25,35,45,60),
  tref = 3,
  accrual = 1,
  tmin = NULL,
  alpha = 0.05,
  n.power = 200,
  power.n = 0.8,
  noncomps = c(0,3,5,10),
  noncomp.n = 10,
  event.name = "Recurrence"
  ) {
  #In this example, 4 plots are drawn on one page, one plot for each
  #combination of noncompliance percentage.  Within a plot, the
  #3-year recurrence % in the control group is on the x-axis, and
  #separate curves are drawn for several % reductions in recurrence
  #with the intervention.  The accrual period is 1y, with all
  #patients followed at least 3y.
  
  if(is.null(tmin)) {
    tmin = tref
  }
  
  par(mfrow=c(2,2),oma=c(3,0,3,0))
  
  
  for(noncomp in noncomps) {
    nc.i <- nc.c <- noncomp
    z <- paste("Drop-in ",nc.c,"%, Non-adherence ",nc.i,"%",sep="")
    plot(0,0,xlim=range(morts),ylim=c(0,1),
         xlab=sprintf("%s-year %s in Control Patients (%%)",tref,event.name),
         ylab="Power",type="n")
    title(z)
    cat(z,"\n")
    lty <- 0
    for(r in red) {
      lty <- lty+1
      power <- morts
      i <- 0
      for(m in morts) {
        i <- i+1
        power[i] <- cpower(tref=tref, n=n.power, mc = m/100, r = r, 
                           accrual = accrual, tmin = tmin, 
                           noncomp.c=nc.c, noncomp.i=nc.i, pr=T)
      }
      lines(morts, power, lty=lty)
      text(x=morts[length(morts)],y=power[length(power)],paste(r,"% red.",sep=""),cex=0.8)
    }
    #if(noncomp==0)legend(x="topright",legend=rev(paste(red,"% reduction",sep="")),bty="n",cex=0.8)
  }
  mtitle("Power vs Non-Adherence for Main Comparison",
         ll=sprintf("alpha=%s, 2-tailed, Total N=%s",alpha,n.power),cex.l=.8)
  #
  # Point sample size requirement vs. recurrence reduction
  # Root finder (uniroot()) assumes needed sample size is between
  # 100 and 4000
  #
  nc.i <- nc.c <- noncomp.n 
  red <- seq(min(red),max(red),by=.25)
  samsiz <- red
  
  
  i <- 0
  for(r in red) {
    i <- i+1
    samsiz[i] <- uniroot(function(x) cpower(tref=tref, n=x, mc = mort.n/100, r = r, 
                                            accrual = accrual, tmin = tmin, 
                                            noncomp.c=nc.c, noncomp.i=nc.i, pr=FALSE) - power.n,
                         c(n.power/50,n.power*50))$root
  }
  
  
  par(mfrow=c(1,1))
  plot(red, samsiz, xlab=sprintf('%% Reduction in %s-Year %s',tref,event.name),
       ylab='Total Sample Size', type='n')
  lines(red, samsiz, lwd=2)
  title(paste(sprintf('Sample Size for Power=%s\nDrop-in %s%%, Non-adherence %s%%',power.n,nc.c,nc.i),
                      sprintf("%s-year %s in Control Patients %s%%",tref,event.name,mort.n),
                      sep="\n"
                ))
  title(sub=sprintf('alpha=%s, 2-tailed',alpha), adj=0)
  
}
