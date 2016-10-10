#------------------------------------------------------
#Supporting Information


#Estimating and Comparing Microbial Diversity in the Presence of Sequencing Errors

#Chun-Huo Chiu and Anne Chao

#Institute of Statistics, National Tsing Hua University, Hsin-Chu, Taiwan, 30043



#Supplemental Text S1: R codes for obtaining estimators of Hill numbers 

#------------------------------------------------------

# R scripts for computing diversity (Hill numbers) profile using individual-based abundance data or sampling-unit-based incidence data.
# In all functions, param x is a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (for incidence data).
# For incidence data, the first entry of x must be the number of sampling units. 
# In all functions, param q is the diversity order; the suggested range for q is [0, 3].

#-----------------------------------------------
# Diversity profile estimator (abundance data)
#-----------------------------------------------
#' Chao_Hill_abu(x, q) is a function of obtaining estimators of Hill numbers of order q based on abundance data.
#' @param x a vector of species sample frequencies. 
#' @param q a numeric or a vector of diversity order; The suggested range for q is [0, 3].
#' @return a numerical vector of diversity. 

Chao_Hill_abu = function(x,q){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      #UE <- sum(x/n*(digamma(n)-digamma(x)))
      #A = 1 - ifelse(f2 > 0, (n-1)*f1/((n-1)*f1+2*f2), (n-1)*f1/((n-1)*f1+2))
      #B = f1/n*(1-A)^(-n+1)*(-log(A)-sum(sapply(1:(n-1),function(k){1/k*(1-A)^k})))
      #exp(UE+B)      
      A <- sum(x/n*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0||p1==1,0,f1/n*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,A^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      (A+B)^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#-----------------------------------------------
# Diversity profile estimator (incidence data)
#-----------------------------------------------
#' Chao_Hill_inc(x, q) is a function of obtaining estimators of Hill numbers of order q based on incidence data.
#' @param x a vector of species incidence-based sample frequencies. The first entry of x must be the number of sampling units.
#' @param q a numeric or a vector of diversity order.
#' @return a numerical vector.

Chao_Hill_inc = function(x,q){
  n = x[1]
  x = x[-1];x = x[x>0]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  p1 = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  r <- 0:(n-1)
  Sub <- function(q){
    if(q==0){
      sum(x>0) + (n-1)/n*ifelse(f2>0, f1^2/2/f2, f1*(f1-1)/2)
    }
    else if(q==1){
      A <- sum(x/U*(digamma(n)-digamma(x)))
      B <- ifelse(f1==0|p1==1,0,f1/U*(1-p1)^(1-n)*(-log(p1)-sum(sapply(1:(n-1), function(r)(1-p1)^r/r))))
      exp(A+B)*U/n
    }else if(abs(q-round(q))==0){
      A <- sum(exp(lchoose(x,q)-lchoose(n,q)))
      ifelse(A==0,NA,((n/U)^q*A)^(1/(1-q)))
    }else {
      sort.data = sort(unique(x))
      tab = table(x)
      term = sapply(sort.data,function(z){
        k=0:(n-z)
        sum(choose(k-q,k)*exp(lchoose(n-k-1,z-1)-lchoose(n,z)))
      })
      A = sum(tab*term)
      B = ifelse(f1==0|p1==1,0,f1/n*(1-p1)^(1-n)*(p1^(q-1)-sum(choose(q-1,r)*(p1-1)^r)))
      ((n/U)^q*(A+B))^(1/(1-q))
    }
  }
  sapply(q, Sub)
}

#' Chao_Hill(x, q,datatype) combines Chao_Hill_abu and Chao_Hill_inc given a specified datatype (either abundance data or incidence data).
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Chao_Hill = function(x,q,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    est = Chao_Hill_abu(x,q)
  }else{
    est = Chao_Hill_inc(x,q)
  }
  return(est)
}

#-----------------------
# The empirical profile 
#-----------------------
#' Hill(x, q, datatype) is a function of obtaining the empirical Hill numbers of order q based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numerical vector.

Hill <- function(x,q,datatype = c("abundance","incidence")){
  if(datatype=="incidence"){x = x[-1]}
  p <- x[x>0]/sum(x)
  Sub <- function(q){
    if(q==0) sum(p>0)
    else if(q==1) exp(-sum(p*log(p)))
    else exp(1/(1-q)*log(sum(p^q)))
  }
  sapply(q, Sub)
}

#-----------------------------------------
# The bootstrap method for obtaining s.e. 
#-----------------------------------------
#' Bt_prob_abu(x) is a function of estimating the species probabilities in the bootstrap assemblage based on abundance data.
#' @param x a vector of species sample frequencies.
#' @return a numeric vector.

Bt_prob_abu = function(x){
  x = x[x>0]
  n = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  C = 1 - f1/n*ifelse(f2>0,(n-1)*f1/((n-1)*f1+2*f2),ifelse(f1>0,(n-1)*(f1-1)/((n-1)*(f1-1)+2),0))
  W = (1-C)/sum(x/n*(1-x/n)^n)
  
  p.new = x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0 = (1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob_inc(x) is a function of estimating the species incidence probabilities in the bootstrap assemblage based on incidence data.
#' @param x a vector of incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @return a numeric vector.

Bt_prob_inc = function(x){
  n = x[1]
  x = x[-1]
  U = sum(x)
  f1 = sum(x==1)
  f2 = sum(x==2)
  A = ifelse(f2>0,2*f2/((n-1)*f1+2*f2),ifelse(f1>0,2/((n-1)*(f1-1)+2),1))
  C=1-f1/U*(1-A)
  W=U/n*(1-C)/sum(x/n*(1-x/n)^n)
  
  p.new=x/n*(1-W*(1-x/n)^n)
  f0 = ceiling(ifelse(f2>0,(n-1)/n*f1^2/(2*f2),(n-1)/n*f1*(f1-1)/2))
  p0=U/n*(1-C)/f0
  p.new=c(p.new,rep(p0,f0))
  return(p.new)
}

#' Bt_prob(x,datatype) combines the two functions Bt_prob_abu and Bt_prob_inc for a specified datatype. 
#' @param x a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @return a numeric vector.

Bt_prob = function(x,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    prob = Bt_prob_abu(x)
  }else{
    prob = Bt_prob_inc(x)
  }
  return(prob)
}

#' Bootstrap.CI(x,q,B,datatype,conf) is a function of calculating the bootsrapping standard error based on abundance data or incidence data.
#' @param x a vector of species sample frequencies (for abundance data) or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param q a numeric or a vector of diversity order.
#' @param B an integer to specify the number of replications in the bootstrap procedure, B = 1000 is suggested for constructing confidence intervals; 
#'  To save running time, use a smaller value (e.g. B = 200)..
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 3 matrices including respectively the difference between the average and lower confidence bound of the B bootstrap estimates, 
#'  the difference between the upper confidence bound and the average of the B bootstrap estimates, and the bootstrap standard error of the diversity estimate.
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

Bootstrap.CI = function(x,q,B = 1000,datatype = c("abundance","incidence"),conf = 0.95){
  datatype = match.arg(datatype,c("abundance","incidence"))
  p.new = Bt_prob(x,datatype)
  n = ifelse(datatype=="abundance",sum(x),x[1])
  # set.seed(456)
  if(datatype=="abundance"){
    data.bt = rmultinom(B,n,p.new)

    I1=which(data.bt==1);data.bt[I1]=rep(1,length(I1));         
    F2=sapply(1:B,function(k) as.double(sum(data.bt[,k]==2)));
    F3=sapply(1:B,function(k) as.double(sum(data.bt[,k]==3)));
    F4=sapply(1:B,function(k) as.double(sum(data.bt[,k]==4)));
    F1=sapply(1:B,function(k) ifelse(F3[k]*F4[k]>0, 4*F2[k]^2/(3*F3[k])-F2[k]*F3[k]/(2*F4[k]), 4*F2[k]^2/(3*(F3[k]+1))-F2[k]*F3[k]/(2*(F4[k]+1))));
    F1=round(F1);Max=max(F1);F1.Mat=matrix(0,ncol=B,nrow=Max);
    sapply(1:B,function(k) F1.Mat[,k][1:length(F1[k])]=rep(1,length(F1[k])) )
    data.bt=rbind(data.bt,F1.Mat);

  }else{
    data.bt = rbinom(length(p.new)*B,n,p.new) 
    data.bt = matrix(data.bt,ncol=B)
    data.bt = rbind(rep(n,B),data.bt)
  }
  
  mle = apply(data.bt,2,function(x)Hill(x,q,datatype))
  pro = apply(data.bt,2,function(x)Chao_Hill(x,q,datatype))
  
  mle.mean = rowMeans(mle)
  pro.mean = rowMeans(pro)
  
  LCI.mle =  -apply(mle,1,function(x)quantile(x,probs = (1-conf)/2)) + mle.mean
  UCI.mle = apply(mle,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - mle.mean
  
  LCI.pro =  -apply(pro,1,function(x)quantile(x,probs = (1-conf)/2)) + pro.mean
  UCI.pro = apply(pro,1,function(x)quantile(x,probs = 1-(1-conf)/2)) - pro.mean
  
  LCI = rbind(LCI.mle,LCI.pro)
  UCI = rbind(UCI.mle,UCI.pro)
  
  sd.mle = apply(mle,1,sd)
  sd.pro = apply(pro,1,function(x)sd(x,na.rm = T))
  se = rbind(sd.mle,sd.pro)
  
  return(list(LCI=LCI,UCI=UCI,se=se))
  
}

#----------------------
# singleton estimator
#----------------------
#' singleton.Est(x) is a function of obtaining the estimator of singleton based on corrected data
#' @param x a vector of species sample frequencies 
#' @return singleton estimate and a numerical vector of corrected data. 
Chao_Hill = function(x,q,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  if(datatype == "abundance"){
    est = Chao_Hill_abu(x,q)
  }else{
    est = Chao_Hill_inc(x,q)
  }
  return(est)
}

singleton.Est=function(dat,datatype = c("abundance","incidence")){
  datatype = match.arg(datatype,c("abundance","incidence"))
  
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)

  if(datatype == "abundance"){
    #convert to double to avoid overflow after product is taken
       f2 <- as.double(sum(dat == 2)) 	#doubleton
       f3 <- as.double(sum(dat == 3))
       f4 <- as.double(sum(dat == 4))
       f1= ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
       I=which(dat==1);dat=dat[-I];
       dat=c(dat,rep(1,round(f1)));
       dat=dat[dat>0];      
  }else{
       N=dat[1];dat=dat[-1];   
       f2 <- as.double(sum(dat == 2)) 	#doubleton
       f3 <- as.double(sum(dat == 3))
       f4 <- as.double(sum(dat == 4))
       f1= ifelse(f3*f4>0, 4*f2^2/(3*f3)-f2*f3/(2*f4), 4*f2^2/(3*(f3+1))-f2*f3/(2*(f4+1)))
       I=which(dat==1);dat=dat[-I];
       dat=c(dat,rep(1,round(f1)));
       dat=dat[dat>0];
       dat=c(N,dat); 
       }      
 return(list(singleton.est=f1,corrected.data=dat)) 
}


#------------------------
# Main function ChaoHill
#------------------------
#' ChaoHill(dat, datatype, from, to, interval, B, conf) is the function of calculating the empirical and the proposed diversity profile, 
#' their bootsrap standard errors and confidance intervals.
#' @param dat a vector of species sample frequencies (for abundance data), or incidence-based sample frequencies (1st entry must be the number of sampling unit).
#' @param datatype a character of data type,"abundance" or "incidence".
#' @param from a numeric number of diversity order q (the start order of profile).
#' @param to a numeric number of diversity order q (the end order of profile).
#' @param interval a numeric number to specify each increment of q from the start to end order.
#' @param B an integer to specify the number of bootstrap replications, B = 1000 is suggested for constructing confidence intervals; 
#'  To save running time, use a smaller value (e.g. B = 200).
#' @param conf a confidence coefficient between 0 and 1.
#' @return a list, consisting of 4 matrices including respectively diversity estimates, bootstrap standard errors, lower confidence bounds, and upper confidence bounds. 
#'  In each matrix, the first row gives the results for the empirical diversity, and the second row gives the results for the proposed diversity estimates.
#'  Columns give the results for different orders of q.

singleton.ChaoHill <- function(dat, datatype=c("abundance", "incidence"), from=0, to=3, interval=0.1, B=1000, conf=0.95){ 
  datatype = match.arg(datatype,c("abundance","incidence"))
  # for real data estimation
  
  if (is.matrix(dat) == T || is.data.frame(dat) == T){
    if (ncol(dat) != 1 & nrow(dat) != 1)
      stop("Error: The data format is wrong.")
    if (ncol(dat) == 1){
      dat <- dat[, 1]
    } else {
      dat <- dat[1, ]
    }
  }
  dat <- as.numeric(dat)

 #------------------
 # correct singleton and raw data
 #------------------   
   temp=singleton.Est(dat,datatype);
   f1=temp$singleton.est;
   dat=temp$corrected.data       
 #-----------------------

   q <- seq(from, to, by=interval)
  
  #-------------
  #Estimation
  #-------------
  MLE=Hill(dat,q,datatype)
  
  qD_pro=Chao_Hill(dat,q,datatype)
  
  qD012=Chao_Hill(dat,q=c(0,1,2),datatype); 

  CI_bound = Bootstrap.CI(dat,q,B,datatype,conf)
  se = CI_bound$se

  #-------------------
  #Confidence interval
  #-------------------
  tab.est=data.frame(rbind(MLE,qD_pro))
  
  LCI <- tab.est - CI_bound$LCI
  UCI <- tab.est + CI_bound$UCI
  
  colnames(tab.est) <- colnames(se) <- colnames(LCI) <- colnames(UCI) <- paste("q = ", q, sep="")    
  rownames(tab.est) <- rownames(se) <- rownames(LCI) <- rownames(UCI) <- c("Observed", "Chao_2013")
  return(list(Singleton.est=f1,
              EST = tab.est, SD = se,
              LCI = LCI,UCI = UCI))
  
}




######Example: Swine-feces data of Table 1 in the main text 


# run all R functions above. 

singleton.ChaoHill.example <- function() {

fi=c(8025,605,129,41,16,8,4,2,1,1,1) # OTU frequency counts of swine-feces data

data=rep(1:length(fi),fi) # transform frequency counts data to species abundances 

output=singleton.ChaoHill(data,"abundance",from=0,to=2,interval=1,B=200);

#estimates of Hill numbers (order q= 0, 1, 2) and corresponding 95% confidence interval


"
output

$Singleton.est
[1] 2831.436

$EST
             q = 0    q = 1    q = 2
Observed   3639.00 3250.379 2741.880
Chao_2013 10261.22 9081.061 6404.025

$SD
             q = 0     q = 1     q = 2
Observed   24.5368  33.72911  54.18838
Chao_2013 278.8550 178.56957 198.79723

$LCI
             q = 0    q = 1    q = 2
Observed  3591.840 3181.323 2630.539
Chao_2013 9696.206 8736.525 6006.244

$UCI
             q = 0    q = 1    q = 2
Observed   3682.89 3304.304 2837.815
Chao_2013 10814.14 9427.188 6768.070
"

##### Rarefaction and Extrapolation of Hill numbers 
#First correct the raw singleton count by running singleton.Est function.
data=singleton.Est(data,"abundance")$corrected.data;

#Then, using R package "iNEXT" to plot the rarefaction and extrapolation curve as shown in Fig. 2 and Fig. 3. 

}
