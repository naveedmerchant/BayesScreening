library(matrixStats)
logmarg.kern=function(y,x,prior=1){
  
  out=optim(.4,loglike.KGauss,method="L-BFGS-B",lower=.0001,upper=5,y=y,x=x)
  
  h=out$par
  
  cons=-out$val
  
  stat=integrate(integrand.Gauss,lower=0.0001,upper=5,y=y,x=x,cons=cons,prior=prior)$val
  
  list(h,cons+log(stat))
  
}


loglike.KGauss=function(h,y,x){
  
  n=length(x)
  
  m=length(y)
  
  nh=length(h)
  
  M=t(matrix(y,m,n))
  
  llike=1:nh
  
  for(j in 1:nh){
    
    M1=(x-M)/h[j]
    
    M1=dnorm(M1)/h[j]
    
    fhat=as.vector(M1 %*% matrix(1,m,1))/m
    
    fhat[fhat<10^(-320)]=10^(-320)
    
    llike[j]=sum(log(fhat))
    
  }
  
  -llike
  
}



integrand.Gauss=function(h,y,x,cons,prior){
  
  n=length(x)
  
  R=quantile(y,probs=c(.25,.75))
  
  R=R[2]-R[1]
  
  beta=R/1.35
  
  beta1=beta*log(2)/sqrt(qgamma(.5,.5,1))
  
  Prior=beta1*exp(-beta1/h)/h^2
  
  if(prior==1) Prior=(2*beta/sqrt(pi))*h^(-2)*exp(-beta^2/h^2)
  
  arg=-loglike.KGauss(h,y,x)-cons
  
  arg[arg>700]=700
  
  f=exp(arg)*Prior
  
  f
  
}

logmarg.kernMC=function(X1,X2,iter = 10000)
{
  require(matrixStats)
  n <- length(X2)
  k <- length(X1)
  sum1 <- c()  
  sum1list <- c()
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  #browser()
  for(i in 1:iter)
  {
    l <- sample(1:k,n, replace = TRUE)
    sum1 <- sum((X2 - X1[l])^2)
    sum1 <- lgamma((n-1)/2) - ((n-1)/2)*log(.5*(2*B^2 + sum1)) 
    sum1list[i] <- sum1
  }
  logMCinteg <-log(B / sqrt(pi)) - (n/2)*log(2*pi) + logSumExp(sum1list) - log(iter)
  return(logMCinteg)
}


logmarg.kernMCimport=function(X1,X2,iter = 50,importsize = 200)
{
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  Loglist <- c()
  k <- length(X1)
  n <- length(X2)
  for(G in 1:iter)
  {
    
    cauchsamp<- rcauchy(importsize)
    
    poscauchsamp <- cauchsamp[cauchsamp > 0] 
    
    importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))
    
    prodlist1<-c()
    
    for(z in 1:length(poscauchsamp))
    {
      prod <- 0
      for(j in 1:n)
      {
        sum <- 0
        for(i in 1:k)
        {
          sum <- sum + exp(-.5*((X2[j]-X1[i])/poscauchsamp[z])^2)
        }
        prod <- prod + log(sum) + log((1/sqrt(2*pi))) - log((k*poscauchsamp[z]))
      }
      prodlist1[z] <- prod + log(importancepart[z])
    }
    Loglist[G] <- logSumExp(prodlist1)
  }
  return(Loglist)
}

logmarg.specialkernMCimport=function(X1,X2,iter = 50,importsize = 200)
{
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  Loglist <- c()
  k <- length(X1)
  n <- length(X2)
  for(G in 1:iter)
  {
    
    cauchsamp<- rcauchy(importsize)
    
    poscauchsamp <- abs(cauchsamp)
    
    importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))
    
    prodlist1<-c()
    
    for(z in 1:length(poscauchsamp))
    {
      prod <- 0
      for(j in 1:n)
      {
        sum <- 0
        for(i in 1:k)
        {
          sum = sum + (((8*pi*exp(1))^.5)*pnorm(1))^(-1)*exp(-.5*(log(1+abs(X2[j] - X1[i])/poscauchsamp[z])^2))
        }
        prod <- prod + log(sum) + log((1/sqrt(2*pi))) - log((k*poscauchsamp[z]))
      }
      prodlist1[z] <- prod + log(importancepart[z])
    }
    Loglist[G] <- logSumExp(prodlist1)
  }
  return(Loglist)
}


logmarg.specialkernMCimport2=function(X1,X2,iter = 50,importsize = 200)
{
  R = quantile(X1,probs=c(.25,.75))
  R=R[2]-R[1]
  R = unname(R)
  B = R / 1.35
  Loglist <- c()
  k <- length(X1)
  n <- length(X2)
  for(G in 1:iter)
  {
    
    cauchsamp<- rt(importsize, df = 4)
    
    poscauchsamp <- abs(cauchsamp)
    
    importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dt(poscauchsamp, df = 4))
    
    prodlist1<-c()
    
    for(z in 1:length(poscauchsamp))
    {
      prod <- 0
      for(j in 1:n)
      {
        sum <- 0
        for(i in 1:k)
        {
          sum = sum + (((8*pi*exp(1))^.5)*pnorm(1))^(-1)*exp(-.5*(log(1+abs(X2[j] - X1[i])/poscauchsamp[z])^2))
        }
        prod <- prod + log(sum) + log((1/sqrt(2*pi))) - log((k*poscauchsamp[z]))
      }
      prodlist1[z] <- prod + log(importancepart[z])
    }
    Loglist[G] <- logSumExp(prodlist1)
  }
  return(Loglist)
}

logmarg.kernH=function(y,x){
  
  out=optim(.4,loglike.KHall,method="L-BFGS-B",lower=.0001,upper=5,y=y,x=x)
  
  h=out$par
  
  cons=-out$val
  
  stat=integrate(integrand.Hall,lower=0.0001,upper=5,y=y,x=x,cons=cons,hhat=h)$val
  
  list(h,cons+log(stat))
  
}



laplace.kernH=function(y,x,hhat){
  
  n=length(x)
  
  start=1.144*(IQR(x)/1.35)*n^(-1/5)
  
  out=optim(start,logintegrand.Hall,method="L-BFGS-B",lower=.0001,upper=5,hessian=T,y=y,x=x,hhat=hhat)
  
  h=out$par
  
  cons=-out$val
  
  hess=abs(out$hessian)
  laplace=0.5*log(2*pi)-0.5*log(hess)+cons
  c(laplace,h,hess)
  
}

laplace.kernH2=function(y,x,hhat){
  
  n=length(x)
  
  #start=1.144*(IQR(x)/1.35)*n^(-1/5)
  
  #out=optim(start,logintegrand.Hall,method="L-BFGS-B",lower=.0001,upper=5,hessian=T,y=y,x=x,hhat=hhat)
  #browser()
  out = hessian(f = function(h) {logintegrand.Hall(h,y=y,x=x,hhat = hhat)}, x = hhat, pert = 10^(-7))
  
  cons = -logintegrand.Hall(hhat, y = y, x = x, hhat = hhat)
  
  hess=abs(out)
  
  laplace=0.5*log(2*pi)-0.5*log(hess)+cons
  
  c(laplace,hhat,hess)
  
}

logintegrand.Gauss=function(h,y,x,hhat){
  
  n=length(x)
  
  beta=hhat
  
  Prior=log(2*beta/sqrt(2*pi))-beta^2/h^2-2*log(h)
  
  f=loglike.KHall(h,y,x)-Prior
  
  f
  
}

laplace.kernG2=function(y,x,hhat){
  #browser()
  n=length(x)
  
  #start=1.144*(IQR(x)/1.35)*n^(-1/5)
  
  #out=optim(start,logintegrand.Hall,method="L-BFGS-B",lower=.0001,upper=5,hessian=T,y=y,x=x,hhat=hhat)
  #browser()
  out = hessian(f = function(h) {logintegrand.Gauss(h,y=y,x=x,hhat = hhat)}, x = hhat, pert = 10^(-5))
  
  cons = -logintegrand.Gauss(hhat, y = y, x = x, hhat = hhat)
  
  hess=abs(out)
  
  laplace=0.5*log(2*pi)-0.5*log(hess)+cons
  
  c(laplace,hhat,hess)
  
}





KHall=function(x){
  
  con=sqrt(8*pi*exp(1))*pnorm(1)
  
  K=exp(-0.5*(log(1+abs(x)))^2)/con
  
  K
  
}



integrand.Hall=function(h,y,x,cons,hhat){
  
  n=length(x)
  
  beta=hhat
  
  Prior=(2*beta/sqrt(2*pi))*exp(-beta^2/h^2)/h^2
  
  arg=-loglike.KHall(h,y,x)-cons
  
  arg[arg>700]=700
  
  f=exp(arg)*Prior
  
  f
  
}



logintegrand.Hall=function(h,y,x,hhat){
  
  n=length(x)
  
  beta=hhat
  
  Prior=log(2*beta/sqrt(2*pi))-beta^2/h^2-2*log(h)
  
  f=loglike.KHall(h,y,x)-Prior
  
  f
  
}


logmarg.specialkernMCimportnewprior=function(X1,X2,hhat,iter = 50,importsize = 200)
{
  B = hhat
  Loglist <- c()
  k <- length(X1)
  n <- length(X2)
  for(G in 1:iter)
  {
    
    cauchsamp<- rcauchy(importsize)
    
    poscauchsamp <- abs(cauchsamp)
    
    importancepart <- ((2*B*(1/poscauchsamp^2)*(exp(-(B^2)/(poscauchsamp)^2)) / sqrt(pi))) / (2*dcauchy(poscauchsamp))
    
    prodlist1<-c()
    
    for(z in 1:length(poscauchsamp))
    {
      prod <- 0
      for(j in 1:n)
      {
        sum <- 0
        for(i in 1:k)
        {
          sum = sum + (((8*pi*exp(1))^.5)*pnorm(1))^(-1)*exp(-.5*(log(1+abs(X2[j] - X1[i])/poscauchsamp[z])^2))
        }
        prod <- prod + log(sum) + log((1/sqrt(2*pi))) - log((k*poscauchsamp[z]))
      }
      prodlist1[z] <- prod + log(importancepart[z])
    }
    Loglist[G] <- logSumExp(prodlist1)
  }
  return(Loglist)
}

#I will speed the lower function up later
#Computes log BF using PolyaTreetest
#Things to specify are inverse function, c, and depth level (leveltot)

PolyaTreetest <- function(datasetX, datasetY, Ginv = NULL, c = NULL, leveltot = NULL)  
{
  #A lot of places can be sped up I think by predefining epsilon list and epsilonlist2.
  #The loops themselves also take some time, maybe we can throw the whole thing into C, but the problem is that
  #We have a jagged matrix that's essentially run with a list rn, changing that to be in C might be trickier than it sounds.
  #But it is something I wanted to look into
  Jointset = c(datasetX, datasetY)
  if(is.null(Ginv))
  {
    meanJ = mean(Jointset)
    sdJ = sd(Jointset)
    Ginv = function(x){qnorm(x, mean = meanJ, sd = sdJ)}
  }
  if(is.null(c))
  {
    c = 1
  } else{
    if(c < 0)
    {
      warning("c must be some positive number, proceeding as if c = 1 was inputted")
      c = 1
    }
  }
  n = max(length(datasetX), length(datasetY))
  if(is.null(leveltot))
  {
    leveltot = ceiling(log2(n))
  }
  else{
    if(leveltot < 1)
    {
      warning("leveltot must be some positive integer that represents depth of tree, proceeding at default choice for leveltot.")
      leveltot = ceiling(leveltot = log2(n))
    }
    leveltot = ceiling(leveltot) #Forces leveltot to be an integer
  }
  alphalist = list()
  for(m in 1:leveltot)
  {
    alphalist[[m]] = rep(c*m^2, 2^m)
  }
  
  epsilonlist = c()
  epsilonlist2 = list()
  splitlist = list()
  for(j in 1:leveltot)
  {
    splitlist[[j]] = rep(0, 2^j)
  }
  
  for(j in 1:length(datasetX))
  {
    thresh = .5
    for(m in 1:leveltot)
    {
      epsilonlist[m] = datasetX[j] > Ginv(thresh)
      multivec = 2^seq(m-1, 0, -1)
      ind = sum(epsilonlist[1:m]*multivec) + 1
      splitlist[[m]][ind] = splitlist[[m]][ind] + 1
      thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
    }
    epsilonlist2[[j]] <- epsilonlist
    epsilonlist = c()
  }
  epsilonlist3 = list()
  splitlist2 = list()
  for(j in 1:leveltot)
  {
    splitlist2[[j]] = rep(0,2^j)
  }
  for(j in 1:length(datasetY))
  {
    thresh = .5
    for(m in 1:leveltot)
    {
      epsilonlist[m] = datasetY[j] > Ginv(thresh)
      multivec = 2^seq(m-1, 0, -1)
      ind = sum(epsilonlist[1:m] * multivec) + 1
      splitlist2[[m]][ind] = splitlist2[[m]][ind] + 1
      thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
    }
    epsilonlist3[[j]] <- epsilonlist
    epsilonlist = c()
  }
  
  bj = rep(0, times = leveltot)
  
  for(j in 1:leveltot)
  {
    k = 1
    for(D in 1:(length(alphalist[[j]]) - 1))
    {
      bj[j] = bj[j] + lbeta(alphalist[[j]][k],alphalist[[j]][k+1]) +
        lbeta(alphalist[[j]][k] + splitlist[[j]][k] +
                splitlist2[[j]][k],alphalist[[j]][k+1] +
                splitlist[[j]][k+1] + splitlist2[[j]][k+1] ) -
        lbeta(alphalist[[j]][k] + splitlist[[j]][k], alphalist[[j]][k+1] + splitlist[[j]][k+1]) -
        lbeta(alphalist[[j]][k] + splitlist2[[j]][k], alphalist[[j]][k+1] + splitlist2[[j]][k+1])
      k = k + 2
      if(k + 1 > (length(alphalist[[j]])))
        break
    }
  }
  return(list(logBF = -sum(bj), logBFcont =  -bj))
}



loglike.KHall=function(h,y,x){
  n=length(x)
  m=length(y)
  nh=length(h)
  M=t(matrix(y,m,n))
  llike=1:nh
  for(j in 1:nh){
    M1=(x-M)/h[j]
    M1=KHall(M1)/h[j]
    fhat=as.vector(M1 %*% matrix(1,m,1))/m
    fhat[fhat<10^(-320)]=10^(-320)
    llike[j]=sum(log(fhat))
  }
  -llike
}



PolyaTreetestlastsplit = function(datasetX,datasetY,Ginv = qnorm,c = 1, leveltot = 9)  
{
  alphalist = list()
  for(m in 1:leveltot)
  {
    alphalist[[m]] = rep(c*m^2,2^m)
  }
  
  #An alternative way to write this?
  epsilonlist = c()
  epsilonlist2 = list()
  splitlist = list()
  for(j in 1:leveltot)
  {
    splitlist[[j]] = rep(0,2^j)
  }
  
  for(j in 1:length(datasetX))
  {
    thresh = .5
    for(m in 1:leveltot)
    {
      epsilonlist[m] = datasetX[j] > Ginv(thresh)
      multivec = 2^seq(m-1, 0, -1)
      ind = sum(epsilonlist[1:m]*multivec) + 1
      splitlist[[m]][ind] = splitlist[[m]][ind] + 1
      thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
    }
    epsilonlist2[[j]] <- epsilonlist
    epsilonlist = c()
  }
  epsilonlist3 = list()
  splitlist2 = list()
  for(j in 1:leveltot)
  {
    splitlist2[[j]] = rep(0,2^j)
  }
  for(j in 1:length(datasetY))
  {
    thresh = .5
    for(m in 1:leveltot)
    {
      epsilonlist[m] = datasetY[j] > Ginv(thresh)
      multivec = 2^seq(m-1, 0, -1)
      ind = sum(epsilonlist[1:m]*multivec) + 1
      splitlist2[[m]][ind] = splitlist2[[m]][ind] + 1
      thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
    }
    epsilonlist3[[j]] <- epsilonlist
    epsilonlist = c()
  }
  #splitlist2
  bj = c()
  
  for(j in 1:leveltot)
  {
    for(D in 1:(length(alphalist[[j]]) - 1))
    {
      bj[j] = lbeta(alphalist[[j]][D],alphalist[[j]][D+1]) +
        lbeta(alphalist[[j]][D] + splitlist[[j]][D] +
                splitlist2[[j]][D],alphalist[[j]][D+1] +
                splitlist[[j]][D+1] + splitlist2[[j]][D+1] ) -
        lbeta(alphalist[[j]][D] + splitlist[[j]][D], alphalist[[j]][D+1] + splitlist[[j]][D+1]) -
        lbeta(alphalist[[j]][D] + splitlist2[[j]][D], alphalist[[j]][D+1] + splitlist2[[j]][D+1])
      
    }
  }
  return(sum(bj))
  
}

laplace.kernH2c = function(y,x,hhat,c){
  
  n=length(x)
  
  #start=1.144*(IQR(x)/1.35)*n^(-1/5)
  
  #out=optim(start,logintegrand.Hall,method="L-BFGS-B",lower=.0001,upper=5,hessian=T,y=y,x=x,hhat=hhat)
  #browser()
  out = hessian(f = function(h) {logintegrand.Hall(h,y=y,x=x,hhat = hhat)}, x = hhat, pert = 10^(-7))
  
  cons = c
  
  hess=abs(out)
  
  laplace=0.5*log(2*pi)-0.5*log(hess)+cons
  
  c(laplace,hhat,hess)
  
}
