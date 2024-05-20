#Laplace approximation functions to evaluate integrals

library(rootSolve)

GaussKernel = function(h,datagen2,x)
{
  sum = 0
  for(i in 1:length(datagen2))
  {
    sum = sum + (1/sqrt(2*pi)) * exp(-.5*((x - datagen2[i])/h)^2)
  }
  return((1/(length(datagen2) * h)) * sum)
}


HallKernel = function(h,datagen2,x)
{
  sum = 0
  for(i in 1:length(datagen2))
  {
    sum = sum + (((8*pi*exp(1))^.5)*pnorm(1))^(-1)*exp(-.5*(log(1+abs(x - datagen2[i])/h)^2))
  }
  return((1/(length(datagen2) * h)) * sum)
}

logpriorused <- function(h,x)
{
  R=quantile(x,probs=c(.25,.75))
  
  R=R[2]-R[1]
  
  R = unname(R)
  
  beta=R/1.35
  
  beta1=beta*log(2)/sqrt(qgamma(.5,.5,1))
  
  Prior = log(2*beta) - .5*log(pi) - 2*log(h) - (beta^2 / h^2)
  
  return(Prior)
  
}

priorused <- function(h,x)
{
  R=quantile(x,probs=c(.25,.75))
  
  R=R[2]-R[1]
  
  R = unname(R)
  
  beta=R/1.35
  
  Prior = log(2*beta) - .5*log(pi) - 2*log(h) - (beta^2 / h^2)
  
  Prior = exp(Prior)
  
  return(Prior)
  
}

logpriorusedlik2 <- function(h,x,hhat)
{
  beta = hhat
  
  Prior = log(2*beta) - .5*log(pi) - 2*log(h) - (beta^2 / h^2)
  
  return(Prior)
  
}