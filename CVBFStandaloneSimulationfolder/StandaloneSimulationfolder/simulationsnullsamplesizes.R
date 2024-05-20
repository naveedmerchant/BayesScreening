install.packages("devtools")
devtools::install_github("naveedmerchant/BayesScreening")
library(BSCRN)

library(stats)
#install.packages("rootSolve")
library(parallel)
library(foreach)
library(doRNG)
library(rootSolve)

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

KHall=function(x){
  
  con=sqrt(8*pi*exp(1))*pnorm(1)
  
  K=exp(-0.5*(log(1+abs(x)))^2)/con
  
  K
  
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

logintegrand.Hall=function(h,y,x,hhat){
  
  n=length(x)
  
  beta=hhat
  
  Prior=log(2*beta/sqrt(2*pi))-beta^2/h^2-2*log(h)
  
  f=loglike.KHall(h,y,x)-Prior
  
  f
  
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


set.seed(100)
samplesizes = c(200,400,800)
#samplesizes = 150
trainsizes = c(50,75,112)
reps = 1500
logCVBFs = matrix(nrow = reps, ncol = length(samplesizes))
logPTBFs = matrix(nrow = reps, ncol = length(samplesizes))
CVBFforAverage = rep(0, times = 30)
logBFs = matrix(nrow = reps*length(samplesizes), ncol = 3)
ncores = 56
for(j in 1:length(samplesizes))
{
  cl <- makeCluster(ncores) 
  registerDoParallel(cl)
  logBFsj = foreach(k = 1:reps, .combine = rbind, .export = c("PolyaTreetest", "CVBFtestrsplit", "HallKernel", "laplace.kernH2c", "KHall", "logpriorused", "logintegrand.Hall", "loglike.KHall", "samplesizes", "trainsizes", "j"), .packages = "rootSolve") %dorng% {
    dataset1 = rnorm(samplesizes[j])
    dataset2 = rnorm(samplesizes[j])
    CVBFforAverage = rep(0, times = 30)
    for(l in 1:30)
    {
      trainset1 = sample(1:samplesizes[j], size = trainsizes[j])
      trainset2 = sample(1:samplesizes[j], size = trainsizes[j])
      CVBFforAverage[l] = CVBFtestrsplit(dataset1 = dataset1, dataset2 = dataset2, train1_ids = trainset1, train2_ids = trainset2, trainsize1 = trainsizes[j], trainsize2 = trainsizes[j])$logBF
    }
    logCVBFs = mean(CVBFforAverage)
    logPTBFs = PolyaTreetest(datasetX = dataset1, datasetY = dataset2, Ginv = qnorm)$logBF
    BFs = c(logCVBFs, logPTBFs)
    BFs
  }
  stopCluster(cl)
  logBFsj = unname(logBFsj)
  logBFsj = logBFsj[1:reps, 1:2]
  logBFs[(1 + (j-1)*reps):(j*reps), 1:2] = logBFsj
  logBFs[(1 + (j-1)*reps):(j*reps), 3] = rep(trainsizes[j], times = reps)
  print(j)
}

#This might take sometime to run. A r data file with the output. Note that this also was run with 50 ish cores, run with the number appropriate for your computer

head(logBFs)




#First col of logBFs is CVBF, second col is PTBF, third row is sample size
#.export = c("CVBFtestrsplit", "HallKernel", "laplace.kernH2c", "KHall", "logpriorused", "logintegrand.Hall", "loglike.KHall")


# cl <- makeCluster(ncores) 
# registerDoParallel(cl)
# logBFlist <- foreach(j=1:p, .combine=c, .export = c("PolyaTreetest")) %dopar% {
#   tempMatrix = PolyaTreetest(datasetX[Class0ind, j], datasetX[-Class0ind, j], Ginv = Ginv, c = c, leveltot = leveltot)$logBF 
#   
#   tempMatrix #Equivalent to logBFlist = c(logBFlist, tempMatrix)
# }
# stopCluster(cl)
# ImportantVars = list(logBFlist = logBFlist, varspicked = which(logBFlist > cutoff))
# #ImportantVars = which(logBFlist > cutoff)

library(ggplot2)
datasetsizes = c(rep(200, times = 1500), rep(400, times = 1500), rep(800, times = 1500))
dfBF = data.frame(AvgCVBF = logBFs[,1], PTBF = logBFs[,2], datasetsize = datasetsizes)
p1 = ggplot(dfBF, aes(x = datasetsize)) + geom_point(aes(y = AvgCVBF)) + geom_hline(yintercept = -log(20), color = "blue") + labs(x = "n", y = "log(CVBF)")
ggsave("logCVBFnullvsn.pdf", plot = p1, device = "pdf")

p1 = ggplot(dfBF, aes(x = datasetsize)) + geom_point(aes(y = PTBF)) + geom_hline(yintercept = -log(20), color = "blue") + labs(x = "n", y = "log(BF)")
ggsave("logPTBFnullvsn.pdf", plot = p1, device = "pdf")




smallestrepcase = subset(dfBF, datasetsize == 200)

medrepcase = subset(dfBF, datasetsize == 400)

largestrepcase = subset(dfBF,datasetsize == 800)

sum(smallestrepcase$AvgCVBF > 0)
#2
sd(smallestrepcase$AvgCVBF)
#1.60
sd(smallestrepcase$PTBF)
#2.05

sum(medrepcase$AvgCVBF < -log(20))
#[1] 1498
sd(medrepcase$AvgCVBF)
#1.95
sd(medrepcase$PTBF)
#2.52
median(medrepcase$PTBF)
#-4.06

sum(medrepcase$AvgCVBF > -log(20))
#[1] 2

sum(smallestrepcase$AvgCVBF > -log(20))
#46

max(largestrepcase$AvgCVBF)
#-7.82

sd(largestrepcase$AvgCVBF)
#2.41
sd(largestrepcase$PTBF)
#3.31
sum(largestrepcase$PTBF > 0) / length(largestrepcase$PTBF)
#Exactly .05

# ggplot(test_data, aes(date)) + 
#   geom_line(aes(y = var0, colour = "var0")) + 
#   geom_line(aes(y = var1, colour = "var1"))