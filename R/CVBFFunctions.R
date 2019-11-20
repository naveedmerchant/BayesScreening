#' Title
#'
#' @param h A bandwidth or vector of bandwidths tot ry out
#' @param y A validation set to evaluate the KDE on
#' @param x A training set to build the Hall kernel density estimate on
#'
#' @return -log likelihood evaluation, where the likelihood is constructed using the training data
#' @export
#'
#' @examples
loglike.KHall=function(h,y,x){
  
  n = length(x)
  m = length(y)
  nh = length(h)
  M = t(matrix(y,m,n))
  llike=1:nh
  
  for(j in 1:nh){
    
    M1 = (x-M)/ h[j]
    M1 = KHall(M1) / h[j]
    fhat=as.vector(M1 %*% matrix(1,m,1))/m
    fhat[fhat<10^(-320)]=10^(-320)
    llike[j]=sum(log(fhat))
  }
  return(-llike)
}

#' Title Hall Kernel
#'
#' @param x The parameter to evaluate the Hall Kernel density on 
#'
#' @return Evaluation of K(X), where $K \propto exp^{-.5log(1 + abs(x))^2}$ and is a valid PDF
#' @export
#'
#' @examples
KHall=function(x){
  
  con=sqrt(8*pi*exp(1))*pnorm(1)
  K=exp(-0.5*(log(1+abs(x)))^2)/con
  return(K)
}



#' Title
#'
#' @param h Bandwidth parameter, this is what the integral is with respect to.
#' @param y Validation set
#' @param x Training set
#' @param cons The integrand is typically too big to easily evaluate, need a constant to help evaluate it. Typically just use the log version of this function instead
#' @param hhat The parameter to center the prior at, recommended to be maximizer of the likelihood
#'
#' @return
#' @export
#'
#' @examples
integrand.Hall=function(h,y,x,cons,hhat){
  
  n=length(x)
  beta=hhat
  Prior=(2*beta / sqrt(2*pi)) * exp(-beta^2/h^2) / h^2
  arg=-loglike.KHall(h,y,x)-cons
  arg[arg>700]=700
  f=exp(arg)*Prior
  return(f)
}



#' Title
#'
#' @param h Bandwidth parameter, this is the variable that is being integrated over
#' @param y Validation set to evaluate the likelihood over
#' @param x Training set to build the KDE
#' @param hhat Parameter that specifies where prior should be centered
#'
#' @return
#' @export
#'
#' @examples
logintegrand.Hall=function(h,y,x,hhat){
  
  n=length(x)
  beta=hhat
  Prior=log(2*beta/sqrt(2*pi))-beta^2/h^2-2*log(h)
  f=loglike.KHall(h,y,x)-Prior
  return(f)
  
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


CVBFtestrsplit = function(dataset1, dataset2, trainsize1, trainsize2, seed = NULL)
{
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  train_ids = sample(1:length(dataset1), size = trainsize1)
  XT1 = dataset1[train_ids]
  XV1 = dataset1[-train_ids]
  train_ids2 = sample(1:length(dataset2), size = trainsize2)
  YT1 = dataset2[train_ids]
  YV1 = dataset2[-train_ids]
  
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
  bwlik2 = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML1 = laplace.kernH2c(y = XT1, x = XV1, hhat = bwlik2$maximum, c = bwlik2$objective)
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = YT1, x = YV1)))}
  bwlikcy = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML2 = laplace.kernH2c(y = YT1, x = YV1, hhat = bwlikcy$maximum, c= bwlikcy$objective)
  
  likveccombc = function(h) {sum(log(HallKernel(h,datagen2 = c(XT1,YT1), x = c(XV1,YV1))))}
  bwlikcombc = optimize(f = function(h){  likveccombc(h)}, lower = 0, upper = 10, maximum = TRUE)
  #ExpectedKernMLcomb = bwlikcombc$objective
  
  ExpectedKernMLcomb = laplace.kernH2c(y = c(XT1,YT1), x = c(XV1,YV1), hhat = bwlikcy$maximum, c= bwlikcombc$objective)
  
  return(logBF = ExpectedKernML1[1] + ExpectedKernML2[1] - ExpectedKernMLcomb[1])
}

