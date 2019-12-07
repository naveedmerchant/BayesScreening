#' Evaluate the log prior on one of the data sets. Called by other functions.
#'
#' @param h A value to evaluate the prior on
#' @param hhat Tuning parameter for the prior
#'
#' @return Evaluates a particular type of prior placed on the bandwidth for CVBF.
#' @export
#'
#' @examples
#' logpriorused(.1,.1)
logpriorused <- function(h,hhat)
{
  beta = hhat
  Prior = log(2*beta) - .5*log(pi) - 2*log(h) - (beta^2 / h^2)
  return(Prior)
}


#' Evaluate a log likelihood using the Hall Kernel
#'
#' @param h A bandwidth or vector of bandwidths tot ry out
#' @param y A validation set to evaluate the KDE on
#' @param x A training set to build the Hall kernel density estimate on
#'
#' @return -log likelihood evaluation, where the likelihood is constructed using the training data
#' @export
#'
#' @examples
#' dataset1 = rnorm(100)
#' DT = dataset1[1:50]
#' DV = dataset1[51:100]
#' loglike.Khall(.01, DT, DV)
loglike.KHall=function(h, y, x){
  
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

#' Title Hall Kernel evaluation
#'
#' @param x The parameter to evaluate the Hall Kernel density on 
#'
#' @return Evaluation of K(X), where K is a heavy tailed kernel. For more information on K(x), see
#'  K_0 on Hall's paper on Kullback Leibler loss (Annals of Statistics 1957)
#' @export
#'
#' @examples
#' Khall(.1)
KHall=function(x){
  
  con = sqrt(8 * pi * exp(1)) * pnorm(1)
  K=exp(-0.5 * (log(1 + abs(x)))^2) / con
  return(K)
}

#' Compute log of the integrand of the Marginal likelihood for CVBF
#'
#' @param h Bandwidth parameter, this is the variable that is being integrated over
#' @param y Validation set to evaluate the likelihood over
#' @param x Training set to build the KDE
#' @param hhat Parameter that specifies where prior should be centered
#'
#' @return Evaluation of the integrand at a particular log likelihood value.
#' @export
#'
#' @examples
#' dataset1 = rnorm(100)
#' DT = dataset1[1:50]
#' DV = dataset1[51:100]
#' logintegrand.Hall(.01, DT, DV, .1)
logintegrand.Hall=function(h, y, x, hhat){
  
  n = length(x)
  beta = hhat
  Prior = log(2 * beta / sqrt(2 * pi)) - beta^2 / h^2 - 2 * log(h)
  f=loglike.KHall(h, y, x) - Prior
  return(f)
  
}


#' Compute Marginal Likelihoods for CVBF
#'
#' @param y Validation Set
#' @param x Training set
#' @param hhat Bandwidth parameter that maximizes the log likelihood
#' @param c A constant that is equal to the log likelihood + log prior evaluated at the maximum
#'
#' @return Evaluation of the CVBF marginal likelihood via Laplace Approximation. Also returns bandwidth that maximized log integrand, and hessian of log integrand at maximum.
#' @export
#'
#' @examples
#' dataset1 = rnorm(100)
#' DT = dataset1[1:50]
#' DV = dataset1[51:100]
#' likvec = function(h) {sum(log(HallKernel(h,datagen2 = DT, x = DT)))}
#' bwlikcy = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
#' ExpectedKernML2 = laplace.kernH2c(y = DT, x = DV, hhat = bwlikcy$maximum, c= bwlikcy$objective + logpriorused(h = bwlikcy$maximum, hhat = bwlikcy$maximum))
#' ExpectedKernML2

laplace.kernH2c = function(y, x, hhat, c){
  
  n = length(x)
  out = hessian(f = function(h) {logintegrand.Hall(h, y=y, x=x, hhat = hhat)}, x = hhat, pert = 10^(-7))
  cons = c
  hess= abs(out)
  laplace = 0.5 * log(2 * pi) - 0.5 * log(hess) + cons
  c(laplace,hhat,hess)
  
}


#' Compute a CVBF, a Bayes factor that checks if two data sets share the same distribution.
#'
#' @param dataset1 One dataset that we want to check if it has the same distribution as another data set
#' @param dataset2 Another dataset that we want to check if it has the same distribution as another data set
#' @param trainsize1 The training set size of dataset 1
#' @param trainsize2 The training set size of dataset 2
#' @param train1_ids Indices for the training set of dataset 1
#' @param train2_ids Indices for the training set of dataset 2
#' @param seed The seed used to generate training set and validation sets for both of the data sets
#'
#' @return A list containing a log BF that tests whether two distributions are the same via CVBF, and the training sets used to generate the CVBF.
#' @export
#'
#' @examples
#' set.seed(100)
#' dataset1 = rnorm(200)
#' dataset2 = rnorm(200) 
#' CVBF1 = CVBFtestrsplit(dataset1, dataset2, trainsize1 = 100, trainsize2 = 100)
#' CVBF1$logBF #Gives back the log Bayes factor
#' CVBF1$train1_ids #Gives back the training set of the first data set
#' CVBF1$train2_ids #Gives back the training set of the second data set
#' 
CVBFtestrsplit = function(dataset1, dataset2, trainsize1, trainsize2, seed = NULL, train1_ids = NULL, train2_ids = NULL)
{
  if(is.null(trainsize1))
  {
    stop("Please enter a trainsize for the dataset of class1.")
  }
  if(is.null(trainsize2))
  {
    stop("Please enter a trainsize for the dataset of class2.")
  }
  if(!is.null(seed))
  {
    set.seed(seed)
  }
  if(!is.null(train1_ids))
  {
    train_ids = train1_ids
  } else{
    train_ids = sample(1:length(dataset1), size = trainsize1)
  }
  if(!is.null(train2_ids))
  {
    train_ids2 = train2_ids
  } else{
    train_ids2 = sample(1:length(dataset2), size = trainsize2)
  }
  
  XT1 = dataset1[train_ids]
  XV1 = dataset1[-train_ids]
  
  YT1 = dataset2[train_ids2]
  YV1 = dataset2[-train_ids2]
  
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
  bwlik2 = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML1 = laplace.kernH2c(y = XT1, x = XV1, hhat = bwlik2$maximum, c = bwlik2$objective + logpriorused(h = bwlik2$maximum, hhat = bwlik2$maximum))
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = YT1, x = YV1)))}
  bwlikcy = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML2 = laplace.kernH2c(y = YT1, x = YV1, hhat = bwlikcy$maximum, c= bwlikcy$objective + logpriorused(h = bwlikcy$maximum, hhat = bwlikcy$maximum))
  
  likveccombc = function(h) {sum(log(HallKernel(h,datagen2 = c(XT1,YT1), x = c(XV1,YV1))))}
  bwlikcombc = optimize(f = function(h){  likveccombc(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernMLcomb = laplace.kernH2c(y = c(XT1,YT1), x = c(XV1,YV1), hhat = bwlikcombc$maximum, c= bwlikcombc$objective + logpriorused(h = bwlikcombc$maximum, hhat = bwlikcombc$maximum))
  
  return(list(logBF = ExpectedKernML1[1] + ExpectedKernML2[1] - ExpectedKernMLcomb[1], train1_ids = train_ids, train2_ids = train_ids2))
}


#' Compute Hall Kernel density estimate, given a training set and a set of values to evaluate on
#'
#' @param h Bandwidth parameter
#' @param datagen2 Training set for Hall KDE
#' @param x Value (or values) to evaluate KDE on. Can be a scalar, a matrix, or a vector.
#'
#' @return An object of the same dimension as x that evaluates a KDE trained on datagen2 on everyone of the points in x
#' @export
#'
#' @examples
#' 
HallKernel = function(h,datagen2,x)
{
  sum = 0
  for(i in 1:length(datagen2))
  {
    sum = sum + (((8*pi*exp(1))^.5)*pnorm(1))^(-1)*exp(-.5*(log(1+abs(x - datagen2[i])/h)^2))
  }
  return((1/(length(datagen2) * h)) * sum)
}



#' Draw bandwidths from CVBF predictive posterior by Metropolis Hasting sampling
#'
#' @param ndraw Number of unique draws desired for the bandwidth parameter from the posterior.
#' @param propsd A tuning parameter, corresponds to what proposal standard deviation should be for when using MH to traverse the posterior. Should be chosen with care to ensure good mixing.
#' @param maxIter The max number of MH iterations to try. Do not set to be too large. It will kick the code out if acceptance rates for MH are small.
#' @param XT1 Training set for a data set
#' @param XV1 Validation set for a data set
#' @param startingbw A value to start the MH chain at. If not provided, starts at posterior mode.
#'
#' @return A list of bandwidths that come from the posterior distribution. This will be larger than ndraw, as some draws will be repeats.
#' @export
#'
#' @examples
#' set.seed(500)
#' datasetsample1 = rnorm(600)
#' trainingindices1 = sample(1:600, size = 300)
#' XT1 = datasetsample1[trainingindices1]
#' XV1 = datasetsample1[-trainingindices1]
#' predbwvec1 = PredCVBFMHbw(ndraw = 500, maxIter = 5000, XT1 = XT1, XV1 = XV1)

PredCVBFMHbw = function(ndraw = 100, propsd = NULL, maxIter = 10000, XT1, XV1, startingbw = NULL)
{
  bwvec = c()
  if(is.null(startingbw))
  {
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
    bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    startingbw = bwvec[1] = bwlik$maximum
  }
  else{
    bwvec[1] = startingbw
  }
  if(is.null(propsd))
  {
    propsd = startingbw * length(XT1) ^ (-1 / 10)
  }
  beta = startingbw
  postcurr = sum(log(HallKernel(bwvec[1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[1]) - (beta^2 / bwvec[1]^2)
  
  acceptances = 0
  for(j in 1:maxIter)
  {
    bwprop = rnorm(1, mean = bwvec[j], sd = propsd)
    if(bwprop <= 0)
    {
      #bandwidths cannot be negative, ignore sample that was drawn.
      bwvec[j+1] = bwvec[j]
    }
    else
    {
      #MH procedure for drawing bandwidths
      postprop = sum(log(HallKernel(bwprop, datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwprop) - (beta^2 / bwprop^2)
      logpostdif = postprop - postcurr
      if(exp(logpostdif) > runif(1))
      {
        bwvec[j+1] = bwprop
        postcurr = sum(log(HallKernel(bwvec[j+1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[j+1]) - (beta^2 / bwvec[j+1]^2)
        acceptances = acceptances + 1
      }
      else
      {
        bwvec[j+1] = bwvec[j]
      }
    }
    if(acceptances >= ndraw)
    {
      break
    }
  }
  if(acceptances < ndraw)
  {
    warning(paste("Terminated because max iterations reached. The number of acceptances is", acceptances, "consider changing max iter or propsd"))
  }
  return(list(predbwsamp = bwvec, acceptancetot = acceptances, drawtot = j))
}


#' Draw bandwidths from CVBF predictive posterior by independent Metropolis Hasting sampling
#'
#' @param ndraw Number of unique draws desired for the bandwidth parameter from the posterior.
#' @param propsd A tuning parameter, corresponds to what proposal standard deviation should be for when using MH to traverse the posterior. We give a decent theoretical default. May need to be altered if performance is bad.
#' @param maxIter The max number of MH iterations to try. Do not set to be too large. It will kick the code out if acceptance rates for MH are small.
#' @param XT1 Training set for a data set
#' @param XV1 Validation set for a data set
#' @param startingbw A value to start the MH chain at. If not provided, starts at posterior mode. All proposals will be drawn from a distribution whose center is startingbw. This is normally a bad idea, but the posterior is some type of unimodal distribution, so this is actually effective.
#'
#' @return A list of bandwidths that come from the posterior distribution. This will be larger than ndraw, as some draws will be repeats.
#' @export
#'
#' @examples
#' set.seed(500)
#' datasetsample1 = rnorm(600)
#' trainingindices1 = sample(1:600, size = 300)
#' XT1 = datasetsample1[trainingindices1]
#' XV1 = datasetsample1[-trainingindices1]
#' predbwvec1 = PredCVBFIndepMHbw(ndraw = 500, maxIter = 5000, XT1 = XT1, XV1 = XV1)
#'
 
PredCVBFIndepMHbw = function(ndraw = 100, propsd = NULL, maxIter = 10000, XT1, XV1, startingbw = NULL)
{
  bwvec = c()
  if(is.null(startingbw))
  {
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
    bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    startingbw = bwvec[1] = bwlik$maximum
  }
  else{
    bwvec[1] = startingbw
  }
  if(is.null(propsd))
  {
    propsd = startingbw * length(XT1) ^ (-1 / 10)
  }
  beta = startingbw
  postcurr = sum(log(HallKernel(bwvec[1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[1]) - (beta^2 / bwvec[1]^2)
  
  acceptances = 0
  for(j in 1:maxIter)
  {
    bwprop = rnorm(1, mean = startingbw, sd = propsd)
    if(bwprop <= 0)
    {
      #bandwidths cannot be negative, ignore sample that was drawn.
      bwvec[j+1] = bwvec[j]
    }
    else
    {
      #MH procedure for drawing bandwidths
      postprop = sum(log(HallKernel(bwprop, datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwprop) - (beta^2 / bwprop^2)
      logpostdif = postprop - postcurr
      if(exp(logpostdif) > runif(1))
      {
        bwvec[j+1] = bwprop
        postcurr = sum(log(HallKernel(bwvec[j+1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[j+1]) - (beta^2 / bwvec[j+1]^2)
        acceptances = acceptances + 1
      }
      else
      {
        bwvec[j+1] = bwvec[j]
      }
    }
    if(acceptances >= ndraw)
    {
      break
    }
  }
  if(acceptances < ndraw)
  {
    warning(paste("Terminated because max iterations reached. The number of acceptances is", acceptances, "consider changing max iter or propsd"))
  }
  return(list(predbwsamp = bwvec, acceptancetot = acceptances, drawtot = j, propmean = startingbw, propsd = propsd))
}

#' Compute a Predictive Posterior function given a sequence of bandwidths from the predictive posterior
#'
#' @param bwvec A vector of bandwidths, can come from either of the methods that draw bandwidths from the posterior.
#' @param XT1 The training set
#'
#' @return A function that evaluates the predictive posterior at particular values
#' @export
#'
#' @examples
#' set.seed(500)
#' datasetsample1 = rnorm(600)
#' trainingindices1 = sample(1:600, size = 300)
#' XT1 = datasetsample1[trainingindices1]
#' XV1 = datasetsample1[-trainingindices1]
#' predbwvec1 = PredCVBFIndepMHbw(ndraw = 500, maxIter = 5000, XT1 = XT1, XV1 = XV1)
#' predpost = PredCVBFDens(predbwvec1, XT1)
#' plot(seq(from = -3, to = 3, by = .1), Predpost(predbwvec1, XT1)(seq(from = -3, to = 3, by = .1)))
#' 
PredCVBFDens = function(bwvec, XT1)
{
  PredictedDens = function(support)
  {
    predpostavg = 0
    for(j in 1:length(bwvec))
    {
      predpostavg = (1/length(bwvec))*HallKernel(bwvec[j], datagen = XT1, support) + predpostavg
    }
    return(predpostavg)
  }
  return(PredictedDens)
}










