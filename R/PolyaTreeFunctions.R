#' Polya Tree test
#'
#' @param datasetX A set of data (or in the case of screening a predictor corresponding to one class), one of the data sets we want to check if its distribution is the same as dataset Y's. 
#' @param datasetY Another set of data (or in the case of screening a predictor corresponding to another class), the other data set we're comparing to datasetX
#' @param Ginv A function that can compute quantiles of some distribution. The default is qnorm. Can specify another function, needs to be able to take in a number between 0 and 1 and return back some positive value. It needs to be a quantile function.
#' @param c A tuning parameter corresponding to how influential the prior should be. Authors recommend to set to 1. Can change from 1. The larger the tuning parameter the more influential the prior. The smaller the tuning parameter the less influential the prior.
#' @param leveltot Total number of levels deep the tree should go. 9 is given as a default. Some authors recommend going to $log_2(sample size)$, but doesn't need to be done. The deeper the tree the more computation that is required. 
#'
#' @return Currently returns a scalar that corresponds to the log BF of the computed test. 
#' @export
#'
#' @examples
PolyaTreetest <- function(datasetX, datasetY, Ginv = qnorm, c = 1, leveltot = 9)  
{
  #A lot of places can be sped up I think by predefining epsilon list and epsilonlist2.
  #The loops themselves also take some time, maybe we can throw the whole thing into C, but the problem is that
  #We have a jagged matrix that's essentially run with a list rn, changing that to be in C might be trickier than it sounds.
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

PolyaTreetest <- function(datasetX, datasetY, Ginv = qnorm, c = 1, leveltot = 9)  
{
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