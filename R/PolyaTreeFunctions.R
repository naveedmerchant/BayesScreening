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
PolyaTreetest <- function(datasetX, datasetY, Ginv = qnorm, c = NULL, leveltot = NULL)  
{
  #A lot of places can be sped up I think by predefining epsilon list and epsilonlist2.
  #The loops themselves also take some time, maybe we can throw the whole thing into C, but the problem is that
  #We have a jagged matrix that's essentially run with a list rn, changing that to be in C might be trickier than it sounds.
  #But it is something I wanted to look into
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
    for(D in 1:(length(alphalist[[j]]) - 1))
    {
      bj[j] = bj[j] + lbeta(alphalist[[j]][D],alphalist[[j]][D+1]) +
        lbeta(alphalist[[j]][D] + splitlist[[j]][D] +
                splitlist2[[j]][D],alphalist[[j]][D+1] +
                splitlist[[j]][D+1] + splitlist2[[j]][D+1] ) -
        lbeta(alphalist[[j]][D] + splitlist[[j]][D], alphalist[[j]][D+1] + splitlist[[j]][D+1]) -
        lbeta(alphalist[[j]][D] + splitlist2[[j]][D], alphalist[[j]][D+1] + splitlist2[[j]][D+1])
      
    }
  }
  return(sum(bj))
  
}

#' Title
#'
#' @param datasetX A dataset to compute the Polya Tree prior on
#' @param Ginv A quantile function of some distribution, use to make bins
#' @param c A scalar. The higher this is, the more influential the prior is on the data set.
#' @param leveltot The number of levels the Polya Tree should go down.
#'
#' @return A list of vectors called alphalist, a list of vectors called splitlist, c, leveltot, and Ginv. Alphalist is used to construct a prior. Splitlist is used to construct a likelihood. Call the collection of these a Polya Tree object.
#' @export
#'
#' @examples
PolyaTreePriorLikCons <- function(datasetX, Ginv = qnorm, c = 1, leveltot = 9)  
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
  return(list(alphalist = alphalist, splitlist = splitlist, c = c, leveltot = leveltot, Ginv = Ginv))
}

#' Title
#'
#' @param PolyaTreePriorLik1 An object constructed by PolyaTreePriorLikCons for a dataset. 
#' @param PolyaTreePriorLik2 Another object constructed by PolyaTreePriorLikCons for another dataset.
#'
#' @return A scalar, which corresponds to the log BF of the test. Another vector is given which corresponds to the contribution of the log BF at each level.
#' @export
#'
#' @examples
PolyaTreeBFcons <- function(PolyaTreePriorLik1, PolyaTreePriorLik2)
{
  if(PolyaTreePriorLik1$leveltot != PolyaTreePriorLik2$leveltot)
  {
    stop("Polya Trees go down to different depths, this test is inappropriate")
  }
  else{
    leveltot = PolyaTreePriorLik1$leveltot
  }
  
  if(PolyaTreePriorLik1$c != PolyaTreePriorLik2$c)
  {
    stop("Polya Trees are generated with different prior, test is possible, but BF is meaningless")
  }
  else{
    c = PolyaTreePriorLik1$c
  }
  #browser()
  bj = rep(0, times = leveltot)
  for(j in 1:leveltot)
  {
    for(D in 1:(length(PolyaTreePriorLik1$alphalist[[j]]) - 1))
    {
      bj[j] = bj[j] + lbeta(PolyaTreePriorLik1$alphalist[[j]][D],PolyaTreePriorLik1$alphalist[[j]][D+1]) +
        lbeta(PolyaTreePriorLik1$alphalist[[j]][D] + PolyaTreePriorLik1$splitlist[[j]][D] +
                PolyaTreePriorLik2$splitlist[[j]][D],PolyaTreePriorLik1$alphalist[[j]][D+1] +
                PolyaTreePriorLik1$splitlist[[j]][D+1] + PolyaTreePriorLik2$splitlist[[j]][D+1] ) -
        lbeta(PolyaTreePriorLik1$alphalist[[j]][D] + PolyaTreePriorLik1$splitlist[[j]][D], PolyaTreePriorLik1$alphalist[[j]][D+1] + PolyaTreePriorLik1$splitlist[[j]][D+1]) -
        lbeta(PolyaTreePriorLik1$alphalist[[j]][D] + PolyaTreePriorLik2$splitlist[[j]][D], PolyaTreePriorLik1$alphalist[[j]][D+1] + PolyaTreePriorLik2$splitlist[[j]][D+1])

    }
  }
  return(list(logBF = sum(bj), logBFcont = bj))
}
 

  

#' Title
#'
#' @param PolyaTreePriorLik An object constructed by PolyaTreePriorLikCons for a dataset. 
#' @param ndraw Number of draws desired from Predictive Polya Tree distribution.
#'
#' @return A list of draws from the Predictive Polya Tree Posterior  
#' @export
#'
#' @examples
PolyaTreePredDraws <- function(PolyaTreePriorLik, ndraw = 2000)
{
  if(ndraw <= 0)
  {
    stop("ndraw should be a positive integer that corresponds to the number of draws you want from the predictive posterior.")
  }
  
  postbetalist = PolyaTreePriorLik$alphalist
  largeqnormlist = seq(from = 0, to = 1, by = .0000001)
  largeqnormlist = PolyaTreePriorLik$Ginv(largeqnormlist)
  
  #This is actually probably inappropriate for some choices of Ginv
  #BUT that is a research topic for ANOTHER paper to address!
  #I'm not actually sure if a package exists that actually constructs and draws from these currently
  #We construct a large list at the beginning to get decent samples by drawing from this list
  #Better to do this then to do rejection sampling an incredibly large amount of times.
  
  for(j in 1:PolyaTreePriorLik$leveltot)
  {
    postbetalist[[j]] = PolyaTreePriorLik$splitlist[[j]] + PolyaTreePriorLik$alphalist[[j]]
    #Leveltot shouldn't actually be that large, so this should be pretty fast
    #Posterior of beta binomial is a beta with above coefficients
  }
  
  posteriorPTdraw = rep(NULL, times = ndraw)
  k = 1
  epsilonlist = rep(0, times = PolyaTreePriorLik$leveltot)
  for(i in 1:ndraw)
  {
    #This might be possible to speed up and avoid running in a for loop. None of the iterations depend on the prev one.
    for(m in 1:PolyaTreePriorLik$leveltot)
    {
      postdraw = rbinom(1, size = 1, prob = postbetalist[[m]][k] / (postbetalist[[m]][k+1]+ postbetalist[[m]][k]))
      if(postdraw == 1)
      {
        epsilonlist[m] = 0
      }
      else
      {
        epsilonlist[m] = 1
      }
      multivec = 2^seq(m, 1, -1)
      ind = sum(epsilonlist[1:m]*multivec) + 1
      k = ind
    }
    #epsilonlist2[[j]] <- epsilonlist
    #need to draw from qnorm((ind - 1) / 2^9)  qnorm(ind / 2^9)
    if(epsilonlist[PolyaTreePriorLik$leveltot] == 0)
    {
      posteriorPTdraw[i] = sample(largeqnormlist[(largeqnormlist < PolyaTreePriorLik$Ginv(ind / 2^(PolyaTreePriorLik$leveltot+1))) & (largeqnormlist > qnorm((ind-1) / 2^(PolyaTreePriorLik$leveltot+1)))], size = 1)
    }
    else
    {
      posteriorPTdraw[i] = sample(largeqnormlist[(largeqnormlist < PolyaTreePriorLik$Ginv((ind+1) / 2^(PolyaTreePriorLik$leveltot+1))) & (largeqnormlist > qnorm((ind) / 2^(PolyaTreePriorLik$leveltot+1)))], size = 1)
    }
    epsilonlist = c()
    #indlist[i] = ind
    k=1
    #print(i)
  }
  return(posteriorPTdraw)
  
}




  


