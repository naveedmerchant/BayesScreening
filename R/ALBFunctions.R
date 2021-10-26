#' Generate draws from a mixture of two normal distributions
#'
#' @param n The number of draws from the mixture of normals
#' @param mean1 The mean of the first normal distribution drawn
#' @param sd1 The standard deviation of the first normal distribution drawn
#' @param mean2 The mean of the second normal distribution drawn
#' @param sd2 The standard deviation of the second normal distribution drawn
#' @param mixingparam The mixing parameter, draws from the first normal with this probability
#'
#' @return Draws from a mixture of two normal distributions
#' @export
#'
#' @examples
#' mixnormdraw(100, mean1 = 1, sd1 = 1 , mean2 = 5, sd2 = 3, mixingparam = .5)
mixnormdraw <- function(n, mean1, sd1, mean2, sd2, mixingparam){
  draws = c()
  for(j in 1:n)
  {
    if(runif(1) < mixingparam)
    {
      draws[j] = rnorm(1, mean = mean1, sd = sd1)
    } else{
      draws[j] = rnorm(1, mean = mean2, sd = sd2)
    }
  }
  return(draws)
}

#' Compute the leave one out likelihood cross validation statistic
#'
#' @param h The bandwidth used to compute the likelihood cross validation statistic 
#' @param x The sample data set used on 
#'
#' @return The leave one out likelihood cross validation statistic, it is the log of the product of
#' kernel density estimates trained on all estimates except 1 and evaluated on the omitted observation
#' @export 
#'
#' @examples
#' LCV_h(h = .01, x = rnorm(100))
LCV_h = function(h,x)
{
  sum = 0
  for(i in 1:length(x))
  {
    sum = sum + log(HallKernel(h,x[-i],x[i]))
  }
  return(sum)
}



#' Compute a ALB for testing if two data sets share the same density with the Hall Kernel
#'
#' @param x A data set
#' @param y Another data set
#' @param bw (optional) A bandwidth to use to compute the Kernel density estimate of the joint data set, if not provided
#' will use either a normal plug in or the one that does best in terms of leave one out cross validation
#' @param objbw (optional) Leave one out likelihood cross validation statistic on the joint data set
#' @param option (optional) If left alone or set to 1, uses likelihood cross validation to compute the ideal bandwidth, 
#' if set to 2 uses the normal plug-in instead
#'
#' @return A list of elements, BF refers to the ALB, bw refers to the bandwidth used for computing the ALB, objbw refers to
#' the value of the leave one out likelihood cross validation statistic on the joint data set 
#' @export
#'
#' @examples
#' set.seed(100)
#' dataset1 = rnorm(200)
#' dataset2 = rnorm(200) 
#' ALB1 = leaveoneoutBFLR2(dataset1, dataset2)
#' ALB1$BF #Gives back the ALB value
leaveoneoutBFLR2 = function(x, y, bw = NULL, objbw = NULL, option = NULL)
{
  if(is.null(bw))
  {
    if(is.null(option))
    {
      bwinfoJ = optimise(f =  LCV_h, x = c(x,y), lower = .001, upper = 5, maximum = TRUE)
      bw = bwinfoJ$maximum
      objbw = bwinfoJ$objective
    }
    if(option == "1")
    {
      bwinfoJ = optimise(f =  LCV_h, x = c(x,y), lower = .001, upper = 5, maximum = TRUE)
      bw = bwinfoJ$maximum
      objbw = bwinfoJ$objective
    }
    if(option == "2")
    {
      n = length(x) + length(y)
      bw = 0.162 * n ^ (-1 / 5) * sd(c(x, y))
      objbw = LCV_h(bw, c(x,y))
    }
  }
  if(is.null(objbw))
  {
    objbw = LCV_h(bw, c(x,y))
  }
  return(list("BF" = (LCV_h(bw, x) + LCV_h(bw, y) - objbw) / (length(x) + length(y)), "bw" = bw, "objbw" = objbw))
}


#' ALB Variable Screening for binary classification
#'
#' @param class1trainpreds Predictors of class 1 
#' @param class2trainpreds Predictors of class 2
#' @param testingset A testing set for evaluating the accuracy of a simple independent classifier
#' @param ALBbwlist A list of bws that were used to compute the ALB, if not given, we will compute 
#' @param ALBlist A list of ALBs, if not given, we will compute
#' @param npermvar Number of permutations to do to compute ALB under null distribution
#' @param alpha An alpha value that corresponds to a type I error rate of each test, set to .05, this explode the familywise
#' type I error rate to be higher, but that's okay for a screening problem.
#'
#' @return A list of variables
#' @export
#'
#' @examples
IndepNWEstprobNpermC = function(class1trainpreds, class2trainpreds, testingset, ALBbwlist = NULL, ALBlist = NULL, npermvar = 3, alpha = .05){
  p = ncol(class1trainpreds)
  ntest = nrow(testingset)
  ntrain = nrow(class1trainpreds) + nrow(class2trainpreds)
  if(is.null(ALBlist))
  {
    ALBlist = rep(0, times = p)
    ALBbwlist = rep(0, times = p)
    objbwlist = rep(0, times = p)
    for(j in 1:p)
    {
      x2 = class1trainpreds[,j]
      y2 = class2trainpreds[,j]
      leaveoneoutinfo = leaveoneoutBFLR3(x = x2, y = y2)
      ALBlist[j] = leaveoneoutinfo$BF
      ALBbwlist[j] = leaveoneoutinfo$bw
      objbwlist[j] = leaveoneoutinfo$objbw
    }
  }
  combclasspreds = rbind(class1trainpreds, class2trainpreds)
  combclassind = c(rep(0, times = nrow(class1trainpreds)), rep(1, times = nrow(class2trainpreds)))
  permbfvec = rep(0, times = npermvar * p)
  z = 1
  for(j in 1:p)
  {
    permclassinds = sample(combclassind, size = length(combclassind))
    class1trainpreds = combclasspreds[permclassinds == 1,]
    class2trainpreds = combclasspreds[permclassinds == 0,]
    for(k in 1:npermvar)
    {
      x2 = class1trainpreds[,j]
      y2 = class2trainpreds[,j]
      leaveoneoutinfo = leaveoneoutBFLR3(x = x2, y = y2, bw = ALBbwlist[j], objbw = objbwlist[j])
      permbfvec[z] = leaveoneoutinfo$BF
      z = z + 1
    }
  }
  cutoff = quantile(x = permbfvec, probs = 1 - alpha)
  
  ALBvars = which(ALBlist > cutoff)
  #testingset = testpreds
  trainingset = rbind(class1trainpreds, class2trainpreds)
  probs = c()
  for(j in 1:ntest)
  {
    prod1 = 0
    prod2 = 0
    for(k in ALBvars)
    {
      prod1 = prod1 + log(HallKernel(h = ALBbwlist[k], datagen2 = class1trainpreds[,k], x = testingset[j,k]))
      prod2 = prod2 + log(HallKernel(h = ALBbwlist[k], datagen2 = class2trainpreds[,k], x = testingset[j,k]))
    }
    prodmin = max(prod1, prod2)
    prod1 = exp(prod1 - prodmin)
    prod2 = exp(prod2 - prodmin)
    probs[j] = (prod1 * (nrow(class1trainpreds)) / ntrain) / (prod1 * (nrow(class1trainpreds)) / ntrain + prod2 * (nrow(class2trainpreds)) / ntrain)
  }
  return(list(probs = probs, numvars = length(ALBvars), ALBlist = ALBlist, ALBbwlist = ALBbwlist, cutoff = cutoff, permbfs = permbfvec))
}