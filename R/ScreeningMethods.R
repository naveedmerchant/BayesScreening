#' Screen a data set for important functions in parallel
#'
#' @param datasetX A matrix containing values that are predictors for the Y values
#' @param datasetY A vector containg the class that each predictor corresponds to. For now can only handle binary responses.
#' @param method A string containing the type of screening to do. Can be "SIS", "KS", "CVBF" or "PT"
#' @param ncores A integer that corresponds to the number of cores to be used for parallelizing computation
#' @param cutoff A real number that corresponds either to an alpha value for testing or a cutoff value on how large the Bayes factor needs to be to conclude a difference exists.
#'
#' @return A list of variables that are interpreted to be important
#' @export
#'
#' @examples
ParScreenVars <- function(datasetX, datasetY, method = "SIS", ncores = 1, cutoff = NULL, train1ids = NULL, trainsize1 = NULL, trainsize2 = NULL, train2ids = NULL, seed = NULL, Ginv = NULL, c = NULL, leveltot = NULL, PTscale = TRUE)
{
  n = nrow(datasetX)
  p = ncol(datasetX)
  if(n != length(datasetY))
  {
    stop("Dimensions of datasetX and datasetY do not match!")
  }
  if(ncores > detectCores())
  {
    stop("ncores exceeds number of available cores on computer!")
  }
  Yfactors = unique(datasetY)
  Class0ind = which(datasetY == Yfactors[1])
  if(length(Yfactors) != 2)
  {
    stop("Program currently only supports screening for binary classification.")
  }
  if(method == "SIS")
  {
    if(is.null(cutoff))
    {
      cutoff = .005
    }
    cl <- makeCluster(ncores) 
    registerDoParallel(cl)
    pvaluelist <- foreach(j=1:p, .combine=c) %dopar% {
      tempMatrix = t.test(datasetX[Class0ind,j], datasetX[-Class0ind,j])$p.value 
      
      
      tempMatrix #Equivalent to pvaluelist = c(pvaluelist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = which(pvaluelist < cutoff)
  } else if(method == "KS"){
    if(is.null(cutoff))
    {
      cutoff = .005
    }
    cl <- makeCluster(ncores) 
    registerDoParallel(cl)
    pvaluelist <- foreach(j=1:p, .combine=c) %dopar% {
      tempMatrix = ks.test(datasetX[Class0ind,j], datasetX[-Class0ind,j])$p.value 
      
      tempMatrix #Equivalent to pvaluelist = c(pvaluelist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = which(pvaluelist < cutoff)
  } else if(method == "CVBF"){
    if(is.null(cutoff))
    {
      cutoff = 0
    }
    if(is.null(trainsize1))
    {
      stop("Please enter a trainsize for the dataset of class1.")
    }
    if(is.null(trainsize2))
    {
      stop("Please enter a trainsize for the dataset of class2.")
    }
    cl <- makeCluster(ncores) 
    registerDoParallel(cl)
    logBFlist <- foreach(j=1:p, .combine=c, .export = c("CVBFtestrsplit", "HallKernel", "laplace.kernH2c", "KHall", "logpriorused", "logintegrand.Hall", "loglike.KHall", "integrand.Hall"), .packages = "rootSolve") %dopar% {
      tempMatrix = CVBFtestrsplit(datasetX[Class0ind, j], datasetX[-Class0ind, j], trainsize1, trainsize2, seed = seed, train1_ids = train1ids, train2_ids = train2ids)$logBF 
      
      tempMatrix #Equivalent to logBFlist = c(logBFlist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = which(logBFlist > cutoff)
    
  } else if(method == "PT"){
    if(is.null(cutoff))
    {
      cutoff = 0
    }
    if(is.null(Ginv))
    {
      Ginv = qnorm
    }
    if(PTscale == TRUE)
    {
      datasetX = scale(datasetX, center = TRUE, scale = TRUE)
    }
    #browser()
    cl <- makeCluster(ncores) 
    registerDoParallel(cl)
    logBFlist <- foreach(j=1:p, .combine=c, .export = c("PolyaTreetest")) %dopar% {
      tempMatrix = PolyaTreetest(datasetX[Class0ind, j], datasetX[-Class0ind, j], Ginv = Ginv, c = c, leveltot = leveltot) 
      
      tempMatrix #Equivalent to logBFlist = c(logBFlist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = which(logBFlist > cutoff)
  }
  else{
    stop("Method entered does not match with screening methods supported. Please recheck spelling or pick an option from 'SIS', 'KS', 'CVBF', or 'PT'")
  }
  return(ImportantVars)
}
