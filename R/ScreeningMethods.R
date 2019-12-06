#' Screen a data set for important functions in parallel
#'
#' @param datasetX A matrix containing values that are predictors for the Y values
#' @param datasetY A vector containg the class that each predictor corresponds to. For now can only handle binary responses.
#' @param method A string containing the type of screening to do. Can be "SIS", "KS", "CVBF" or "PT"
#' @param ncores A integer that corresponds to the number of cores to be used for parallelizing computation
#' @param cutoff A real number that corresponds either to an alpha value for testing or a cutoff value on how large the Bayes factor needs to be to conclude a difference exists.
#' @param train1ids A vector of ids that correspond to which observations to use for the training set for the first data set
#' @param trainsize1 Size of the training set for one of the classes for CVBF
#' @param trainsize2 Size of the training set for the other one of the classes for CVBF
#' @param train2ids A vector of ids that correspond to which observations to use for the training set for the second data set
#' @param seed A seed for CVBF based screening, can use this to reproduce results instead of train_ids.
#' @param Ginv A function to compute quantiles with for Polya tree.
#' @param c Tuning parameter for Polya tree, signifies how impactful prior should be. 
#' @param leveltot Depth of Polya tree to construct if Polya tree based screening is type of screening chosen
#' @param PTscale A True / false variable. Should columns be standardized before proceeding with Polya tree based screening? Default is to screen as recommended by authors.
#'
#' @return A list of variables that are interpreted to be important
#' @export
#'
#' @examples
#' data(gisettetrainlabs)
#' data(gisettetrainpreds)
#' nworkers = detectCores()
#' ImpVarsSIS1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "SIS", ncores = nworkers / 2)
#' length(ImpVarsSIS1$varspicked)
#' ImpVarsKS1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "KS", ncores = nworkers / 2)
#' length(ImpVarsKS1$varspicked)
#' ImpVarsPT1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "PT", ncores = nworkers / 2, c = 1, leveltot = 12, Ginv = qnorm, PTscale = TRUE)
#' #Only do on first 500
#' length(ImpVarsPT1$varspicked)
#' hist(ImpVarsPT1$logBFlist)
#' ImpVarsCVBF1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "CVBF", ncores = nworkers / 2, trainsize1 = 2960, trainsize2 = 2960, seed = 200)
#' length(ImpVarsCVBF1$varspicked)
#' 
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
    ImportantVars = list(pvaluelist = pvaluelist , varspicked = which(pvaluelist < cutoff))
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
    ImportantVars = list(pvaluelist = pvaluelist , varspicked = which(pvaluelist < cutoff))
  } else if(method == "CVBF"){
    if(is.null(cutoff))
    {
      cutoff = 0
    }
    cl <- makeCluster(ncores) 
    registerDoParallel(cl)
    logBFlist <- foreach(j=1:p, .combine=c, .export = c("CVBFtestrsplit", "HallKernel", "laplace.kernH2c", "KHall", "logpriorused", "logintegrand.Hall", "loglike.KHall", "integrand.Hall"), .packages = "rootSolve") %dopar% {
      tempMatrix = CVBFtestrsplit(datasetX[Class0ind, j], datasetX[-Class0ind, j], trainsize1, trainsize2, seed = seed, train1_ids = train1ids, train2_ids = train2ids)$logBF 
      
      tempMatrix #Equivalent to logBFlist = c(logBFlist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = list(logBFlist = logBFlist, varspicked = which(logBFlist > cutoff))
    
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
      tempMatrix = PolyaTreetest(datasetX[Class0ind, j], datasetX[-Class0ind, j], Ginv = Ginv, c = c, leveltot = leveltot)$logBF 
      
      tempMatrix #Equivalent to logBFlist = c(logBFlist, tempMatrix)
    }
    stopCluster(cl)
    ImportantVars = list(logBFlist = logBFlist, varspicked = which(logBFlist > cutoff))
    #ImportantVars = which(logBFlist > cutoff)
  }
  else{
    stop("Method entered does not match with screening methods supported. Please recheck spelling or pick an option from 'SIS', 'KS', 'CVBF', or 'PT'")
  }
  return(ImportantVars)
}
