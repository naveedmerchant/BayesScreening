#' Screen a data set for important functions in parallel
#'
#' @param datasetX A matrix containing values that are predictors for the Y values
#' @param datasetY A vector containg the class that each predictor corresponds to
#'
#' @return A list of variables that are interpreted to be important
#' @export
#'
#' @examples
ParScreenVars <- function(datasetX, datasetY, method = "SIS")
{
  n = nrow(datasetX)
  p = ncol(datasetX)
  if(n != length(datasetY))
  {
    stop("Dimensions of datasetX and datasetY do not match!")
  }
  if(method == "SIS")
  {
    
  } else if(method == "KS"){
    
  } else if(method == "CVBF"){
    
  } else if(method == "PT"){
    
  }
  return(ImportantVars)
}