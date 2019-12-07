#' Labels for data that represent 4 or 9 in image recognition. It is unknown which label corresponds to which image. This data set was given as a blind classification problem by NIPS in 2003. The competition involves guessing what some validation labels were. These were the training labels ued for statisticains to train their models before submission.
#' 
#' @format A named matrix with 6000 rows and 1 column, each value is -1 or 1 and corresponds to whether or not an image is a 4 or 9 in the training data set.
#'
#' @source  \url{https://archive.ics.uci.edu/ml/datasets/Gisette}
#'
"gisettetrainlabs"


#' Predictors of the labels in gisettetrainlabs. Some of these correspond to pixel intensity measurements of some images, while others are random noise. It is uncertain which are noisy probes and which are pixel measurements.
#' 
#' @format A named matrix with 6000 rows and 5000 columns. It is zero inflated, as not every pixel is touched when a person writes a number.
#'
#' @source  \url{https://archive.ics.uci.edu/ml/datasets/Gisette}
#'
"gisettetrainpreds"

#' More labels for data that represent 4 or 9 in image recognition. It's unknown which label corresponds to which image. This data set was given as a blind classification problem by NIPS in 2003. The competition involves guessing what the validation labels were. These were released after the competition.
#' 
#' @format A named matrix with 1000 rows and 1 column. 
#'
#' @source  \url{https://archive.ics.uci.edu/ml/datasets/Gisette}
#'
"gisettevalidlabs"

#' Predictors of the labels in gisettevalidlabs. Some of these correspond to pixel intensity measurements of some images, while others are random noise. It is uncertain which are noisy probes and which are pixel measurements.
#' 
#' @format A named matrix with 1000 rows and 5000 columns. 
#'
#' @source  \url{https://archive.ics.uci.edu/ml/datasets/Gisette}
#'
"gisettevalidpreds"