% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CVBFFunctions.R
\name{PredCVBFDens}
\alias{PredCVBFDens}
\title{Compute a Predictive Posterior function given a sequence of bandwidths from the predictive posterior}
\usage{
PredCVBFDens(bwvec, XT1)
}
\arguments{
\item{bwvec}{A vector of bandwidths, can come from either of the methods that draw bandwidths from the posterior.}

\item{XT1}{The training set}
}
\value{
A function that evaluates the predictive posterior at particular values
}
\description{
Compute a Predictive Posterior function given a sequence of bandwidths from the predictive posterior
}
\examples{
set.seed(500)
datasetsample1 = rnorm(600)
trainingindices1 = sample(1:600, size = 300)
XT1 = datasetsample1[trainingindices1]
XV1 = datasetsample1[-trainingindices1]
predbwvec1 = PredCVBFIndepMHbw(ndraw = 500, maxIter = 5000, XT1 = XT1, XV1 = XV1)
predpost = PredCVBFDens(predbwvec1, XT1)
plot(seq(from = -3, to = 3, by = .1), Predpost(predbwvec1, XT1)(seq(from = -3, to = 3, by = .1)))

}
