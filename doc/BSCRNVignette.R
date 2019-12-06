## ---- include = FALSE---------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ----setup and Polya tree-----------------------------------------------------
library(BSCRN)
set.seed(500)
datasetsample1 = rnorm(600)
PTsample1 = PolyaTreePriorLikCons(datasetX = datasetsample1, Ginv = qnorm)
names(PTsample1)


## ----drawing samples from the Polya tree--------------------------------------
sampledraws = PolyaTreePredDraws(PTsample1, ndraw = 500)


## ----plot, fig.width = 7, fig.height= 5---------------------------------------

plot(density(sampledraws), main = "Comparison of simple predictive posteriors", xlab = "Red line is true density, black line is approximate Polya Tree predictive posterior", ylim = c(-.01,.41))
lines(density(datasetsample1), col = "blue")
curve(dnorm(x), add=TRUE, col = "red")

## ----Drawing bandwidth samples from CVBF, fig.width = 7, fig.height= 5--------
trainingindices1 = sample(1:600, size = 300)
XT1 = datasetsample1[trainingindices1]
XV1 = datasetsample1[-trainingindices1]
predbwvec1 = PredCVBFIndepMHbw(ndraw = 500, maxIter = 5000, XT1 = XT1, XV1 = XV1)

predbwvec2 = PredCVBFMHbw(ndraw = 500, maxIter = 4000, XT1 = XT1, XV1 = XV1)
 
plot(predbwvec1$predbwsamp)
 
plot(predbwvec2$predbwsamp)

#Independence Metropolis acceptance rate

predbwvec1$acceptancetot / predbwvec1$drawtot

#Metropolis algorithm acceptance rate

predbwvec2$acceptancetot / predbwvec2$drawtot

## ----Drawing from the Predictive posterior of CVBF, fig.width = 7, fig.height= 5----
predpostsamp = PredCVBFDens(predbwvec1$predbwsamp, XT1 = XT1)

predpostsamp2 = PredCVBFDens(predbwvec2$predbwsamp, XT1 = XT1)

plot(seq(from = min(datasetsample1), to = max(datasetsample1), length.out = 200) , predpostsamp(seq(from = min(datasetsample1), to = max(datasetsample1), length.out = 200)), main = "Comparison of simple predictive posteriors", xlab = "Red line is true density, black dots is approximate CVBF predictive posterior", ylab = "density")
lines(density(datasetsample1), col = "blue")
curve(dnorm(x), add=TRUE, col = "red")



plot(seq(from = min(datasetsample1), to = max(datasetsample1), length.out = 200) , predpostsamp2(seq(from = min(datasetsample1), to = max(datasetsample1), length.out = 200)), main = "Comparison of simple predictive posteriors", xlab = "Red line is true density, black dots is approximate CVBF predictive posterior", ylab = "density")
lines(density(datasetsample1), col = "blue")
curve(dnorm(x), add=TRUE, col = "red")



## ----Applying Polya Tree bayesian procedures for testing----------------------
#Method done by constructing two Polya tree objects, and testing the two for differences
datasetsample2 = rnorm(600)
PTsample2 = PolyaTreePriorLikCons(datasetX = datasetsample2, Ginv = qnorm)
logBF = PolyaTreeBFcons(PTsample1, PTsample2)

#Method directly computes the log BF for two data sets. The same values are available

altsamelogBF = PolyaTreetest(datasetX = datasetsample1, datasetY = datasetsample2, Ginv = qnorm)

logBF
altsamelogBF

datasetsample3 = rnorm(600, mean = 2, sd = 1)
PTsample3 = PolyaTreePriorLikCons(datasetX = datasetsample3, Ginv = qnorm)

logBF2 = PolyaTreeBFcons(PTsample1, PTsample3)
altsamelogBF2 = PolyaTreetest(datasetX = datasetsample3, datasetY = datasetsample1, Ginv = qnorm)
logBF2
altsamelogBF2


## ----Applying CVBF bayesian procedures for testing----------------------------

CVBF1 = CVBFtestrsplit(dataset1 = datasetsample1, dataset2 = datasetsample2, trainsize1 = 300, trainsize2 = 300, train1_ids = trainingindices1, train2_ids = trainingindices1)

CVBF2 = CVBFtestrsplit(dataset1 = datasetsample1, dataset2 = datasetsample3, trainsize1 = 300, trainsize2 = 300, train1_ids = trainingindices1, train2_ids = trainingindices1)

CVBF1$logBF

CVBF2$logBF


## ----Applying These tests as a screening method for classification------------
data(gisettetrainlabs)
data(gisettetrainpreds)

#There are also validation labels and validation predictors
#Can be included via
#data(gisettevalidlabs)
#data(gisettevalidpreds)

dim(gisettetrainpreds)
head(gisettetrainlabs)
dim(gisettetrainlabs)
nworkers = detectCores()

ImpVarsSIS1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "SIS", ncores = nworkers / 2)

length(ImpVarsSIS1$varspicked)

ImpVarsKS1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "KS", ncores = nworkers / 2)

length(ImpVarsKS1$varspicked)

ImpVarsPT1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "PT", ncores = nworkers / 2, c = 1, leveltot = 12, Ginv = qnorm, PTscale = TRUE)
 #Only do on first 500
length(ImpVarsPT1$varspicked)
hist(ImpVarsPT1$logBFlist)

#
ImpVarsCVBF1 = ParScreenVars(datasetX = gisettetrainpreds[, 1:500], datasetY = gisettetrainlabs[,1], method = "CVBF", ncores = nworkers / 2, trainsize1 = 2960, trainsize2 = 2960, seed = 200)
length(ImpVarsCVBF1$varspicked)
hist(ImpVarsCVBF1$logBFlist)



