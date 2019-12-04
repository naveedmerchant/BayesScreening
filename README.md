# BayesScreening Package
Naveed Merchant
  - [Installation](#installation)
  - [Usage](#usage)
      - [1. Parallel Screening of data for binary classification]
         (#Parallel Screening of data for binary classification)
      - [2. Testing two data sets to see if they share the same distribution
](#Testing two data sets to see if they share the same distribution)
      - [3. Generating samples from the predictive posterior of the two previous tests](#generate-data)

## Introduction

BayesScreening is an R package for screening important variables for data sets with large amounts of predictors. This screening is intended for binary classification. It does this by applying Bayesian tests that check if two or more samples have the same distribution. There are also other screening methods to compare these methods to in this package. This package also supports running this screening in parallel.

The tests themselves may also be useful alone. The functions to run the tests themselves are in the package. As the tests are nonparametric bayesian tests, a person may also be interested in sampling from the predictive posterior of these bayesian tests. There is functionality available to do this as well. The GISETTE dataset is also supplied through this package. If desired, downloading the package will give access to the GISETTE dataset.

Intended users of this package are people that either want to use the screening methods, want to use the tests to check if the distributions are the same, or are people that want to see bayesian predictive posteriors that estimate the distribution of the data under few assumptions.


## Installation

``` r
devtools::install_github("naveedmerchant/BayesScreening")
```

## Usage

``` r
library(BSCRN)
```

## Examples

### 1\. Parallel Screening of data for binary classification

For data sets with a large number of predictors, its rare that all predictors are truly useful. This function takes
a data set that contains a column corresponding to classes, and a matrix containing predictors. It returns a 
list of indices that correspond to the predictors that are believed to be useful for screening. It does so in parallel. 

``` r
nworkers = detectCores()
#We'll use nworkers - 1 to be careful not to overwork the PC
data(gisettetrainpreds)
data(gisettetrainlabs)
#SIS is a classic method
SISvars = ParScreenVars(datasetX = gisettetrainpreds, datasetY = gisettetrainlabs[,1], method = "SIS", ncores = nworkers - 1)
#Return number of vars that SIS deemed important
length(SISvars)

#PT or Polya tree is a bayesian method
PTvars = ParScreenVars(datasetX = gisettetrainpreds, datasetY = gisettetrainlabs[,1], method = "PT", ncores = nworkers - 1)
#May take a little more time to run
#Return number of vars that PT deemed important
length(PTvars)
```

### 2\. Testing two data sets to see if they share the same distribution

Testing to see if two data sets share the same distribution is a classic problem in statistics, but it can be tricky to do if there isn't a distribution specified on the two data sets. This package comes with two bayesian functions that perform this test, but do it without specifying a distribution. The value returned by these tests are log Bayes factors, not p-values. If they are positive and not close to 0, that should imply that the two distributions are different. If they are negative and not close to 0, that should imply the two distributions are the same. If they are close to 0, the test is inconclusive. 

``` r
set.seed(100)
# generate noise with same distribution
dataset1 = rnorm(200)
dataset2 = rnorm(200) 
logPTBF1 = PolyaTreetest(datasetX = dataset1, datasetY = dataset2)
logPTBF1
logCVBF1 = CVBFtestrsplit(dataset1, dataset2, trainsize1 = 100, trainsize2 = 100)
logCVBF1

#generate noise with different distribution 
dataset3 = rnorm(200, mean = 0, sd = 4)
logPTBF1 = PolyaTreetest(datasetX = dataset1, datasetY = dataset3)
logPTBF1
logCVBF1 = CVBFtestrsplit(dataset1, dataset3, trainsize1 = 100, trainsize2 = 100)
logCVBF1

```

### 3\. Generating samples from the predictive posterior of the two previous tests

For a bayesian non-parametric test to work, it should adopt some idea as to what the underlying distribution of the data set really is.
Examining the predictive posterior of these tests gives an idea as to what the bayesian non-parametric procedures finds the underlying distribution of the data to be.

``` r
#Polya tree predictive posterior
dataset1 = rnorm(200)
sampPT1 = PolyaTreePriorLikCons(datasetX = dataset1)
sampledraws = PolyaTreePredDraws(sampPT1, ndraw = 200)
plot(density(sampledraws))

#CVBF predictive posterior
XT1 = datasetX[1:100]
XV1 = datasetX[101:200]
predbwvec1 = PredCVBFIndepMHbw(ndraw = 200, propsd = 0.01, maxIter = 1000, XT1 = XT1, XV1 = XV1)
predpostsamp = PredCVBFDens(predbwvec1$predbwsamp, XT1 = XT1)
plot(seq(from = min(dataset1), to = max(dataset1), length.out = 100) , predpostsamp(seq(from = min(dataset1), to = max(dataset1), length.out = 100)))

plot(seq(from = min(dataset1), to = max(dataset1), length.out = 100) , dnorm(seq(from = min(dataset1), to = max(dataset1), length.out = 100)))

```
