library(parallel)
library(ggplot2)
setwd("C:/Users/navee/OneDrive/Documents/NonparamCVKernelBasedBF/NonparamCVKernelBasedBF/StandaloneSimulationfolder")
source("MarginalLikIntfunctions.R")
source("Laplacefunction.R")
set.seed(1000)

dlength = c(200,500,1000)
dlength = rep(dlength, each = 500)

i = 0
NormalLapRes = data.frame(LaplaceApprox = rep(0, times = length(dlength)),
           QuadApprox = rep(0, times = length(dlength))
           )

CauchyLapRes = data.frame(LaplaceApprox = rep(0, times = length(dlength)),
                          QuadApprox = rep(0, times = length(dlength))
)

while(i < length(dlength))
{
  dataset1 <- rnorm(dlength[i])
  
  XT1 <- dataset1[1:(length(dataset1)*.25)]
  XV1 <- dataset1[-(1:(length(dataset1)*.25))]
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
  bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML1 = laplace.kernH2(y = XT1, x = XV1, hhat = bwlik$maximum)
  
  QuadResult = logmarg.kernH(y = XT1, x = XV1)
  NormalLapRes$LaplaceApprox[i] = ExpectedKernML1[1]
  NormalLapRes$QuadApprox[i] = QuadResult[[2]][1]
  
  dataset1 <- rcauchy(dlength[i])
  
  XT1 <- dataset1[1:(length(dataset1)*.25)]
  XV1 <- dataset1[-(1:(length(dataset1)*.25))]
  
  likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
  bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
  ExpectedKernML1 = laplace.kernH2(y = XT1, x = XV1, hhat = bwlik$maximum)
  
  QuadResult = logmarg.kernH(y = XT1, x = XV1)
  CauchyLapRes$LaplaceApprox[i] = ExpectedKernML1[1]
  CauchyLapRes$QuadApprox[i] = QuadResult[[2]][1]
  print(i)
  i = i + 1
}

CauchyLapRes$n = dlength
NormalLapRes$n = dlength
save.image("C:/Users/navee/OneDrive/Documents/NonparamCVKernelBasedBF/NonparamCVKernelBasedBF/StandaloneSimulationfolder/QuadratureVLapRes.RData")

load.image("C:/Users/navee/OneDrive/Documents/NonparamCVKernelBasedBF/NonparamCVKernelBasedBF/StandaloneSimulationfolder/QuadratureVLapRes.RData")

plot(x = CauchyLapRes$LaplaceApprox[CauchyLapRes$n == 500], y = CauchyLapRes$QuadApprox[CauchyLapRes$n == 500],
     xlab = "Quadrature", ylab = "Laplace Approximation", main = "Cauchy Data")

plot(x = NormalLapRes$LaplaceApprox[NormalLapRes$n == 500], y = NormalLapRes$QuadApprox[NormalLapRes$n == 500],
     xlab = "Quadrature", ylab = "Laplace Approximation", main = "Normal Data")

CauchyLapRes$QuadLapDiff = (CauchyLapRes$LaplaceApprox - CauchyLapRes$QuadApprox) / CauchyLapRes$QuadApprox

NormalLapRes$QuadLapDiff = (NormalLapRes$LaplaceApprox - NormalLapRes$QuadApprox) / NormalLapRes$QuadApprox

aggregate(QuadLapDiff ~ n, data = NormalLapRes, summary)

aggregate(QuadLapDiff ~ n, data = CauchyLapRes, summary)
