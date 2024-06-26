HiggsCSV = file.choose()

#setwd("C:/Users/Naveed/Dropbox/NonparamCVKernelBasedBF")
fulltrans = read.csv(file = HiggsCSV, colClasses = c(NA,rep("NULL",times = 21),rep(NA,times = 7)), header = FALSE)

#fulltrans = read.csv(file = "C:/Users/Naveed/Documents/HIGGS.csv.gz", colClasses = c(rep("NULL",times = 22),rep(NA,times = 7)), header = FALSE)
dim(fulltrans)

noisedataind = fulltrans[,1]== 0

X = fulltrans[noisedataind,]
#This is "background noise"
Y = fulltrans[!noisedataind,]


source("MarginalLikIntfunctions.R")
source("Laplacefunction.R")

usedsamptrans = fulltrans2[1:20000,]
noisedataind3 = usedsamptrans[,1]== 0
X2 = usedsamptrans[noisedataind3,]
#This is "background noise"
Y2 = usedsamptrans[!noisedataind3,]


plot(density(Y[,8], n = 511, from = 0, to = 4), xlab = "x", main = "")
lines(density(X[,8], n = 511, from = 0, to = 4), col = "blue")



plot(density(Y[,2], n = 511, from = 0, to = 3), xlab = "x", main = "")
lines(density(X[,2], n = 511, from = 0, to = 3), col = "blue")

fulltrans2 = fulltrans[1:1000000,]

noisedataind2 = fulltrans2[,1]== 0
X2 = fulltrans2[noisedataind,]
#This is "background noise"
Y2 = fulltrans2[!noisedataind,]

source("MarginalLikIntfunctions.R")
source("Laplacefunction.R")


####Code to create plots for 23rd col of Higgs Boson

usedsamptrans = fulltrans2[1:30000,]
noisedataind3 = usedsamptrans[,1]== 0
X2 = usedsamptrans[noisedataind3,]


X2 = X2[1:10000,]
#This is "background noise"
Y2 = usedsamptrans[!noisedataind3,]

Y2 = Y2[1:10000,]
iter = 20
trainsizelist = seq(from = 1000, by = 1000, to = 5000)
#trainsizelist = ceiling(exp(seq(from = log(1000), to = log(5000), length.out = 20)))
#validsize = 5000:10000
#logBFmat3 = matrix(nrow = length(trainsizelist), ncol = iter)
#The below code does essentially give the correct BFs that were provided below. Make prior BFs with same code
set.seed(1000)

#logBFmat6 = matrix(nrow = length(trainsizelist), ncol = iter)
logBFmat3 = matrix(nrow = length(trainsizelist), ncol = iter)

#logBFmat8 = matrix(nrow = length(trainsizelist), ncol = 3)
for(j in 1:trainsizelist)
{
  for(k in 1:iter)
  {
    trainindX = sample(1:nrow(X2) ,size = trainsizelist[j])
    trainindY = sample(1:nrow(Y2) ,size = trainsizelist[j])
    XT1 <- X2[trainindX,2]
    XV1 <- X2[-(trainindX),2]
    YT1 <- Y2[(trainindY),2]
    YV1 <- Y2[-(trainindY),2]
    
    #This is almost arbitrarily chosen, we proceed like this for now...
    tic()
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
    bwlik2 = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    ExpectedKernML1 = laplace.kernH2c(y = XT1, x = XV1, hhat = bwlik2$maximum, c = bwlik2$objective)
    
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = YT1, x = YV1)))}
    
    bwlikcy = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    
    ExpectedKernML2 = laplace.kernH2c(y = YT1, x = YV1, hhat = bwlikcy$maximum, c= bwlikcy$objective)
    
    likveccombc = function(h) {sum(log(HallKernel(h,datagen2 = c(XT1,YT1), x = c(XV1,YV1))))}
    bwlikcombc = optimize(f = function(h){  likveccombc(h)}, lower = 0, upper = 10, maximum = TRUE)
    
    ExpectedKernMLcomb = laplace.kernH2c(y = c(XT1,YT1), x = c(XV1,YV1), hhat = bwlikcombc$maximum, c= bwlikcombc$objective)
    
    logBFmat3[j,k] = ExpectedKernML1[1] + ExpectedKernML2[1] - ExpectedKernMLcomb[1]
    
    
    toc()
    print(j)
  }
}
logBFmeans3 = rowMeans(logBFmat3)
logBFmedains3 = rowMedians(logBFmat3)
logBFranges3 = rowRanges(logBFmat3)

dfmat3 <- data.frame(trainsize=trainsizelist, min=logBFranges3[,1], max=logBFranges3[,2], mean = logBFmeans3)
library(ggplot2)
ggplot(dfmat3, aes(x=trainsize))+
  geom_linerange(aes(ymin=min,ymax=max),linetype=2,color="blue")+
  geom_point(aes(y=min),size=3,color="red")+
  geom_point(aes(y=max),size=3,color="red")+
  geom_point(aes(y=mean),size = 4, color = "green")+
  geom_line(aes(x = trainsize, y = mean), size = 3, color = "purple", linetype = 2)
theme_bw()

-PolyaTreetest(X2[,2],Y2[,2])

###Code for figure of 29th col

logBFmat3 = matrix(nrow = length(trainsizelist), ncol = iter)

#logBFmat8 = matrix(nrow = length(trainsizelist), ncol = 3)
for(j in 1:trainsizelist)
{
  for(k in 1:iter)
  {
    trainindX = sample(1:nrow(X2) ,size = trainsizelist[j])
    trainindY = sample(1:nrow(Y2) ,size = trainsizelist[j])
    XT1 <- X2[trainindX,8]
    XV1 <- X2[-(trainindX),8]
    YT1 <- Y2[(trainindY),8]
    YV1 <- Y2[-(trainindY),8]
    
    #This is almost arbitrarily chosen, we proceed like this for now...
    tic()
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
    bwlik2 = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    ExpectedKernML1 = laplace.kernH2c(y = XT1, x = XV1, hhat = bwlik2$maximum, c = bwlik2$objective)
    
    likvec = function(h) {sum(log(HallKernel(h,datagen2 = YT1, x = YV1)))}
    
    bwlikcy = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
    
    ExpectedKernML2 = laplace.kernH2c(y = YT1, x = YV1, hhat = bwlikcy$maximum, c= bwlikcy$objective)
    
    likveccombc = function(h) {sum(log(HallKernel(h,datagen2 = c(XT1,YT1), x = c(XV1,YV1))))}
    bwlikcombc = optimize(f = function(h){  likveccombc(h)}, lower = 0, upper = 10, maximum = TRUE)
    
    ExpectedKernMLcomb = laplace.kernH2c(y = c(XT1,YT1), x = c(XV1,YV1), hhat = bwlikcombc$maximum, c= bwlikcombc$objective)
    
    logBFmat3[j,k] = ExpectedKernML1[1] + ExpectedKernML2[1] - ExpectedKernMLcomb[1]
    
    
    toc()
    print(j)
  }
}
logBFmeans3 = rowMeans(logBFmat3)
logBFmedains3 = rowMedians(logBFmat3)
logBFranges3 = rowRanges(logBFmat3)

dfmat3 <- data.frame(trainsize=trainsizelist, min=logBFranges3[,1], max=logBFranges3[,2], mean = logBFmeans3)
library(ggplot2)
ggplot(dfmat3, aes(x=trainsize))+
  geom_linerange(aes(ymin=min,ymax=max),linetype=2,color="blue")+
  geom_point(aes(y=min),size=3,color="red")+
  geom_point(aes(y=max),size=3,color="red")+
  geom_point(aes(y=mean),size = 4, color = "green")+
  geom_line(aes(x = trainsize, y = mean), size = 3, color = "purple", linetype = 2)
theme_bw()

-PolyaTreetest(X2[,8],Y2[,8])

################## DO ABOVE PLOTS IN PARALLEL IF POSSIBLE TO SAVe TIME (should help a lot)

library(tictoc)

trainsizelist = seq(from = 1000, by = 1000, to = 5000)

set.seed(1000)
bwvec2 = c()

trainindX = sample(1:nrow(X2) ,size = trainsizelist[5])
trainindY = sample(1:nrow(Y2) ,size = trainsizelist[5])
XT1 <- X2[trainindX,2]
XV1 <- X2[-(trainindX),2]

likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)

beta = bwlik$maximum
bwvec2[1] = bwlik$maximum
postcurr = sum(log(HallKernel(bwvec2[1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec2[1]) - (beta^2 / bwvec2[1]^2)
#Lets try a normal proposal first 
#N(0,.0005)
acceptances = 0
tic()
iter = 3000
for(j in 1:iter)
{
  bwprop = rnorm(1, mean = bwvec2[j], sd = .0025)
  if(bwprop <= 0)
  {
    bwvec2[j+1] = bwvec2[j]
  }
  else
  {
    postprop = sum(log(HallKernel(bwprop, datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwprop) - (beta^2 / bwprop^2)
    logpostdif = postprop - postcurr
    if(exp(logpostdif) > runif(1))
    {
      bwvec2[j+1] = bwprop
      postcurr = sum(log(HallKernel(bwvec2[j+1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec2[j+1]) - (beta^2 / bwvec2[j+1]^2)
      acceptances = acceptances + 1
    }
    else
    {
      bwvec2[j+1] = bwvec2[j]
    }
  }
  print(j)
}
toc()
plot(bwvec2)
#Its gonna take a lot more time than PT, this can't easily be sped up, but can be parallelized with Independence Metrop sampler
#Maybe approximate KDE with less points and proceed (this sounds like an approximate techinique though idk)
#Independence sampler also PROBABLY has higher acceptance rates since predictive posterior of the BW is close to a 
#Normal distribution with the mean being around the place where the bw maximizes the likelihood.
#Mixing seems nice enough

tic()
sum(log(HallKernel(bwvec2[j+1], datagen2 = XT1, x = XV1)))
toc()

#Roughly 2.75 seconds to evaluate


library(coda)

autocorr.plot(bwvec2, auto.layout = FALSE)
#Maybe take every 10th sample?

thinnedbwvec2 = bwvec2[seq(1,length(bwvec2),10)]

predpostavg2 = 0
for(j in 1:length(thinnedbwvec2))
{
  predpostavg2 = (1/length(thinnedbwvec2))*HallKernel(thinnedbwvec2[j], datagen = XT1, seq(from = 0, to = 4, length.out = 10000)) + predpostavg2
}

#################################Code for fig 3



set.seed(500)
datalist = X2[,2]

alphalist = list()
leveltot = 9
c = 1
for(m in 1:leveltot)
{
  alphalist[[m]] = rep(c*m^2,2^m)
}

#An alternative way to write this?
epsilonlist = c()
epsilonlist2 = list()
splitlist = list()
for(j in 1:leveltot)
{
  splitlist[[j]] = rep(0,2^j)
}

for(j in 1:length(datalist))
{
  thresh = .5
  for(m in 1:leveltot)
  {
    epsilonlist[m] = datalist[j] > qnorm(thresh)
    multivec = 2^seq(m-1, 0, -1)
    ind = sum(epsilonlist[1:m]*multivec) + 1
    splitlist[[m]][ind] = splitlist[[m]][ind] + 1
    thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
  }
  epsilonlist2[[j]] <- epsilonlist
  epsilonlist = c()
}
splitlist

# epsilonlist3 = list()
# splitlist2 = list()
# for(j in 1:leveltot)
# {
#   splitlist2[[j]] = rep(0,2^j)
# }
# for(j in 1:length(datalist2))
# {
#   thresh = .5
#   for(m in 1:leveltot)
#   {
#     epsilonlist[m] = datalist2[j] > qnorm(thresh)
#     multivec = 2^seq(m-1, 0, -1)
#     ind = sum(epsilonlist[1:m]*multivec) + 1
#     splitlist2[[m]][ind] = splitlist2[[m]][ind] + 1
#     thresh = thresh + (2*epsilonlist[m] - 1) / 2^{m+1}
#   }
#   epsilonlist3[[j]] <- epsilonlist
#   epsilonlist = c()
# }
# splitlist2
# 

largeqnormlist = seq(from = 0, to = 1, by = .0000001)
largeqnormlist = qnorm(largeqnormlist)


#Want a sample from distribution that has polya tree prior

postbetalist = list()
for(j in 1:leveltot)
{
  postbetalist[[j]] = splitlist[[j]] + alphalist[[j]]
}

#posteriorPTdrawlist = list()

posteriorPTdraw = c()
leveltot = 9
ndraw = 1000
k = 1
indlist = c()
for(i in 1:ndraw)
{
  for(m in 1:leveltot)
  {
    postdraw = rbinom(1,size = 1, prob = (postbetalist[[m]][k]) / (postbetalist[[m]][k+1] + postbetalist[[m]][k]))
    if(postdraw == 1)
    {
      epsilonlist[m] = 0
      #k = k 
    }
    else
    {
      epsilonlist[m] = 1
      #k = k + 2
    }
    #epsilonlist[m] = datalist[j] > qnorm(thresh)
    multivec = 2^seq(m, 1, -1)
    ind = sum(epsilonlist[1:m]*multivec) + 1
    
    #print(postbetalist[[m]][k])
    #print(postbetalist[[m]][k+1])
    #print(postdraw)
    #print(ind)
    k = ind
    
    #splitlist[[m]][ind] = splitlist[[m]][ind] + 1
  }
  #epsilonlist2[[j]] <- epsilonlist
  #need to draw from qnorm((ind - 1) / 2^9)  qnorm(ind / 2^9)
  #VFY below line
  posteriorPTdraw[i] = sample(largeqnormlist[(largeqnormlist < qnorm((ind+1) / 2^(leveltot+1))) & (largeqnormlist > qnorm((ind) / 2^(leveltot+1)))], size = 1)
  epsilonlist = c()
  indlist[i] = ind
  k=1
  print(i)
}

plot(density(fulltrans[noisedataind,2], from = 0, to = 4, n = 10000), lwd = 2)
lines(seq(from = 0, to = 4, length.out = 10000), predpostavg2, col = "purple", main = "Predictive Posteriors vs True Density Estimate", xlab = "23rd column values of Noise", ylab = "Density", lwd = 2)

lines(density(posteriorPTdraw, from = 0, to = 4, n = 10000), col = "red", lwd = 2)




