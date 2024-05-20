source("C:/Users/Naveed/Dropbox/NonparamCVKernelBasedBF/CVKernelBF/Screening_Method/ArceneAnalysis/Important R functions/MarginalLikIntfunctions.R")
source("C:/Users/Naveed/Dropbox/NonparamCVKernelBasedBF/CVKernelBF/Screening_Method/ArceneAnalysis/Important R functions/Laplacefunction.R")

###The goal of these simulations is to try to see what the ideal size r, is.

###We generate data from the compact support case

###We combine the data sets and try to see what training size produces Bayes factors that are good in a 
###"Type I error sense", (pick log BF such that in the null case (which we can simulate by combining data sets))
### 

p = rbeta(500,.5,.5)
dlength = 400

dataset1 <- runif(dlength)
dataset2 <- rbeta(dlength, shape1 = 4, shape2 = 4)


dataset1 <- runif(dlength)
dataset2 <- rbeta(dlength, shape1 = 4, shape2 = 4)

logBFlist = c()
for(j in 1:20)
{
  logBFlist[j] = BSCRN::CVBFtestrsplit(dataset1 = dataset1, dataset2 = dataset2, trainsize1 = 140, trainsize2 = 140)$logBF
}
  library(BSCRN)

set.seed(1000)

jointset = c(dataset1, dataset2)
numreps = 100
propr = seq(from = 80, by = 20, to = 200)
m = length(dataset1)
n = length(dataset2)
logBFmatrs = matrix(nrow = numreps, ncol = length(propr))
logBFdf = data.frame("r" = NA, "logBF" = NA)
z = 1
for(j in 1:numreps)
{
  permdataset1IDs = sample(1:length(jointset), size = m)
  permdataset1 = jointset[permdataset1IDs]
  permdataset2 = jointset[-permdataset1IDs]
  for(k in 1:length(propr))
  {
    s = ceiling((m / n) * propr[k])
    CVBFoutput = CVBFtestrsplit(dataset1 = permdataset1, dataset2 = permdataset2, trainsize1 = propr[k], trainsize2 = s)$logBF
    logBFmatrs[j, k] = CVBFoutput
    logBFdf[z, ] = c(propr[k], logBFmatrs[j, k])
    z = z + 1
  }
  print(j)
}

library(ggplot2)

plot1 = ggplot(data = logBFdf, mapping = aes(x = r, y = logBF)) + geom_point()

plot2 = ggplot(data = logBFdf, mapping = aes(x = logBF, y = r)) + geom_point()

#We can try this setup on one of the cols in the Higgs Boson data set as well.

HiggsCSV = "C:/Users/Naveed/Documents/HIGGS.csv.gz"

fulltrans = read.csv(file = HiggsCSV, colClasses = c(NA,rep("NULL",times = 21),rep(NA,times = 7)), header = FALSE, nrows = 1000000)
dim(fulltrans)

fulltrans2 = fulltrans[1:1000000,]

usedsamptrans = fulltrans2[1:20000,]

noisedataind3 = usedsamptrans[,1]== 0
X2 = usedsamptrans[noisedataind3,]
#This is "background noise"
Y2 = usedsamptrans[!noisedataind3,]

dataset1 <- X2[,8]
dataset2 <- Y2[,8]

library(BSCRN)

set.seed(1000)

jointset = c(dataset1, dataset2)
numreps = 10
propr = seq(from = 1000, by = 1000, to = 5000)
m = length(dataset1)
n = length(dataset2)
logBFmatrs = matrix(nrow = numreps, ncol = length(propr))
logBFdf = data.frame("r" = NA, "logBF" = NA)
z = 1
for(j in 1:numreps)
{
  permdataset1IDs = sample(1:length(jointset), size = m)
  permdataset1 = jointset[permdataset1IDs]
  permdataset2 = jointset[-permdataset1IDs]
  for(k in 1:length(propr))
  {
    s = ceiling((m / n) * propr[k])
    CVBFoutput = CVBFtestrsplit(dataset1 = permdataset1, dataset2 = permdataset2, trainsize1 = propr[k], trainsize2 = s)$logBF
    logBFmatrs[j, k] = CVBFoutput
    logBFdf[z, ] = c(propr[k], logBFmatrs[j, k])
    z = z + 1
    print(k)
  }
  print(j)
}

library(ggplot2)

logBFmeans = aggregate(data = logBFdf, logBF ~ r, FUN = mean)

logBFmins = aggregate(data = logBFdf, logBF~ r, FUN = min)

logBFmaxs = aggregate(data = logBFdf, logBF~ r, FUN = max)

#logBFsummary = cbind(logBFmeans)



logBFsummary = data.frame(r = logBFmeans$r, min = logBFmins$logBF, max = logBFmaxs$logBF, mean = logBFmeans$logBF)

plot3 = ggplot(data = logBFsummary, mapping = aes(x = r, y = mean)) + geom_point(data = logBFdf, aes(x = r, y = logBF)) +
  geom_line(data = logBFsummary, mapping = aes(x = r, y = mean), size = 3, color = "purple") +
  geom_linerange(data = logBFsummary, aes(x = r, ymin = min, ymax = max), linetype=2, color = "blue")

#ggsave("BayesSimPlots/Col29Nullrplots.pdf",  plot = plot3)

plot4 = ggplot(data = logBFdf, mapping = aes(x = logBF, y = r)) + geom_point()









HiggsCSV = "C:/Users/Naveed/Documents/HIGGS.csv.gz"

fulltrans = read.csv(file = HiggsCSV, colClasses = c(NA,rep("NULL",times = 21),rep(NA,times = 7)), header = FALSE, nrows = 1000000)
dim(fulltrans)

fulltrans2 = fulltrans[1:1000000,]

usedsamptrans = fulltrans2[1:20000,]

noisedataind3 = usedsamptrans[,1]== 0
X2 = usedsamptrans[noisedataind3,]
#This is "background noise"
Y2 = usedsamptrans[!noisedataind3,]

dataset1 <- X2[,2]
dataset2 <- Y2[,2]

library(BSCRN)

set.seed(1000)

jointset = c(dataset1, dataset2)
numreps = 10
propr = seq(from = 1000, by = 1000, to = 5000)
m = length(dataset1)
n = length(dataset2)
logBFmatrs = matrix(nrow = numreps, ncol = length(propr))
logBFdf = data.frame("r" = NA, "logBF" = NA)
z = 1
for(j in 1:numreps)
{
  permdataset1IDs = sample(1:length(jointset), size = m)
  permdataset1 = jointset[permdataset1IDs]
  permdataset2 = jointset[-permdataset1IDs]
  for(k in 1:length(propr))
  {
    s = ceiling((m / n) * propr[k])
    CVBFoutput = CVBFtestrsplit(dataset1 = permdataset1, dataset2 = permdataset2, trainsize1 = propr[k], trainsize2 = s)$logBF
    logBFmatrs[j, k] = CVBFoutput
    logBFdf[z, ] = c(propr[k], logBFmatrs[j, k])
    z = z + 1
    print(k)
  }
  print(j)
}

library(ggplot2)

logBFmeans = aggregate(data = logBFdf, logBF ~ r, FUN = mean)

logBFmins = aggregate(data = logBFdf, logBF~ r, FUN = min)

logBFmaxs = aggregate(data = logBFdf, logBF~ r, FUN = max)

#logBFsummary = cbind(logBFmeans)



logBFsummary = data.frame(r = logBFmeans$r, min = logBFmins$logBF, max = logBFmaxs$logBF, mean = logBFmeans$logBF)

plot32 = ggplot(data = logBFsummary, mapping = aes(x = r, y = mean)) + geom_point(data = logBFdf, aes(x = r, y = logBF)) +
  geom_line(data = logBFsummary, mapping = aes(x = r, y = mean), size = 3, color = "purple") +
  geom_linerange(data = logBFsummary, aes(x = r, ymin = min, ymax = max), linetype=2, color = "blue")


ggplot(data = logBFsummary, mapping = aes(x = r, y = mean)) +
  geom_linerange(data = logBFsummary, aes(x = r, ymin = min, ymax = max), linetype=2, color = "blue")
#ggsave("BayesSimPlots/Col23Nullrplots.pdf",  plot = plot3)

plot42 = ggplot(data = logBFdf, mapping = aes(x = logBF, y = r)) + geom_point()

