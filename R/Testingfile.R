datasetX = rnorm(500)

datasetY = rnorm(500)


sampPT1 = PolyaTreePriorLikCons(datasetX = datasetX)

sampPT2 = PolyaTreePriorLikCons(datasetX = datasetY)

PolyaTreetest(datasetX = datasetX, datasetY = datasetY)

#This is verified as correct

PolyaTreeBFcons(sampPT1, sampPT2)

#The two match!

#We're good to go!
 
sampledraws = PolyaTreePredDraws(sampPT1, ndraw = 200)



length(datasetX)
datasetX = rnorm(500)
XT1 = datasetX[1:250]
XV1 = datasetX[251:500]

predbwvec1 = PredCVBFIndepMHbw(ndraw = 600, propsd = 0.01, maxIter = 1000, XT1 = XT1, XV1 = XV1)

predbwvec2 = PredCVBFMHbw(ndraw = 100, propsd = 0.1, maxIter = 1000, XT1 = XT1, XV1 = XV1)

plot(predbwvec1$predbwsamp)

plot(predbwvec2$predbwsamp)

predpostsamp = PredCVBFDens(predbwvec1$predbwsamp, XT1 = XT1)

plot(seq(from = -4, to = 4, length.out = 100) , predpostsamp(seq(from = -4, to = 4, length.out = 100)))


#toc()

#Mixing seems nice enough

# library(coda)
# plot(predbwvec2$predbwsamp)
# autocorr.plot(predbwvec2$predbwsamp, auto.layout = FALSE)
# 
# plot(predbwvec1$predbwsamp)
# autocorr.plot(predbwvec1$predbwsamp, auto.layout = FALSE)
# 
# #Maybe take every 10th sample?
# 
# thinnedbwvec = bwvec[seq(1,length(bwvec),10)]
# #thinnedbwvec = thinnedbwvec[-1]
# 
# predpostavg = 0
# for(j in 1:length(thinnedbwvec))
# {
#   predpostavg = (1/length(thinnedbwvec))*HallKernel(thinnedbwvec[j], datagen = XT1, seq(from = -4, to = 4, length.out = 1000)) + predpostavg
# }
# 
# plot(seq(from = -4, to = 4, length.out = 1000), predpostavg, col = "purple")
# lines(seq(from = -4, to = 4, length.out = 1000), dnorm(seq(from = -4, to = 4, length.out = 1000)))
# 
# plot(seq(from = -10, to = 10, length.out = 1000), dnorm(seq(from = -10, to = 10, length.out = 1000)))
# lines(seq(from = -10, to = 10, length.out = 1000), predpostavg, col = "purple")
# 
# 
# 
# 
# 
# XT1 = datasetX[1:250]
# XV1 = datasetX[251:500]
# likvec = function(h) {sum(log(HallKernel(h,datagen2 = XT1, x = XV1)))}
# bwlik = optimize(f = function(h){  likvec(h)}, lower = 0, upper = 10, maximum = TRUE)
# 
# bwvec = c()
# bwvec[1] = bwlik$maximum
# beta = bwlik$maximum
# postcurr = sum(log(HallKernel(bwvec[1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[1]) - (beta^2 / bwvec[1]^2)
# #Lets try a normal proposal first 
# #N(0,.05)
# acceptances = 0
# #tic()
# iter = 10000
# for(j in 1:iter)
# {
#   bwprop = rnorm(1, mean = bwvec[j], sd = .05)
#   if(bwprop <= 0)
#   {
#     bwvec[j+1] = bwvec[j]
#   }
#   else
#   {
#     postprop = sum(log(HallKernel(bwprop, datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwprop) - (beta^2 / bwprop^2)
#     logpostdif = postprop - postcurr
#     if(exp(logpostdif) > runif(1))
#     {
#       bwvec[j+1] = bwprop
#       postcurr = sum(log(HallKernel(bwvec[j+1], datagen2 = XT1, x = XV1))) + log(2*beta) - .5*log(pi) - 2*log(bwvec[j+1]) - (beta^2 / bwvec[j+1]^2)
#       acceptances = acceptances + 1
#     }
#     else
#     {
#       bwvec[j+1] = bwvec[j]
#     }
#   }
# }
# #toc()
# plot(bwvec)
# 
# predpostsamp2 = PredCVBFDens(bwvec = bwvec, XT1 = XT1)
# 
# plot(seq(from = -4, to = 4, length.out = 100) , predpostsamp2(seq(from = -4, to = 4, length.out = 100)))


#head(unname(gisettetrainpreds))

data(gisettetrainlabs)
data(gisettetrainpreds)

unname(gisettetrainpreds[1:10,1:10])



# ImpVarsSIS1 = ParScreenVars(datasetX = gisettetrainpreds, datasetY = gisettetrainlabs[,1], method = "SIS", ncores = 1)
# 
# ImpVarsKS1 = ParScreenVars(datasetX = gisettetrainpreds, datasetY = gisettetrainlabs[,1], method = "KS", ncores = 1)
# 
# ImpVarsPT1 = ParScreenVars(datasetX = gisettetrainpreds[,1:500], datasetY = gisettetrainlabs[,1], method = "PT", ncores = 10, c = 1, leveltot = 12)
# #Only do on first 500
# length(ImpVarsPT1)
# 
# ImpVarsCVBF1 = ParScreenVars(datasetX = gisettetrainpreds[,1:100], datasetY = gisettetrainlabs[,1], method = "CVBF", ncores = 10, trainsize1 = 2960, trainsize2 = 2960, seed = 200)
# 
# ImpVarsCVBF2 = ParScreenVars(datasetX = gisettetrainpreds[,1:40], datasetY = gisettetrainlabs[,1], method = "CVBF", ncores = 10, trainsize1 = 2960, trainsize2 = 2960, seed = 100)
# 
#Suprisingly (or unsuprisingly) changing seed alters which vars are picked
