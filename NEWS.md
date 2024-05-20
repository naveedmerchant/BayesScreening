# CVBF files added

* The intention of this update is to provide code used for the CVBF paper, and better explain what happens in the simulation section of CVBF. 
* This package currently provides functions to more easily implement CVBF and Polya Tree. We defer to the package's documentation and the README of the repo for guidance on how to implement those methods.

## Simulations in CVBF paper

* All simulations can be seen in the CVBFStandaloneSimulationfolder in this github. It contains many Rmarkdown scripts that generate plots for the paper.
* The intention of the simulations is to see how the the CVBF performs when checking if two distributions are the same, or different.

### Simulations for how our test does when testing if two distributions differ and they really differ

* We'll start by exploring the cases when distributions are different, there are a lot of ways this can occur and we tried to explore hard variations of these. We'll call one of these distributions $f$ and the other distribution $g$
  1. Distributions can be different because of a scale difference in tails. One distribution could be $N(0,1)$ ($f$) and the other one $N(0,4)$ ($g$)
  2. They could be different due to a difference in location. One distribution could be a standard Cauchy ($f$) and the other one a Cauchy where the median location is 1 instead ($g$.
  3. They could be different in tail behavior instead. One distribution is Normal ($f$) and the other is Cauchy ($g$).
    * In our paper (and simulations), we deliberately picked the parameters on $f$ and $g$ so that they have the same IQR.
  4. We supposed compact support in our proof, so we looked at a compact support case. $f$ could be a standard uniform and $g$ could be a Beta distribution. We picked the parameters that $\alpha = 4, \beta = 4$ so that the mean and median of the distributions would be the same, but they still have big structural differences.
* We want to explore a variety of distributions with varying degrees on how far apart the two distributions are. To do this, we compared $f$ to $(1-p)f + pg$. We drew $p$ from a Beta distribution for our simulations. When $p$ is closer to 1, the two distributions are more different, where as when $p$ is close to 0, the distributions are more the same.
* Simulations showcasing 1. have the title  'BayesSimShortTailvsShortTailSims'
* Simulations showcasing 2. have the title  'BayesSimLongTailvsLongTailSims' 
* Simulations showcasing 3. have the title  'BayesSimShortTailvsLongTailSims'
* Simulations showcasing 4. have the title  'CompactSupportBayesSims'
* We mentioned in our paper that we require specification of training and testing set. It is possible to run our test across many different splits or choices of training and testing. We ended up testing the performance of our test where we averaged the log BFs across 30 splits, and where we only ran our test across one split.
  * The suffix AVG is given if the log BFs computed averaged across the 30 splits. If there is no suffix, all plots computed, instead are only off of 1 split.
  * The pdfs of plots currently stored in the GitHub are all based off of 30 splits.
* We compared our test with Holmes' Polya Tree (with different choices of base distribution (Normal and Cauchy)) and the Kolmogorov-Smirnov test. Plots of each of the tests against $p$ alone are available too. We defer to the Rmarkdown files themselves for documentation on those.
* To see all the test results plotted against each other, see the pdfs 'jointlogBF' plots in the BayesSimPlots folder. It has two suffixes. One contains the type of simulation indicated ('ShortvShort', 'Compact', 'ShortvLong', 'LongvLong'), and the other contains what info is on the plot.
* A legend is included in plots with the 'legend' suffix, plots with the 'full' suffix contain all the test statistics plotted, and all plots contain lines that smooth out the trajectory of the test statistics over p.

  
### Simulations for how our test does when testing if two distributions differ and they really differ and they're really the same

* We only tested how our test did against Polya Tree, where data is the same, by looking at two standard Normal distributions. We varied the number of samples as (200, 400, 800). We then plotted all the test stats against the number of samples. We did this for 1500 data sets.
* The code that contains this is 'simulationsnullsamplesizes.R'. The plots made from it can be seen in BayesSimPlots as 'logCVBFnullvsn' and 'logPTBFnullvsn'. 
* We ran this pretty early and it took a while. When running the same file, we encourage setting the number of cores used to a much smaller number and letting it run over the course of a day.

### Application of CVBF on Higgs Boson

* The intention with using the Higgs Boson data set was to apply our method on a large data set where the distributions appear to be the same for some and different for others.
* The data set is huge. It can be downloaded here: https://archive.ics.uci.edu/dataset/280/higgs.
* We took a subset of the data and applied CVBF to some of the columns. The code that applied it (and used the Rpackage) can be seen in `HiggsBosonCol23andCol29withRpackage.rmd`.
* We also used the same data set to compute predictive posteriors of the densities. It's one way to compare our method to Polya Trees. That can be seen in `HiggsBosonRfileWithPackage.rmd`.
* As a bit of a historical relic, versions of this that we made, before using the Rpackage, are also available to see. Diagnostic plots for assessing that good behavior of the MCMC (less auto correlation between points, decent rejection rates) for the predictive posterior (for CVBF) are present. Such a procedure isn't needed for Polya Tree as their predictive posterior is closed (MCMC isn't needed for it). 

