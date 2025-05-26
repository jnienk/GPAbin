###################################################
## Commands for: GPAbin: unifying visualizations of multiple mputations for 
## missing value
## Communications in Statistics: Simulation and Computation
## Authors: J Nienkemper-Swanepoel, NJ le Roux & S Gardner-Lubbe
## Corresponding author: J Nienkemper-Swanepoel
## Centre for Multi-dimensional data visualisation (MuViSU), 
## Department of Statistics and Actuarial Science, Stellenbosch University.
## 2023
###################################################

###required libraries:
library(missMDA)
library(FactoMineR)
library(ca)

###example data sets:

comp.dat <- read.table(file="Data\\comp.dat.txt",sep=",")
miss.dat <- read.table(file="Data\\miss.dat.txt",sep=",")

###functions:
#biplFig(), CLPred(), compMeas(), delCL(), df2fact(),
#FormatDat(), FormatDimNam(), FormatImpList(), GPA(), GPAbin(),
#indcol(), indmat(), is.integer0(), MIimpute(), myestim_ncpMCA(),
#myimputeMCA(), myMIMCA(), myOPA(), rmOneCL()

source("Code\\GPAbinfunctions.R")

######################################FUNCTION CALLS######################################

##Multiple imputation and construction of MCA biplots for each imputation
#depending on computer capactiy, should compute for a few minutes, 
#due to the estimation of the number of dimensions for the reconstruction of MIMCA

MI.out <- MIimpute(datNA=miss.dat,seed=123,imps=3)

##aligning multiple imputation biplots and calculation of GPAbin biplot

GPAbin.out <- GPAbin(CLP.list=MI.out[[1]],Z.list=MI.out[[2]])

##MCA of complete data

mca.comp <- mjca(comp.dat, indicator="lambda")

##visual interpretations in 2D

pred.comp <- CLPred(comp.lvls=as.matrix(mca.comp[[6]]),datIN=comp.dat,CLPs=mca.comp[[23]],Zs=mca.comp[[16]],nsamples=100,pvar=5)
pred.imp <- CLPred(comp.lvls=as.matrix(mca.comp[[6]]),datIN=miss.dat,CLPs=GPAbin.out[[2]],Zs=GPAbin.out[[1]],nsamples=100,pvar=5)

##measures of comparison
measures <- compMeas(Target=mca.comp[[23]], Testee=GPAbin.out[[2]], dim="2D", pred.dat1=pred.comp, pred.dat2=pred.imp,nsamples=100,pvar=5)

#construction of biplots (selected illustrations)
biplFig(CLPs=mca.comp[[23]], Zs=mca.comp[[16]], Lvls=mca.comp[[6]], title="Fig1")
biplFig(CLPs=GPAbin.out[[2]], Zs=GPAbin.out[[1]], title="Fig5")

#real application Figure 12
data(TitanicNA)
real.dat <- TitanicNA
MI.real <- MIimpute(datNA=real.dat,seed=123,imps=10)
GPAbin.real <- GPAbin(CLP.list=MI.real[[1]],Z.list=MI.real[[2]])
biplFig(CLPs=GPAbin.real[[2]],Zs=GPAbin.real[[1]],title="Fig12")