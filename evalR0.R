source('evalModules.R')
library("ggplot2")

path <- "longRun"



#define basic parameters
N <- 400**2
nShort <- N/100
avgRecoveryTime <- 6
sdRecoveryTime <- 2
sChoice <- 'sReal'

ageDistIdx <- 1
sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
sAgeDist2 <- c(0.9, 0.6, 0.4, 0.6, 0.8, 0.9)
sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) #susceptibility depending on age
sAgeDist <- sAgeDistArray[ageDistIdx,]

sDistFactor <- 1
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 

epidemicThreshold <- 0.02

#--------------

evalR0 <- function(R0Choice = "", params){
  #!! check params immunity <- params[1]
  R0MeanDf <- read.table(file = paste(c(path,"/","R0",R0Choice,"Mean_", params, ".txt"),sep="", collapse=""))
  # how to treat 0 standard deviation values
  if (length(which((R0MeanDf[,2]==0)))>0){R0MeanDf <- R0MeanDf[-which((R0MeanDf[,2]==0),arr.ind = TRUE),]}
  R0Mean <- sum(R0MeanDf[,1]*(R0MeanDf[,2]**(-2)),na.rm = TRUE)/sum(R0MeanDf[,2]**(-2),na.rm = TRUE)
  R0Se <- sqrt(1/sum(R0MeanDf[,2]**(-2),na.rm = TRUE))
  meanAgeDist <- weighted.mean(x = sAgeDist/sDistFactor, w = c(3, 15, 25, 28, 22, 7))
  meanRecoveryTime <- avgRecoveryTime
  varAgeDist <- sum(c(0.03, 0.15, 0.25, 0.28, 0.22, 0.07)*(sAgeDist/sDistFactor-meanAgeDist)**2)
  varRecTime <- sdRecoveryTime**2
  
  R0calc <- (3+2*nShort/N)*meanAgeDist*meanRecoveryTime*(1-immunity)
  
  R0calcSd <- (3+2*nShort/N)*(1-immunity)*sqrt(varAgeDist*varRecTime + varAgeDist*meanRecoveryTime**2+ varRecTime*meanAgeDist**2)
  #http://www.odelama.com/data-analysis/Commonly-Used-Math-Formulas/
  plot(x=c(1:nrow(R0MeanDf)), y=R0MeanDf[,1], ylim =c(0,5))
  arrows(c(1:nrow(R0MeanDf)), R0MeanDf[,1]-R0MeanDf[,2], c(1:nrow(R0MeanDf)), R0MeanDf[,1]+R0MeanDf[,2], length=0.05, angle=90, code=3)
  abline(h = (R0Mean+R0Se))
  abline(h = R0Mean)
  print(R0calc)
  abline(h = (R0Mean-R0Se))
  abline(h = R0calc)#!! other color 
  
}

evalR0(R0Choice = "2",params=params)
