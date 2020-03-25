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
calcR0 <- function(nInf0, nInf1, dt=1){
  if ((nInf0 == 0) | (dt == 0)){
    stop()
  }
  else {
    fractChange <- nInf1/nInf0#1 angesteckt, 1 recovered: fractChange = 1 -->  R0 = 1 sollte stimmen! (nInf1-nInf0)/nInf0#???
    R0 <- fractChange**(1/dt)
    return(R0)
  }
}

calcR0File <- function(params){
  #get file numberInfected
  nInfDf <- read.table(file = paste(c(path,"/","numberInfected","_", params, ".txt"),sep="", collapse=""))
  R0Df <- array(NA, dim=c(nrow(nInfDf)-1, ncol(nInfDf)-1))

  #get dt values (only take the first, is the same anyways)
  dt <- nInfDf[2,1]
  for (runIdx in c(1:(ncol(R0Df)))){ 
    nInf <- nInfDf[,runIdx+1]
    if (anyNA(nInf)){nInf <- nInf[c(1:(which(is.na(nInf)==TRUE)[1]-1))]}
    #define steps from 0 to 1
    nInf0 <- nInf[c(1:(length(nInf)-1))]
    nInf1 <- nInf[c(2:length(nInf))]
    #calcR0
    R0Df[c(1:length(nInf0)),runIdx] <- mapply(calcR0, nInf0, nInf1, dt=dt) 
  }
  #average runs for each timestep
  R0MeanDf <- nInfDf[c(1:nrow(R0Df)),c(1,2,3)]#time mean err 
  R0MeanDf[,2] <- rowMeans(R0Df, na.rm = TRUE)
  R0MeanDf[,3] <- apply(X=R0Df, MARGIN = 1, FUN=bootstrap)
  return(R0MeanDf)
}

R0MeanDf <- calcR0File(params = params)

#----------------------------------------

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

#evalR0(R0Choice = "2",params=params)
