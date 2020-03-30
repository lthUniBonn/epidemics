source('evalModules.R')
library("ggplot2")

path <- "LaviBlackie"



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

sDistFactor <- 4
immunity <- 0.4

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
  nInfDf <- nInfDf[seq(1,nrow(nInfDf),1), ]
  R0Df <- array(NA, dim=c(nrow(nInfDf)-1, ncol(nInfDf)-1))

  #get dt values (only take the first, is the same anyways)
  dt <- nInfDf[2,1]
  runCounter <- 0
  for (runIdx in c(1:(ncol(R0Df)))){ 
    nInf <- nInfDf[,runIdx+1]
    if (anyNA(nInf)){
      whichNA <- which(is.na(nInf)==TRUE)[1]
      if (nInf[whichNA-1] != 0) {
        nInf[whichNA] <- 0
        nInf <- nInf[c(1:whichNA)]
      }
      else{
        nInf <- nInf[c(1:(which(is.na(nInf)==TRUE)[1]-1))]
      }
    }
    if (length(nInf)==1){next()}#ignore runs that die out in the first timestep
    runCounter <- runCounter + 1
    #define steps from 0 to 1
    nInf0 <- nInf[c(1:(length(nInf)-1))]
    nInf1 <- nInf[c(2:length(nInf))]
    #calcR0
    R0Df[c(1:length(nInf0)),runCounter] <- mapply(calcR0, nInf0, nInf1, dt=dt) 
  }
  if (anyNA(R0Df[1,])){R0Df <- R0Df[,c(1:(which(is.na(R0Df[1,])==TRUE)[1]-1))]}
  
  #average runs for each timestep
  R0MeanDf <- nInfDf[c(1:nrow(R0Df)),c(1,2,3,4)]#time mean err usedRuns
  R0MeanDf[,1] <- R0MeanDf[,1]+dt/2 #shift time to middle of step for better visualization
  R0MeanDf[,2] <- rowMeans(R0Df, na.rm = TRUE)
  R0MeanDf[,3] <- apply(X=R0Df, MARGIN = 1, FUN=bootstrap)
  
  for(i in c(1:nrow(R0MeanDf))){
    R0MeanDf[i,4] <- length(which(is.na(R0Df[i,])==FALSE))/(ncol(nInfDf)-1)
  }
  return(R0MeanDf)
}

R0MeanDf <- calcR0File(params = params)
ylim = c(0,3)
par(mar = c(3,3,3,3))
plot(R0MeanDf[,c(1,2)], ylim = ylim,, xlab = '', ylab = '', cex = 1, xaxt = 'n', yaxt = 'n')
arrows(R0MeanDf[,1], R0MeanDf[,2]-R0MeanDf[,3], R0MeanDf[,1], R0MeanDf[,2]+R0MeanDf[,3], length=0.05, angle=90, code=3)
#mtext(text = paste("sqrt(N):", fixedParams[1],"  nShort:", fixedParams[2], "  immunity:", fixedParams[3],"  recoveryTime:", fixedParams[4],"+-", fixedParams[5],"  suscDist:", fixedParams[6],"  suscFactor:", fixedParams[7], sep = " "),side = 3)
par(new = TRUE, mar = c(3,3,3,3))
plot(R0MeanDf[,c(1,4)], ylim = ylim, col='red', xlab = '', ylab = '', cex = 1, xaxt = 'n', yaxt = 'n')

title(ylab = 'R0*', line = 1.7, cex.lab = 1.5)
title(xlab = "t", line = 1.7, cex.lab = 1.5)
axis(2, mgp=c(3, .5, 0))
axis(4, mgp=c(3, .5, 0))
axis(1, mgp=c(3, .5, 0))
mtext(side = 4, line = 1.7, 'portion of runs used', col = 'red', cex = 1.5)
mtext(text = paste( "r:", fixedParams[3], sep = " "),side = 3, padj = 3.5, adj = 0.9, col = 'black', cex = 1.5)
mtext(text = paste( "D:", fixedParams[7], sep = " "),side = 3, padj = 5, adj = 0.9, col = 'black', cex = 1.5)
#mtext(text = paste( "%sim", sep = " "),side = 3, padj = 6.5, adj = 0.9, col = 'red', cex = 1.5)
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
