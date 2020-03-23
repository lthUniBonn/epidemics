source('evalModules.R')
library("ggplot2")

path <- "data"



#define basic parameters
N <- 400**2
nShort <- N/100
avgRecoveryTime <- 6
sdRecoveryTime <- 2
sChoice <- 'sReal'

ageDistIdx <- 1
#sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#sAgeDist2 <- c(0.9, 0.6, 0.4, 0.6, 0.8, 0.9)
#sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) #susceptibility depending on age


sDistFactor <- 3
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 


epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------
#plot obs vs var


#heatmap of immunity and sDistFactor | p(ME) or accInfections
#R0calc vs R0 measured 
  #bootstrapping for error?!
#(R0calc und p(ME)) vs sDistFactor

#-----------------------
#get the existing filenames
fileNames <- array(unlist(strsplit(list.files(path)[], "_")), dim=c(9,length(list.files(path))))
paramList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files(path))))
nameList <- fileNames[1,]
#alt

#a <- readAllObs(nStat=50)


#plotting mean vs x 
comp <- F
params <- fixedParams #set these in evalModules.R
obs <- 'accInfections'
df <- read.table(paste(c(path,"/",obs, '_',paste(params, sep="", collapse="_"), '.txt'), sep="", collapse="")) 
meanPlot(obs, params, df, compare = comp)

#plotting 
a <- c(1,2)

#meanPlot(evalThisNames[a[1]], evalThisParams[,a[1]], dfList[[a[1]]], compare = comp, evalThisNames[a[2]], evalThisParams[,a[2]], dfList[[a[2]]])



#R0

evalR0 <- function(R0Choice = "", params){
  #!! check params immunity <- params[1]
  R0MeanDf <- read.table(file = paste(c(path,"/","R0",R0Choice,"Mean_", params, ".txt"),sep="", collapse=""))
  # how to treat 0 standard deviation values
  R0MeanDf <- R0MeanDf[-which((R0MeanDf[,2]==0),arr.ind = TRUE),]
  R0Mean <- sum(R0MeanDf[,1]*(R0MeanDf[,2]**(-2)),na.rm = TRUE)/sum(R0MeanDf[,2]**(-2),na.rm = TRUE)
  R0Se <- sqrt(1/sum(R0MeanDf[,2]**(-2),na.rm = TRUE))
  meanAgeDist <- weighted.mean(x = sAgeDist/sDistFactor, w = c(3, 15, 25, 28, 22, 7))
  meanRecoveryTime <- avgRecoveryTime
  varAgeDist <- sum(c(0.03, 0.15, 0.25, 0.28, 0.22, 0.07)*(sAgeDist/sDistFactor-meanAgeDist)**2)
  varRecTime <- sdRecoveryTime**2
  
  R0calc <- (3+2*nShort/N)*meanAgeDist*meanRecoveryTime*(1-immunity)
  
  R0calcSd <- (3+2*nShort/N)*(1-immunity)*sqrt(varAgeDist*varRecTime + varAgeDist*meanRecoveryTime**2+ varRecTime*meanAgeDist**2)
  #http://www.odelama.com/data-analysis/Commonly-Used-Math-Formulas/
  plot(x=c(1:19), y=R0MeanDf[,1], ylim =c(0,5))
  arrows(c(1:19), R0MeanDf[,1]-R0MeanDf[,2], c(1:19), R0MeanDf[,1]+R0MeanDf[,2], length=0.05, angle=90, code=3)
  abline(h = (R0Mean+R0Se))
  abline(h = R0Mean)
  abline(h = (R0Mean-R0Se))
  abline(h = R0calc)#!! other color 
  
}

findObsvsParams <- function(obs='numberInfected', parIdx=3, params, checkSpecific=c()){
  #browser()
  # find all files with obs 
  parList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files(path))))
  obsList <- fileNames[1,]
  parList <- paramList[,which(obsList == obs)]
  #browser()
  
  if (length(checkSpecific) == 0){
    notParIdxVec <- c(1:7)
    notParIdxVec <- notParIdxVec[-which(notParIdxVec == parIdx)]
    
    #select param config, all but one fixed
    fixedParListCheck <- rep(TRUE, ncol(parList))
    for (idx in notParIdxVec){
      fixedParListCheck <- fixedParListCheck & (parList[idx,]==params[idx])
    }
    parList <- array(parList[,which(fixedParListCheck==TRUE)], dim = c(7, length(which(fixedParListCheck==TRUE))))
  }
  else {
    notParIdxVec <- c(1:7)
    notParIdxVec <- notParIdxVec[-which(notParIdxVec == parIdx)]
    #select param config, all but one fixed
    fixedParListCheck <- rep(TRUE, ncol(parList))
    for (idx in notParIdxVec){
      fixedParListCheck <- fixedParListCheck & (parList[idx,]==params[idx])
    }
    fixedParListCheck <- fixedParListCheck & (parList[parIdx,] %in% checkSpecific)
    print(length(which(fixedParListCheck == TRUE)))
    parList <- array(parList[,which(fixedParListCheck==TRUE)], dim = c(7, length(which(fixedParListCheck==TRUE))))
  }
  #read files into list of data frames 
  dfList <- list(length=ncol(parList))
  for(idx in c(1:ncol(parList))){
    readParams <- paste(c(parList[1,idx], parList[2,idx], parList[3,idx], parList[4,idx], parList[5,idx], parList[6,idx], parList[7,idx], sChoice), sep="", collapse="_") 
    df <- read.table(file = paste(c(path,"/", obs, "_", readParams, ".txt"),sep="", collapse=""))
    dfList[[idx]] <- df
  }
  
  #meanPlot(name = 'numberInfected', params = parList[,1], df = dfList[[1]], compare = TRUE, name2 = 'numberInfected', params2 = parList[,3], df2 = dfList[[3]])
  
  return(list(parList, dfList))
}
evalMax <- function(dfList, obs){
  maxVal <- numeric(length(dfList))
  maxSd <- numeric(length(dfList))
  x <- 0
  if( obs == 'accInfections'){
    for (df in dfList){
      x <- x + 1
      lastVal <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
      thisObsMean <- mean(lastVal)#rowMeans(df[,c(2:ncol(df))], na.rm = TRUE)
      thisObsSd <- sd(lastVal)#apply(X = df[,c(2:ncol(df))],MARGIN = 1, FUN = sd, na.rm =TRUE)
      if (anyNA(thisObsSd))(View(df[,c(2:ncol(df))]))
      maxVal[x] <- max(thisObsMean)
      maxSd[x] <- thisObsSd[which.max(thisObsMean)]
    }  
  }
  else if(obs == 'largeOverTotal'){
    for (df in dfList){
      x <- x + 1
      thisObsMean <- rowMeans(df[,c(2:ncol(df))], na.rm = TRUE)
      thisObsSd <- apply(X = df[,c(2:ncol(df))],MARGIN = 1, FUN = sd, na.rm =TRUE)
      maxVal[x] <- max(thisObsMean)
      maxSd[x] <- thisObsSd[which.max(thisObsMean)]
    }
  }
  else{
  for (df in dfList){
    x <- x + 1
    thisObsMean <- rowMeans(df[,c(2:ncol(df))], na.rm = TRUE)
    thisObsSd <- apply(X = df[,c(2:ncol(df))],MARGIN = 1, FUN = sd, na.rm =TRUE)
    if (anyNA(thisObsSd))(View(df[,c(2:ncol(df))]))
    maxVal[x] <- max(thisObsMean)
    maxSd[x] <- thisObsSd[which.max(thisObsMean)]
  }
  }
  return(array(data=c(maxVal, maxSd), dim=c(length(dfList),2)))
}

plotObsvsParam <- function(){
  
}
#param names accInfections
parIdx <- 3
tmpList <- findObsvsParams(obs='accInfections',parIdx = parIdx, params = paramList[,5])[]
dfList <- tmpList[[2]]
parList <- tmpList[[1]]
maxVal <- evalMax(dfList, obs = 'accInfections')
parVal <- parList[parIdx,]

plot(x = parVal, y = maxVal[,1])
arrows(parVal, maxVal[,1]-maxVal[,2], parVal, maxVal[,1]+maxVal[,2], length=0.05, angle=90, code=3)
R0Params <- which(nameList == "R02Mean")


#define major epidemic: over 100 people infected = 100/160000

#plot prob of outbreak:
outbreakProb <- function(dfList, prob=T){
  p <- numeric(length(dfList))
  lastValMean <- numeric(length(dfList))
  #x <- 0
  for (dfIdx in c(1:length(dfList))){
    #x <- x + 1
    #print(df)
    df <- dfList[[dfIdx]]
    x <- dfIdx
    lastVal <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
    lastValMean[x] <- mean(lastVal)
    p[x] <- length(which(lastVal>epidemicThreshold))/(ncol(df)-1)
    
  } 
  if(prob){
    return(p)
  } else { return(lastValMean)}
}


checkThisImmunity <- unique(c(seq(0,0.7,0.1)))#,seq(0,0.9,0.1)))
checkThisSDistFac <- unique(c(seq(1,10,1)))#,seq(1,10,1)))


p <- expand.grid(immunities = checkThisImmunity, sDistFactors = checkThisSDistFac)
p$prob <- NA
counter <- 0
for(x in c(1:length(checkThisImmunity))){
  # sqrt(N) nShort immunity avgRecTime sdRecTime sAgeDist  sDistFactor sChoice
  params <-  c(sqrt(N), nShort, checkThisImmunity[x], avgRecoveryTime, sdRecoveryTime, i, sDistFactor, sChoice)
  #find files
  tmpList <- findObsvsParams(obs='accInfections',parIdx = 7, params = params, checkSpecific= checkThisSDistFac)[]
  dfList <- tmpList[[2]]
  parList <- tmpList[[1]] 
  outbreakP <- outbreakProb(dfList, F)
  for (idx in c(1:ncol(parList))){  
    counter <- counter + 1
    p[(p$sDistFactors == parList[7,idx])&(p$immunities == checkThisImmunity[x]),]$prob <- outbreakP[idx]
  }
}

ggplot(data = as.data.frame(p), mapping = aes(x=sDistFactors, y=immunities  )) + geom_tile(aes(fill = prob)) 

#plot(a, ylab = "immunity", xlab = "social distancing")

#function plot observable against selected parameter

#choose some params 


#finish R0 eval 
#epidemic prob vs par  
#

#
#max(infected) -> R0 --> major epidemic 

#max(infected) vs sDistFact --> major epidemic (percolation threshold)
# same w/ immunity

#plot SIR vs time
# accumulatedplots of actual spreading

evalTheseIdx <- which(nameList == "numberInfected")
evalThisParams <- paramList[,evalTheseIdx]
evalThisNames <- nameList[evalTheseIdx]
dfList <- list(length=length(evalTheseIdx))


