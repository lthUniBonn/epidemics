source('evalModules.R')



#get the existing filenames
fileNames <- array(unlist(strsplit(list.files("data")[], "_")), dim=c(9,length(list.files("data"))))
paramList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files("data"))))
nameList <- fileNames[1,]

#read all specified files
evalTheseIdx <- which(nameList == "numberInfected")
evalThisParams <- paramList[,evalTheseIdx]
evalThisNames <- nameList[evalTheseIdx]
dfList <- list(length=length(evalTheseIdx))
for (idx in c(1:length(evalTheseIdx))){
  params <- paste(c(evalThisParams[1,idx], evalThisParams[2,idx], evalThisParams[3,idx], evalThisParams[4,idx], evalThisParams[5,idx], evalThisParams[6,idx], evalThisParams[7,idx], sChoiceNames[sChoice]), sep="", collapse="_") 
  df <- read.table(file = paste(c("data/", evalThisNames[idx], "_", params, ".txt"),sep="", collapse=""))
  dfList[[idx]] <- df
}


#plotting / evaluation
a <- c(1,2)
comp <- T
meanPlot(evalThisNames[a[1]], evalThisParams[,a[1]], dfList[[a[1]]], compare = comp, evalThisNames[a[2]], evalThisParams[,a[2]], dfList[[a[2]]])



#R0

evalR0 <- function(R0Choice = "", params){
  #!! check params immunity <- params[1]
  R0MeanDf <- read.table(file = paste(c("data/R0",R0Choice,"Mean_", params, ".txt"),sep="", collapse=""))
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

findObsvsParams <- function(obs='numberInfected', parIdx=3, params){
  #browser()
  # find all files with obs 
  parList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files("data"))))
  obsList <- fileNames[1,]
  parList <- paramList[,which(obsList == obs)]
  notParIdxVec <- c(1:7)
  notParIdxVec <- notParIdxVec[-which(notParIdxVec == parIdx)]

  #select param config, all but one fixed
  fixedParListCheck <- rep(TRUE, ncol(parList))
  for (idx in notParIdxVec){
    fixedParListCheck <- fixedParListCheck & (parList[idx,]==params[idx])
  }
  parList <- array(parList[,which(fixedParListCheck==TRUE)], dim = c(7, length(which(fixedParListCheck==TRUE))))
  #read files into list of data frames 
  dfList <- list(length=ncol(parList))
  for(idx in c(1:ncol(parList))){
    readParams <- paste(c(parList[1,idx], parList[2,idx], parList[3,idx], parList[4,idx], parList[5,idx], parList[6,idx], parList[7,idx], sChoiceNames[sChoice]), sep="", collapse="_") 
    df <- read.table(file = paste(c("data/", obs, "_", readParams, ".txt"),sep="", collapse=""))
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
outbreakProb <- function(dfList){
  p <- numeric(length(dfList))
  x <- 0
  for (df in dfList){
    x <- x + 1
    lastVal <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
    p[x] <- length(which(lastVal>epidemicThreshold))/(ncol(df)-1)
    
  } 
  return(p)
}
a <- outbreakProb(dfList)
plot(parVal, a)
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


