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
#eval max vs immunity / social distancing 

#max number cluster vs social distancing 
#max number cluster vs immunity

#obs vs param

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

