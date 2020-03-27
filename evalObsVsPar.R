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
#sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#sAgeDist2 <- c(0.9, 0.6, 0.4, 0.6, 0.8, 0.9)
#sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) #susceptibility depending on age


sDistFactor <- 1
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 

epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------

#t till epidemic dies out vs par (heatMap)


evalObs <- function(dfList, obs){
  maxVal <- numeric(length(dfList))
  maxErr <- numeric(length(dfList))
  x <- 0
  if( obs == 'accInfections'){
    thisObs <- outbreakMeasure(dfList, prob = F)
    maxVal <- thisObs[[1]]
    maxErr <- thisObs[[2]]
  }
  else if((obs == 'numberInfected') | (obs == 'numberCluster')){
    for (df in dfList){
      x <- x + 1
      thisObsMax <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
      if(obs == 'numberInfected'){thisObsMax <- thisObsMax/N}#quota of pop
      maxVal[x] <- mean(thisObsMax)
      maxErr[x] <- bootstrap(thisObsMax)
    }
  }
  else if(obs == 'duration'){
    for (df in dfList){
      x <- x + 1
      duration <- numeric(length=(ncol(df)-1))
      for (col in c(2:ncol(df))){
        isNa <- is.na(df[,col])
        if(isNa[nrow(df)]==TRUE){
          whichLast <- which(isNa == TRUE)[1]-1
        } else {
          whichLast <- nrow(df)
        }
        
        duration[col-1] <- df[whichLast,1]
      }
      #browser()
      durationMean <- mean(duration)
      durationErr <- bootstrap(durationMean)
      maxVal[x] <- durationMean
      maxErr[x] <- durationErr
    }
  }
  else{
    print('not specified what to do')
    stop()
  }
  return(list(maxVal, maxErr))
}

parIdx <- 7
if(parIdx == 7){
  xlab <- 'social distancing factor D'
} else if(parIdx == 3){
  xlab <- 'immunity R'
}

#evaluate observable to plot -> y
#ylab <- 'total infected'
#obs <- c('accInfections','accInfections')


ylab <- expression('maximal concurrently infected I'['max'])
obs <- c('numberInfected', 'numberInfected')
 
#ylab <- 'number of separate clusters'
#obs <- c('numberCluster','numberCluster')


#ylab <- 'duration'
#obs <- c('accInfections', 'duration')

tmpList <- findObsvsParams(obs=obs[1],parIdx = parIdx, params = fixedParams)[]
dfList <- tmpList[[2]]
parList <- tmpList[[1]]
y <- evalObs(dfList, obs = obs[2])
#corresponding parameter values -> x
parVal <- parList[parIdx,]
yVal <- y[[1]]
yErr <- y[[2]]


plot(x = parVal, y = yVal, xlab = "", ylab = "", pch = 19, cex = 1,xaxt = 'n', yaxt = 'n')
title(ylab = ylab, line = 1.7, cex.lab = 1.5)
title(xlab = parNames[parIdx], line = 1.7, cex.lab = 1.5)
axis(2, mgp=c(3, .5, 0))
axis(1, mgp=c(3, .5, 0))
arrows(parVal, yVal-yErr, parVal, yVal+yErr, length=0.05, angle=90, code=3)





