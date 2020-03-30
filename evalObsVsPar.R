source('evalModules.R')
library("ggplot2")

#path <- "longRun"
path <- "forgottenData"


#define basic parameters
N <- 400**2
nShort <- N/100
avgRecoveryTime <- 6
sdRecoveryTime <- 2
sChoice <- 'sReal'

ageDistIdx <- 2
#sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#sAgeDist2 <- c(0.9, 0.6, 0.4, 0.6, 0.8, 0.9)
#sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) #susceptibility depending on age


sDistFactor <- 1
immunity <- 0.2
sDistFactor2 <- 3
immunity2 <- 0.0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
fixedParams2 <- c(sqrt(N), nShort, immunity2, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor2, sChoice)
#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 

epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------

#t till epidemic dies out vs par (heatMap)


evalObs <- function(dfList, obs){
  maxVal <- numeric(length(dfList))
  maxErr <- numeric(length(dfList))
  x <- 0
  if( obs == 'accInfections'){
    thisObs <- outbreakMeasure(dfList, prob = T)
    maxVal <- thisObs[[1]]
    maxErr <- thisObs[[2]]
  }
  else if((obs == 'numberInfected') | (obs == 'numberCluster') | (obs == 'maxWeight')){
    for (df in dfList){
      x <- x + 1
      thisObsMax <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
      #if(obs == 'numberInfected'){thisObsMax <- thisObsMax/N}#quota of pop
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
  xlab <- 'D'
  xlim <- c(1,7)
  constant1 <- paste('r:', fixedParams[3], sep = "", collapse = "")
  constant2 <- paste('r:', fixedParams2[3], sep = "", collapse = "")
} else if(parIdx == 3){
  xlab <- 'r'
  xlim <- c(0,0.5)
  constant1 <- paste('D:', fixedParams[7], sep = "", collapse = "")
  constant2 <- paste('D:', fixedParams2[7], sep = "", collapse = "")
}

#evaluate observable to plot -> y
ylab <- expression('P'[ME])
obs <- c('accInfections','accInfections')


#ylab <- expression('I'['max'])
#obs <- c('numberInfected', 'numberInfected')
 
#ylab <- 'number of separate clusters'
#obs <- c('numberCluster','numberCluster')


#ylab <- 'duration'
#obs <- c('accInfections', 'duration')

#ylab <- 'max maxWeight'
#obs <- c('maxWeight', 'maxWeight')

tmpList <- findObsvsParams(obs=obs[1],parIdx = parIdx, params = fixedParams)[]
dfList <- tmpList[[2]]
parList <- tmpList[[1]]
y <- evalObs(dfList, obs = obs[2])
#corresponding parameter values -> x
parVal <- parList[parIdx,]
yVal <- y[[1]]
yErr <- y[[2]]

tmpList2 <- findObsvsParams(obs=obs[1],parIdx = parIdx, params = fixedParams2)[]
dfList2 <- tmpList2[[2]]
parList2 <- tmpList2[[1]]
y2 <- evalObs(dfList2, obs = obs[2])
#corresponding parameter values -> x
parVal2 <- parList2[parIdx,]
yVal2 <- y2[[1]]
yErr2 <- y2[[2]]


plot(x = parVal, y = yVal, xlab = "", ylab = "", xlim = xlim,ylim = c(0,1), pch = 19, cex = 1,xaxt = 'n', yaxt = 'n')
arrows(parVal, yVal-yErr, parVal, yVal+yErr, length=0.05, angle=90, code=3)
title(ylab = ylab, line = 1.7, cex.lab = 1.5)
title(xlab = xlab, line = 1.7, cex.lab = 1.5)
axis(2, mgp=c(3, .5, 0))
axis(1, mgp=c(3, .5, 0))
par(new=TRUE)
plot(x = parVal2, y = yVal2, xlab = "", ylab = "", xlim = xlim, ylim = c(0,1), pch = 19, cex = 1,xaxt = 'n', yaxt = 'n', col = 'red')

mtext(text = constant1, side = 3, padj = 3.5, adj = 0.9, cex = 1.5)
mtext(text = constant2, side = 3, padj = 2, adj = 0.9, cex = 1.5, col = 'red')

arrows(parVal2, yVal2-yErr2, parVal2, yVal2+yErr2, length=0.05, angle=90, code=3, col = 'red')





