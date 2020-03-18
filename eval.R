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
  R0MeanDf <- read.table(file = paste(c("data/R0",R0Choice,"Mean_", params, ".txt"),sep="", collapse=""))
  # how to treat 0 standard deviation values
  R0MeanDf <- R0MeanDf[-which((R0MeanDf[,2]==0),arr.ind = TRUE),]
  R0Mean <- sum(R0MeanDf[,1]*(R0MeanDf[,2]**(-2)),na.rm = TRUE)/sum(R0MeanDf[,2]**(-2),na.rm = TRUE)
  R0Se <- sqrt(1/sum(R0MeanDf[,2]**(-2),na.rm = TRUE))
  meanAgeDist <- weighted.mean(x = sAgeDist/sDistFactor, w = c(3, 15, 25, 28, 22, 7))
  meanRecoveryTime <- avgRecoveryTime
  varAgeDist <- sum(c(0.03, 0.15, 0.25, 0.28, 0.22, 0.07)*(sAgeDist/sDistFactor-meanAgeDist)**2)
  varRecTime <- sdRecoveryTime**2
  
  R0calc <- (3+2*nShort/N)*meanAgeDist*meanRecoveryTime
  
  R0calcSd <- (3+2*nShort/N)*sqrt(varAgeDist*varRecTime + varAgeDist*meanRecoveryTime**2+ varRecTime*meanAgeDist**2)
  #http://www.odelama.com/data-analysis/Commonly-Used-Math-Formulas/
  plot(x=c(1:19), y=R0MeanDf[,1], ylim =c(0,5))
  arrows(c(1:19), R0MeanDf[,1]-R0MeanDf[,2], c(1:19), R0MeanDf[,1]+R0MeanDf[,2], length=0.05, angle=90, code=3)
  abline(h = (R0Mean+R0Se))
  abline(h = R0Mean)
  abline(h = (R0Mean-R0Se))
  abline(h = R0calc)#!! other color 
  
}


R0Params <- which(nameList == "R02Mean")

#max(infected) -> R0 --> major epidemic 

#max(infected) vs sDistFact --> major epidemic (percolation threshold)
# same w/ immunity

#plot SIR vs time
# accumulatedplots of actual spreading

evalTheseIdx <- which(nameList == "numberInfected")
evalThisParams <- paramList[,evalTheseIdx]
evalThisNames <- nameList[evalTheseIdx]
dfList <- list(length=length(evalTheseIdx))


