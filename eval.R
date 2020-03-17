
N <- 400**2
nShort <- N/100
immunity <- 0
avgRecoveryTime <- 3
sdRecoveryTime <- 1
i <- 1
sAgeDist0 <- c(1, 0.7, 0.5, 0.7, 0.8, 1)
sDistFactor <- 3
sChoice <- c(F,F,F,T)
sChoiceNames <- c("sBool", "sFixed", "sNot", "sReal")
sBool <- sChoice[1] # if True the susceptibility is 1 or 0 // other poss like gaussian with age etc
sFixed <- sChoice[2]
sNot <- sChoice[3]
sReal <- sChoice[4]

params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, i, sDistFactor, sChoiceNames[sChoice]), sep="", collapse="_") 
print(params)

#read data from file 
maxWeightDf <- read.table(file = paste(c("data/maxWeight_", params, ".txt"),sep="", collapse=""))
# numberClusterDf <- read.table(file = paste(c("data/numberCluster_", params, ".txt"),sep="", collapse=""))
# numberInfectedDf <- read.table(file = paste(c("data/numberInfected_", params, ".txt"),sep="", collapse=""))
# largeOverTotalDf <- read.table(file = paste(c("data/largeOverTotal_", params, ".txt"),sep="", collapse=""))
# accInfectionsDf <- read.table(file = paste(c("data/accInfections_", params, ".txt"),sep="", collapse=""))


maxWeightMean <- rowMeans(maxWeightDf[,c(2:ncol(maxWeightDf))], na.rm = TRUE)
maxWeightSd <- apply(X = maxWeightDf[,c(2:ncol(maxWeightDf))],MARGIN = 1, FUN = sd, na.rm =TRUE)


plot(maxWeightDf[,1], maxWeightMean, xlab = "T",ylab = "MaxWeightMean")
arrows(maxWeightDf[,1], maxWeightMean-maxWeightSd, maxWeightDf[,1], maxWeightMean+maxWeightSd, length=0.05, angle=90, code=3)

R0MeanDf <- read.table(file = paste(c("data/R04Mean_", params, ".txt"),sep="", collapse=""))
R0Mean <- sum(R0MeanDf[,1]*(R0MeanDf[,2]**(-2)),na.rm = TRUE)/sum(R0MeanDf[,2]**(-2),na.rm = TRUE)
R0Se <- sqrt(1/sum(R0MeanDf[,2]**(-2),na.rm = TRUE))
R0calc <- (3+2*nShort/N)*weighted.mean(x = sAgeDist0/sDistFactor, w = c(3, 15, 25, 28, 22, 7))*avgRecoveryTime
R0calcSd <- #http://www.odelama.com/data-analysis/Commonly-Used-Math-Formulas/

plot(numberClusterDf[,1],numberClusterDf[,2], xlab = "T",ylab = "Number of Clusters")
plot(numberInfectedDf[,1],numberInfectedDf[,2], xlab = "T",ylab = "Number Infected")
plot(largeOverTotalDf[,1],largeOverTotalDf[,2], xlab = "T",ylab = "Largest Cluster over Total Infected")
plot(accInfectionsDf[,1],accInfectionsDf[,2], xlab = "T",ylab = "Attack Rate")
plot(R0MeanDf[,1],R0MeanDf[,2], xlab = "T",ylab = "R0 Mean")
#!! maybe add a line for R0calc and lines for its variances