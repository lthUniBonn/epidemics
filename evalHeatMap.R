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
#heatmap of immunity and sDistFactor | p(ME) or accInfections
checkThisImmunity <- unique(c(seq(0,0.7,0.1)))#,seq(0,0.9,0.1)))
checkThisSDistFac <- unique(c(seq(1,10,1)))#,seq(1,10,1)))

p <- expand.grid(immunities = checkThisImmunity, sDistFactors = checkThisSDistFac)
p$prob <- NA
counter <- 0
for(x in c(1:length(checkThisImmunity))){
  params <-  c(sqrt(N), nShort, checkThisImmunity[x], avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
  #find files
  tmpList <- findObsvsParams(obs='accInfections',parIdx = 7, params = params, checkSpecific= checkThisSDistFac)[]
  dfList <- tmpList[[2]]
  parList <- tmpList[[1]] 
  outbreakP <- outbreakMeasure(dfList, prob = FALSE)[[1]]
  for (idx in c(1:ncol(parList))){  
    counter <- counter + 1
    p[(p$sDistFactors == parList[7,idx])&(p$immunities == checkThisImmunity[x]),]$prob <- outbreakP[idx]
  }
}

ggplot(data = as.data.frame(p), mapping = aes(x=sDistFactors, y=immunities  )) + geom_tile(aes(fill = prob)) 



