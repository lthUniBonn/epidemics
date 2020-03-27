source('evalModules.R')
library("ggplot2")
library("viridis")
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


sDistFactor <- 3
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 


epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------
#heatmap of immunity and sDistFactor | p(ME) or accInfections
checkThisImmunity <-seq(0,0.5,0.02)
checkThisSDistFac <-round(seq(1,7,0.2),1)#need to round to compare seq[8] == 2.4 to TRUE
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


m <- ggplot(as.data.frame(p), aes(sDistFactors, immunities, fill= prob), ) + 
  geom_tile() +
  scale_fill_gradient2(low = 'white', mid = 'green', high = 'black', midpoint = 0.5, limits=c(0,1)) +
  theme(panel.background = element_blank(), plot.background = element_blank()) +
  labs(x='social distancing D',y= 'immunity') +
  theme(axis.title.x = element_text(color = "black", size = 14, face = "bold"),
        axis.title.y = element_text(color = "black", size = 14, face = "bold"),
        legend.position = c(0.7,0.6)) +
  labs(fill='P(ME)')
  #guides(fill=guide_legend(title="P(ME)"))
print(m)


