source('evalModules.R')
library("ggplot2")

path <- "longRun"

#define basic parameters of the evaluated simulation
N <- 400**2
nShort <- N/100
avgRecoveryTime <- 6
sdRecoveryTime <- 2
sChoice <- 'sReal'

ageDistIdx <- 1
#sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
#sAgeDist2 <- c(0.9, 0.6, 0.4, 0.6, 0.8, 0.9)
#sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) #susceptibility depending on age


sDistFactor <- 1#3
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)
#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice), sep="", collapse="_") 


epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------

#plot means of obs vs t
obs <- 'numberInfected' # 'accInfections'
comp <- T
params <- fixedParams 
df <- read.table(paste(c(path,"/",obs, '_',paste(params, sep="", collapse="_"), '.txt'), sep="", collapse="")) 
if (comp){
  params2 <- params
  params2[7] <- 3
  df2 <- read.table(paste(c(path,"/",obs, '_',paste(params2, sep="", collapse="_"), '.txt'), sep="", collapse="")) 
  meanPlot('Infected I(t)', params, df, compare = comp, name2 = obs, params2, df2)
} else{meanPlot(obs, params, df, compare = comp)}
  
