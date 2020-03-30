source('evalModules.R')
library("ggplot2")

path <- "../data"

#define basic parameters of the evaluated simulation
N <- 400**2
nShort <- N/100
avgRecoveryTime <- 6
sdRecoveryTime <- 2
sChoice <- 'sReal'

ageDistIdx <- 1


sDistFactor <- 1
immunity <- 0

fixedParams <- c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, ageDistIdx, sDistFactor, sChoice)


epidemicThreshold <- 0.02
#-------------------------------------------------------------------------------

#plot means of obs vs t
obs <- 'numberInfected' # 'accInfections'
comp <- T
params <- fixedParams 
df <- read.table(paste(c(path,"/",obs, '_',paste(params, sep="", collapse="_"), '.txt'), sep="", collapse="")) 
if (comp){
  params2 <- params
  params2[3] <- 0.2
  df2 <- read.table(paste(c(path,"/",obs, '_',paste(params2, sep="", collapse="_"), '.txt'), sep="", collapse="")) 
  meanPlot('I(t)', params, df, compare = comp, name2 = obs, params2, df2)
} else{meanPlot(obs, params, df, compare = comp)}
  
