library('plot.matrix')
source('modules.R')


set.seed(1)


#-----------------------------run params 

path <- "../data/"

#observables
   #plotting
plotEvery <- 10
plotAccumulated <- F
   #clusters and general eval
checkCluster <- T
clusterEvery <- 5


nStatRun <- 2 # how many times is the same thing done (statistical simulation check)

#-----------------------------lattice 

#network
N <- 400**2 # number of people
nShort <- N/100 # how many shortcuts are created // not impossibly many, could lead to ~inf loop #!!
periodicBoundaries <- TRUE

lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N)))

if(periodicBoundaries){# these are the possible connections in 2D lattice
  possConn <- array(0,dim=c((2*N+nShort),2)) 
} else {
  possConn <- array(0,dim=c((2*N-2*sqrt(N)+nShort),2))
}

#find all possible connections and their bond probability 
#instead of using findConn any other list of bonds can be inserted into possConn
findConn()
nCon <- nrow(possConn)




#age distribution of germany drawn for each individual
ageDistribution <- sample(c(3, 20,40,60,80, 100), replace=TRUE, size=N,prob = c(3, 15, 25, 28, 22, 7))


#-----------------------------epidemic params
immunity <- seq(0,0.5,0.02) #ratio of initially immune people 

#recovery distribution 
avgRecoveryTimeVec <- c(6) #mean 
sdRecoveryTimeVec <- c(2) #sd

#susceptibility distributions
sChoice <- c(F,F,F,T)
sChoiceNames <- c("sBool", "sFixed", "sNot", "sReal")
sBool <- sChoice[1] # if True the susceptibility is 1 or 0 
sFixed <- sChoice[2] 
sNot <- sChoice[3]
sReal <- sChoice[4] #real: used for major simulation, refers to age distribution
#age distribution susceptibility 
sAgeDist1 <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)
sAgeDist2 <- c(0.98, 0.7, 0.45, 0.7, 0.85, 1)
sAgeDistArray <- rbind(sAgeDist1, sAgeDist2) 

#transmissibility set to 1, modified by social distancing factor
sDistFactorVec <- seq(1,7,0.2) # social distancing factor



#----------------------------- start running
#set globals
infected <- logical(length = N)
infectionTime <- numeric(length=N)
recovered <- logical(length=N)
people <- c(1:N)
weight <- numeric(N) 

x <- 0

timeGauge <- length(immunity)*nrow(sAgeDistArray)*length(sdRecoveryTimeVec)*length(avgRecoveryTimeVec)*length(sDistFactorVec)*nStatRun
runTime <- 0

obsNames <- c("x","maxWeight", "numberCluster", "numberInfected", "largeOverTotal","accInfections")

nPar <- 1
nParConfigs <- timeGauge/nStatRun
paramConfigs <- array(data=NA, dim=c(nParConfigs,5))
for(imm in immunity){  
  for(i in c(1:nrow(sAgeDistArray))){
    for(sdRecoveryTime in sdRecoveryTimeVec){
      for(avgRecoveryTime in avgRecoveryTimeVec){ 
        for(sDistFactor in sDistFactorVec){
            paramConfigs[nPar,] <- c( i,sdRecoveryTime,avgRecoveryTime, sDistFactor, imm) 
            nPar <- nPar + 1
          }
      }
    }
  }
}
for (parConf in c(1:nParConfigs)){
  i <- paramConfigs[parConf,1]
  sdRecoveryTime <- paramConfigs[parConf, 2]
  avgRecoveryTime <- paramConfigs[parConf, 3]
  sDistFactor <- paramConfigs[parConf, 4]
  immunity <- paramConfigs[parConf, 5]
    
  sAgeDist <- sAgeDistArray[i,]
  #set susceptibility distribution in population
  if(sBool){
    sDistribution <- sample(c(0,1), replace=TRUE, size=N,prob = c(immunity,1-immunity))
  }
  if(sFixed){
    sDistribution <- sample(c(0,0.25,0.5,0.75,1), replace=TRUE, size=N,prob = c(immunity,(1-immunity)/8,  (1-immunity)/2, (1-immunity)/4, (1-immunity)/8))
  }
  if(sNot){
    sDistribution <- sample(c(0.25,0.5,0.75,1), replace=TRUE, size=N,prob = c((1-immunity)/8,  (1-immunity)/2, (1-immunity)/4, (1-immunity)/8))
  }
  if(sReal){
    
    sDistribution <- unlist(lapply(ageDistribution,FUN = sAge))
    sDistribution[sample(c(1:N), size=N*immunity, replace = FALSE)] <- 0
  }
  
  #save Distribution of initially susceptible people for eval
  initialsDistribution <- sDistribution
  
  #parameters displayed in filename
  params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, i, sDistFactor, sChoiceNames[sChoice]), sep="", collapse="_") 
  print(params)
  
  evalArray <- array(NA, dim = c(10000,nStatRun,6), dimnames = list(c(),c(),obsNames))
  for(statRun in c(1:nStatRun)){
    simulationRun(statRun)
    runTime <- runTime + 1
    print(paste(c(runTime,timeGauge), sep="", collapse = " / "))
  }
  writeEval(evalArray, params = params)
}  


