# https://www.healthknowledge.org.uk/public-health-textbook/research-methods/1a-epidemiology/epidemic-theory
library('plot.matrix')
source('modules.R')
source('runFunctions.R')

#profile <- profvis({

#set.seed(1)




#-----------------------------run params 
path <- "data/"

#observables
   #plotting
plotIt <- FALSE
plotEvery <- 10
plotAccumulated <- FALSE
   #clusters
checkCluster <- TRUE
clusterEvery <- 5


nStatRun <- 50 # how many times is the same thing done (statistical simulation check)
checkR0Here <- 10 # after how many recoveries is R0 evaluated 

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
findConn()
nCon <- nrow(possConn)





ageDistribution <- sample(c(3, 20,40,60,80, 100), replace=TRUE, size=N,prob = c(3, 15, 25, 28, 22, 7))

#-----------------------------epidemic params



#vaccination #!!
immunity <- 0.3 #ratio of immune people 

#recovery
avgRecoveryTimeVec <- c(4)
sdRecoveryTimeVec <- c(1)

#susceptibility
sChoice <- c(F,F,F,T)
sChoiceNames <- c("sBool", "sFixed", "sNot", "sReal")
sBool <- sChoice[1] # if True the susceptibility is 1 or 0 // other poss like gaussian with age etc
sFixed <- sChoice[2]
sNot <- sChoice[3]
sReal <- sChoice[4]

sAgeDist0 <- c(1, 0.7, 0.5, 0.7, 0.8, 1)
sAgeDist1 <- c(1, 0.7, 0.7, 0.7, 0.8, 1)
sAgeDist2 <- c(2, 2, 2, 2, 2, 2)
sAgeDistArray <- rbind(sAgeDist0)#,sAgeDist1,sAgeDist2) #susceptibility depending on age
#transmissibility
sDistFactorVec <- seq(1,20, 0.1) # social distancing factor



#----------------------------- start running
#set globals
R0OverInfectiousPeriod <- numeric(length=N)
R0Mean <- numeric()
R0Sd <- numeric()
R0Mean2 <- numeric()
R0Sd2 <- numeric()
R0Mean3 <- numeric()
R0Sd3 <- numeric()
R0Mean4 <- numeric()
R0Sd4 <- numeric()
R0Mean7 <- numeric()
R0Sd7 <- numeric()

infected <- logical(length = N)
infectionTime <- numeric(length=N)
recovered <- logical(length=N)
people <- c(1:N)
weight <- numeric(N) 

x <- 0

timeGauge <- nrow(sAgeDistArray)*length(sdRecoveryTimeVec)*length(avgRecoveryTimeVec)*length(sDistFactorVec)*nStatRun
runTime <- 0

obsNames <- c("x","maxWeight", "numberCluster", "numberInfected", "largeOverTotal","accInfections", "R0Mean")
#!! write a list from in which every param is listed, then it is possible to go through the list, stop and continue if something goes wrong

nPar <- 1
nParConfigs <- timeGauge/nStatRun
paramConfigs <- array(data=NA, dim=c(nParConfigs,4))
for(i in c(1:nrow(sAgeDistArray))){
  for(sdRecoveryTime in sdRecoveryTimeVec){
    for(avgRecoveryTime in avgRecoveryTimeVec){ 
      for(sDistFactor in sDistFactorVec){
          paramConfigs[nPar,] <- c(i,sdRecoveryTime,avgRecoveryTime, sDistFactor) 
          nPar <- nPar + 1
        }
    }
  }
}

for (parConf in c(1:nParConfigs)){
  i <- paramConfigs[parConf,1]
  sdRecoveryTime <- paramConfigs[parConf, 2]
  avgRecoveryTime <- paramConfigs[parConf, 3]
  sDistFactor <- paramConfigs[parConf, 4]
  
    
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
  
  R0calc <- (3+2*nShort/N)*weighted.mean(x = sAgeDist/sDistFactor, w = c(3, 15, 25, 28, 22, 7))*avgRecoveryTime
  evalArray <- array(NA, dim = c(10000,nStatRun,7), dimnames = list(c(),c(),obsNames))
  evalR0 <- array(NA, dim = c(nStatRun,2), dimnames = list(c(), c("R0Mean", "R0MeanSd"))) 
  evalR02 <- array(NA, dim = c(nStatRun,2), dimnames = list(c(), c("R0Mean", "R0MeanSd")))
  evalR03 <- array(NA, dim = c(nStatRun,2), dimnames = list(c(), c("R0Mean", "R0MeanSd")))
  evalR04 <- array(NA, dim = c(nStatRun,2), dimnames = list(c(), c("R0Mean", "R0MeanSd")))
  evalR07 <- array(NA, dim = c(nStatRun,2), dimnames = list(c(), c("R0Mean", "R0MeanSd")))
  for(statRun in c(1:nStatRun)){
    simulationRun(statRun)
    runTime <- runTime + 1
    print(paste(c(runTime,timeGauge), sep="", collapse = " / "))
  }
  writeEval(evalArray, evalR0, params = params)
  write.table(x = evalR02, file = paste(c(path, "R02Mean_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
  write.table(x = evalR03, file = paste(c(path, "R03Mean_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
  write.table(x = evalR04, file = paste(c(path, "R04Mean_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
  write.table(x = evalR07, file = paste(c(path, "R07Mean_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
}  

#})

#-----------------------------------------------------------------------------------
#observed that for large (>100) lattices the boundaries are basically irrelevant for the largest cluster size

#population definieren (2-3) (different sDistFactor --> FLATTEN THE CURVE)
# verschieden viele nShort N=400**2 [N/100, N/50, N/10 ] [example: 0, boolsche susc]
  # berechne R0 erwartet | #!![3 + 1*p(nShort)] (3+2*nShort/N) Verbindungen pP | p(spread) = E(s*sDistFactor) | R0 = (3+2*n...)*p(spread)*trec = (3.2*weighted.mean(x = sAgeDist/sDistFactor, w = c(3, 15, 25, 28, 22, 7))*3.55)
  # R0 erwartet >>1 <<1 =1 [] 
# in file schreiben? 

# ausf√ºhren f√ºr verschiedene kofigs (susc / bond prob / recovery times) || param variation -> errors?? test
# verschiedene susc (verteilung/faktor), powerBond, recTime | jeweils mehrmals laufen lassen -> statistische auswertung
# auswerten! (neues programm) -define epidemic outbreak - errors!! ( stat analysis)
# !! laufzeit absch‰tzen!!

# -> aussuchen von 3-4 krankheiten mit verschiedenen besonderheiten (aussterben, starke verbreitung...)

#mit ausghesuchten krankehtien -> vaccination analyse
# vaccination offset
# immunity von 0 bis 0.9 rel kleinschrittig durchsim
# social Distancing can be changed during epidemic (time delay with seed to compare?)
# erst groﬂschrittig und dann bei  evtl. kippunkt  
#auswert: R != R0

#fﬂr andere pop testen --> unterschiede/gleich? --> statistische fehler der aussagen prﬂfen/abschﬂtzen


# recovery time vom alter abh√§ngig
# vaccination function nach sDist erstellung (nur f√ºr nicht alte / babies vaccination?)
# try to find a function  which describes R0 as a function of the parameters set
# bond prob dep on age dist
# nBonds <- numeric(N)
# for(i in c(1:N)){
#   nBonds[i] <- length(which(possConn == i))
# }
# print(mean(nBonds))

