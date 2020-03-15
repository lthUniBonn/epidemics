# https://www.healthknowledge.org.uk/public-health-textbook/research-methods/1a-epidemiology/epidemic-theory
library('plot.matrix')
source('modules.R')

#startTime <- proc.time()
#---------- Parameters to be set ----------------
#profile <- profvis({

#set.seed(1)

immunity <- 0 #ratio of immune people 
bondOccProb <- 0.9

params <- paste(c(immunity, bondOccProb, avgRecoveryTime), sep="", collapse="_")
#!! match recovery with age?? 
avgRecoveryTime <- 3.55
sdRecoveryTime <- 1

sAgeDist <- c(1, 0.7, 0.5, 0.7, 0.8, 1)/3.5 #susceptibility depending on age
#!! bond prob dep on age dist
#network
N = 400**2 # number of people
nShort <- 0 # how many shortcuts are created // not impossibly many, could lead to ~inf loop #!!
periodicBoundaries <- TRUE #if FALSE the opposing edges of the lattice are connected (periodic)

sBool <- FALSE # if True the susceptibility is 1 or 0 // other poss like gaussian with age etc
sFixed <- FALSE
sNot <- FALSE
sReal <- TRUE
#observables
   #plotting
plotIt <- FALSE
plotEvery <- 100
plotAccumulated <- TRUE
   #clusters
checkCluster <- TRUE
clusterEvery <- 10

count <- 0
auswertungsVector <- array(0, dim=c(1000,8), dimnames = list(c(),c("time","largestCluster", "numberCluster","numberInfected","largeOverTotal","largeOverRest","accInfections", "R0")))
R0OverInfectiousPeriod <- numeric(length=N)
R0Mean <- numeric(10000)
# pFrom <- 0 # which canonical Q(p) are created
# pTo <- 1
# nProb <- 40 # how many datapoints are calculated in the above range
# nTest <- 10 # how many times is the same thing done


#-----------initialise network---------------------
lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N)))

if(periodicBoundaries){# these are the possible connections in 2D lattice
  possConn <- array(0,dim=c((2*N+nShort),3)) 
} else {
  possConn <- array(0,dim=c((2*N-2*sqrt(N)+nShort),3))
}
#set Infection status
infected <- logical(length = N)
infectionTime <- numeric(length=N)
recovered <- logical(length=N)
#set susceptibility distribution in population
if(sBool){
  sDistribution <- sample(c(0,1), replace=TRUE, size=N,prob = c(immunity,1-immunity))
}
if(sFixed){
  sDistribution <- sample(c(0,0.25,0.5,0.75,1), replace=TRUE, size=N,prob = c(immunity,(1-immunity)/8,  (1-immunity)/2, (1-immunity)/4, (1-immunity)/8))
}
if(sNot){
  sDistribution <- sample(c(0.25,0.5,0.75,1), replace=TRUE, size=N,prob = c((1-immunity)/8,  (1-immunity)/2, (1-immunity)/4, (1-immunity)/8))
#!! implement a scaling for the susc  // age
}
if(sReal){
  ageDistribution <- sample(c(3, 20,40,60,80, 100), replace=TRUE, size=N,prob = c(3, 15, 25, 28, 22, 7))
  sDistribution <- unlist(lapply(ageDistribution,FUN = sAge))
}






#save Distribution of initially susceptible people for eval
initialsDistribution <- sDistribution

#find all possible connections and their bond probability
findConn()
nCon <- nrow(possConn)

# infect Patient 0
patientZero()


x <- 0
while (TRUE) {
  R0 <- timesteps()
  # returns NA if no more infected
  if(is.na(R0)){ break}
  # accumulated plot of infected | currently infected in red
  if((x %% plotEvery == 0) && (plotAccumulated == TRUE)){
    visibleLattice <- array(0, dim= c(sqrt(N),sqrt(N)))
    visibleLattice[which(sDistribution  != initialsDistribution)] <- 1
    plot(which(visibleLattice==1, arr.ind = TRUE)[,1], which(visibleLattice==1, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '.', pty = "s")
    par(new=TRUE)
    visibleLattice[which(infected == TRUE)] <- 2
    plot(which(visibleLattice==2, arr.ind = TRUE)[,1], which(visibleLattice==2, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '*',col = 'red', pty = "s")
  }
  # evaluation of disease parameters
  if ((x %% clusterEvery == 0) && (checkCluster == TRUE)){
    count <- count +1
    people <- c(1:N)
    weight <- numeric(N) 
    # now identifiy clusters use Newman Algo with infConn as used bonds
    infPeople <- which(infectionTime != 0)
    weight[infPeople] <- 1
    infConn <- possConn[which((possConn[,1] %in% infPeople) & (possConn[,2] %in% infPeople)),c(1,2)]
    if (!is.array(infConn)){infConn <- array(data = infConn, dim = c(1,2))}
    if(nrow(infConn) != 0){
      for(i in c(1:nrow(infConn))){
        
        addConnection(infConn[i,1],infConn[i,2])
      }
    }
    if (x == 0){#write headers + overwrite old files to start appending to fresh file
      write(x = c("#x","maxWeight"), file = paste(c("data/maxWeight_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
      write(x = c("#x","numberCluster"), file = paste(c("data/numberCluster_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
      write(x = c("#x","numberInfected"), file = paste(c("data/numberInfected_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
      write(x = c("#x","largeOverTotal"), file = paste(c("data/largeOverTotal_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
      write(x = c("#x","accInfections"), file = paste(c("data/accInfections_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
      write(x = c("#x","R0Mean"), file = paste(c("data/R0Mean_", params, ".txt"),sep="", collapse=""), append = FALSE, sep = "\t", ncolumns = 2)
    }
    write(x = c(x, weight[which.max(weight)]), file = paste(c("data/maxWeight_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")
    write(x = c(x, length(which(weight!=0))), file = paste(c("data/numberCluster_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")
    write(x = c(x, sum(weight)), file = paste(c("data/numberInfected_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")
    write(x = c(x, weight[which.max(weight)]/sum(weight)), file = paste(c("data/largeOverTotal_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")
    write(x = c(x, (length(which(sDistribution != initialsDistribution))+length(infPeople))/N), file = paste(c("data/accInfections_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")
    if (x != 0){write(x = c(x, R0Mean[x]), file = paste(c("data/R0Mean_", params, ".txt"),sep="", collapse=""), append = TRUE,sep = "\t")}
    # auswertungsVector[count,1] <- x
    # auswertungsVector[count,2] <- weight[which.max(weight)]
    # auswertungsVector[count,3] <- length(which(weight!=0))
    # auswertungsVector[count,4] <- sum(weight)
    # auswertungsVector[count,5] <- weight[which.max(weight)]/sum(weight)
    # auswertungsVector[count,6] <- weight[which.max(weight)]/sum(weight[-which.max(weight)]) #! is this right?
    # auswertungsVector[count,7] <- (length(which(sDistribution != initialsDistribution))+length(infPeople))/N
    # auswertungsVector[count,8] <- R0
  }
  
  x <- x + 1
}

#read data from file 
maxWeightDf <- read.table(file = paste(c("data/maxWeight_", params, ".txt"),sep="", collapse=""))
numberClusterDf <- read.table(file = paste(c("data/numberCluster_", params, ".txt"),sep="", collapse=""))
numberInfectedDf <- read.table(file = paste(c("data/numberInfected_", params, ".txt"),sep="", collapse=""))
largeOverTotalDf <- read.table(file = paste(c("data/largeOverTotal_", params, ".txt"),sep="", collapse=""))
accInfectionsDf <- read.table(file = paste(c("data/accInfections_", params, ".txt"),sep="", collapse=""))
R0MeanDf <- read.table(file = paste(c("data/R0Mean_", params, ".txt"),sep="", collapse=""))

plot(maxWeightDf[,1], maxWeightDf[,2], xlab = "T",ylab = "largest Cluster")
plot(numberClusterDf[,1],numberClusterDf[,2], xlab = "T",ylab = "Number of Clusters")
plot(numberInfectedDf[,1],numberInfectedDf[,2], xlab = "T",ylab = "Number Infected")
plot(largeOverTotalDf[,1],largeOverTotalDf[,2], xlab = "T",ylab = "Largest Cluster over Total Infected")
plot(accInfectionsDf[,1],accInfectionsDf[,2], xlab = "T",ylab = "Attack Rate")
plot(R0MeanDf[,1],R0MeanDf[,2], xlab = "T",ylab = "R0 Mean")

# auswertungsVector <- auswertungsVector[c(1:count),]
# plot(auswertungsVector[,1],auswertungsVector[,2], xlab = "T",ylab = "largest Cluster")
# plot(auswertungsVector[,1],auswertungsVector[,3], xlab = "T",ylab = "Number of Clusters")
# plot(auswertungsVector[,1],auswertungsVector[,4], xlab = "T",ylab = "Number Infected")
# plot(auswertungsVector[,1],auswertungsVector[,5], xlab = "T",ylab = "Largest Cluster over Total Infected")
# plot(auswertungsVector[,1],auswertungsVector[,6], xlab = "T",ylab = "Largest Cluster over Smaller Clusters")
# plot(auswertungsVector[,1],auswertungsVector[,7], xlab = "T",ylab = "Attack Rate")
# plot(auswertungsVector[,1],auswertungsVector[,8], xlab = "T",ylab = "R0")

#endTime <- proc.time()
#runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")


#})

#-----------------------------------------------------------------------------------
#observed that for large (>100) lattices the boundaries are basically irrelevant for the largest cluster size



# ausführen für verschiedene susc / bond prob / recovery times 
# recovery time vom alter abhängig
# vaccination function nach sDist erstellung (nur für nicht alte / babies vaccination?)



