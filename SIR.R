# https://www.healthknowledge.org.uk/public-health-textbook/research-methods/1a-epidemiology/epidemic-theory
library('plot.matrix')
source('modules.R')

#startTime <- proc.time()
#---------- Parameters to be set ----------------
profile <- profvis({
#set.seed(1)

immunity <- 0.3 #ratio of immune people 
bondOccProb <- 0.5

#network
N = 400**2 # number of people
nShort <- 100 # how many shortcuts are created // not impossibly many, could lead to ~inf loop #!!
openBoundaries <- FALSE # if FALSE the opposing edges of the lattice are connected (periodic)

sBool <- FALSE # if True the susceptibility is 1 or 0 // other poss like gaussian with age etc
sFixed <- FALSE
sNot <- TRUE

#observables
#plotting
plotIt <- TRUE
plotEvery <- 10

#clusters
checkCluster <- TRUE
clusterEvery <- 10

count <- 0
auswertungsVector <- array(0, dim=c(100,7), dimnames = list(c(),c("time","largestCluster", "numberCluster","numberInfected","largeOverTotal","largeOverRest","accInfections")))

# pFrom <- 0 # which canonical Q(p) are created
# pTo <- 1
# nProb <- 40 # how many datapoints are calculated in the above range
# nTest <- 10 # how many times is the same thing done


#-----------initialise network---------------------
lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N),3),dimnames = list(c(),c() , c("id", "I", "S")))

if(!openBoundaries){# these are the possible connections in 2D lattice
  possConn <- array(0,dim=c((2*N+nShort),3)) 
} else {
  possConn <- array(0,dim=c((2*N-2*sqrt(N)+nShort),3))
}

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
initialsDistribution <- sDistribution
lattice[,,3] <- sDistribution # set a uniform suceptibility 
#at the beginning nobody is infected
lattice[,,2] <- 0
infectionTime <- numeric(length=N)

#find all possible connections and their likelihoods
findConn()

nCon <- nrow(possConn)

#initialise observables
#binoms <- array(0, dim=c(nCon,nProb))
#percolTest <- array(0, dim = c(nCon,nTest))
#largestWeight <- array(0, dim = c(nCon,nTest))
#initializeBinoms()

# infect Patient 0
patientZero()

timesteps <- function(){
  #infected people
  infPeople <- which(infectionTime != 0)
  if (length(infPeople) == 0){
    stop(... = "No more infected", call. = TRUE)
  }
  #infected connections
  infConn <- possConn[which(possConn[,1] %in% infPeople | possConn[,2] %in% infPeople),c(1,2)]
  
  #infect people
  possiblyInf <- infConn[which(!(infConn %in% infPeople))]
  possiblyInf <- possiblyInf[which(!(possiblyInf %in% infPeople))]
  randVec <- runif(n = length(possiblyInf))
  susc <- sDistribution[possiblyInf]  
  newlyInf <- possiblyInf[susc>=randVec]
  
  #recovery // remove people
  recoveryCheck <- findRecBool(infectionTime[infPeople])
  recPeople <- infPeople[recoveryCheck]
  infPeople <- infPeople[!recoveryCheck]
  infectionTime[recPeople] <<- 0
  sDistribution[recPeople] <<- 0
  
  #increase infection time 
  infectionTime[infPeople] <<- infectionTime[infPeople] + 1
  infectionTime[newlyInf] <<- 1
}



findRecBool <- function(t){
  p <- 0.5-1/t
  return(p>=runif(length(t))) #!! more realistic: sample gaussian? 
}

x <- 0
while (TRUE) {

    
  timesteps()
  if((x %% plotEvery == 0) && (plotIt == TRUE)){
    visibleLattice <- array(0, dim= c(sqrt(N),sqrt(N)))
    visibleLattice[which(infectionTime != 0)] <- 1
    plot(which(visibleLattice==1, arr.ind = TRUE)[,1], which(visibleLattice==1, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)))
    #plot(visibleLattice)
  }

  if ((x %% clusterEvery == 0) && (checkCluster == TRUE)){
    count <- count +1
    people <- c(1:N)
    weight <- numeric(N) #rep(1,N)
    # now identifiy clusters as before, but instead of the connections sampled from possConn use infConn
    infPeople <- which(infectionTime != 0)
    weight[infPeople] <- 1
    infConn <- possConn[which((possConn[,1] %in% infPeople) & (possConn[,2] %in% infPeople)),c(1,2)]
    if (!is.array(infConn)){infConn <- array(data = infConn, dim = c(1,2))}
    if(nrow(infConn) != 0){
      for(i in c(1:nrow(infConn))){
        
        addConnection(infConn[i,1],infConn[i,2])
      }
    }
    auswertungsVector[count,1] <- x
    auswertungsVector[count,2] <- weight[which.max(weight)]
    auswertungsVector[count,3] <- length(which(weight!=0))
    auswertungsVector[count,4] <- sum(weight)
    auswertungsVector[count,5] <- weight[which.max(weight)]/sum(weight)
    auswertungsVector[count,6] <- weight[which.max(weight)]/sum(weight[-which.max(weight)]) #! is this right?
    auswertungsVector[count,7] <- length(which(sDistribution != initialsDistribution))/N
    # print(paste("largest cluster:" , weight[which.max(weight)], sep = " "))
    # print(paste("Number cluster:" , length(which(weight!=0)), sep = " "))
    # print(paste("Infected:" , sum(weight), sep = " "))
    # print(paste("Largest cluster over total infected:" , weight[which.max(weight)]/sum(weight), sep = " "))
    # print(paste("Largest cluster over smaller clusters:" , weight[which.max(weight)]/sum(weight[-which.max(weight)]), sep = " "))
    #!! interessantes in auswertungsvektor schreiben? 
     
  }
  
  x <- x + 1
}

# snapshot <- function(){
#   identify clusters
#   measure observables
# }

#!! cluster identification (rework)
# for (a in c(1:nTest)){
#   people <- c(1:N)
#   weight <- rep(1,N)
#   maxWeight <- 1
#   connections <- possConn[sample(nCon,nCon,replace=FALSE,prob = possConn[,3]),c(1,2)] #choose bonds which will be occupied gradually
#   for(x in c(1:nCon)){ 
#     if(checkLargestCluster){ # this if is really just for comfort, can remove it if it takes too long
#       newWeight <- addConnection(connections[x,2], connections[x,1])
#       if(newWeight > maxWeight){
#         maxWeight <- newWeight
#       }
#       largestWeight[x,a] <- maxWeight
#       #largestWeight[x,a] <- weight[which.max(weight)] # this is sligthly faster than max(weight)
#     }else {
#       addConnection(connections[x,2], connections[x,1])
#     }
#   }
#   
# }

#endTime <- proc.time()
#runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")


})

#-----------------------------------------------------------------------------------
#observed that for large (>100) lattices the boundaries are basically irrelevant for the largest cluster size
# bond occupation probablitiy include? in infection / in possCOnn.. ?? at all?? ->  transmissibility? 
# transmissiblity: the same as S?