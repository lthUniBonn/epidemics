# create 2d lattice 
# add shortcuts (density of these is important says newman)
# SIR model - suscetible, infectious, removed
# S: prob that an exposed individual contract
# transmissiblity: the same as S?

#binomial koeffizienten in file schreiben 

#right now: the largest cluster size is the maximum spread of the disease 
# this will happen if patient zero is within this cluster
# if patient zero is in a smaller cluster, that means (s)he does not have contact to as many people and the disease dies
# largest cluster gives a sort of worst case

# visualization? 
library('plot.matrix')
library('gmp')
library('profvis')
library('Brobdingnag')

#startTime <- proc.time()
#---------- Parameters to be set ----------------
#profile <- profvis({
N = 100**2 # number of people
immunity <- 0.4
nShort <- 0 # how many shortcuts are created
pFrom <- 0 # which canonical Q(p) are created
pTo <- 1
nProb <- 40 # how many datapoints are calculated in the above range
nTest <- 10 # how many times is the same thing done
checkPercolation <- FALSE
checkLargestCluster <- TRUE
openBoundaries <- FALSE # if FALSE the opposing edges of the lattice are connected
sBool <- TRUE # if True the susceptibility is 1 or 0
#observed that for large (>100) lattices the boundaries are basically irrelevant for the largest cluster size

#-----------global declarations---------------------
lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N),3),dimnames = list(c(),c() , c("id", "I", "S")))

#set.seed(1)


merges <- 0 # how many connections were already made
if(!openBoundaries){
  possConn <- array(0,dim=c((2*N+nShort),3)) # these are the possible connections
} else {
  possConn <- array(0,dim=c((2*N-2*sqrt(N)+nShort),3))
}
topIndices <- seq(1, N-sqrt(N)+1, by = sqrt(N))
botIndices <- seq(sqrt(N), N, by = sqrt(N))
leftIndices <- c(1:sqrt(N))
rightIndices <- c((N-sqrt(N)+1):N)


nCon <- nrow(possConn)
percolTest <- array(0, dim = c(nCon,nTest))
largestWeight <- array(0, dim = c(nCon,nTest))

conProb <- numeric(nCon)
#--------------functions-------------------------
# find possible connections in 2d case
findConn <- function(){
  counter <- 0
  for(j in c(1:(sqrt(N)-1))){
    for(i in c(1:(sqrt(N)-1))){
      #poss con add C(lattice[i,j], lattice [i+1,j]) to the bottom and to the right
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,j,1]  
      possConn[counter,2] <<- lattice [i+1,j,1]
      possConn[counter,3] <<- transProb(lattice[i,j,3],lattice [i+1,j,3])
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j,1] 
      possConn[counter,2] <<- lattice [i,j+1,1]
      possConn[counter,3] <<- transProb(lattice[i,j,3],lattice [i,j+1,3])
    }
    #poss con add C(lattice[i,j], lattice [i+1,j]) to the bottom and to the right
    counter <- counter + 1
    possConn[counter,1] <<- lattice[sqrt(N),j,1]  
    possConn[counter,2] <<- lattice [sqrt(N),j+1,1] 
    possConn[counter,3] <<- transProb(lattice[sqrt(N),j,3],lattice [sqrt(N),j+1,3])
  }
  for (i in c(1:(sqrt(N)-1))){
    counter <- counter + 1
    possConn[counter,1] <<- lattice[i,sqrt(N),1]  
    possConn[counter,2] <<- lattice [i+1,sqrt(N),1]
    possConn[counter,3] <<- transProb(lattice[i,sqrt(N),3],lattice[i+1,sqrt(N),3])
  }
  
  if(!openBoundaries){
    for(i in c(1:length(topIndices))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[1,i,1]
      possConn[counter,2] <<- lattice[sqrt(N),i,1]
      possConn[counter,3] <<- transProb(lattice[1,i,3],lattice[sqrt(N),i,3])
    }
    for(i in c(1:length(leftIndices))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,1,1]
      possConn[counter,2] <<- lattice[i,sqrt(N),1]
      possConn[counter,3] <<- transProb(lattice[i,1,3],lattice[i,sqrt(N),3])
    }
  }
  
  if(nShort != 0){
    counter <- counter + 1
    fromRow <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
    fromCol <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
    toRow   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
    toCol   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
    possConn[c(counter:(counter+nShort-1)),1] <<- lattice[fromRow,fromCol,1]
    possConn[c(counter:(counter+nShort-1)),2] <<- lattice[toRow,toCol,1]
    possConn[c(counter:(counter+nShort-1)),3] <<- transProb(lattice[fromRow,fromCol,3],lattice[toRow,toCol,3])
    # get duplicates and self connections
  }
}

transProb <- function(x,y){ #this function determines the likelihood of transmission along the connection from x to y
  return(x*y) #simplest function for now
}


findRoot <- function(startIndex){
  #find root plus path compression
  root <- startIndex
  while (people[root] != root){
    root <- people[root]
  }
  while (people[startIndex] != root){
    parent <- people[startIndex]
    people[startIndex] <<- root # path compression:  -> always point to root
    startIndex <- parent
  }
  return(root)
}

addConnection <- function(from, to){ # does this work with shortcuts? I think yes, but not 100% sure
  #definetly get a problem when people can recover, as this should interrupt the tree
  #this cannot be taken care of in hindsight, as the trees do not contain this information anymore
  # maybe sacrifice the path compression?
  #find roots for both parts of conn (A,B) [path compression]
  fromRoot <- findRoot(from)#index of root node 
  toRoot <- findRoot(to)
  if (fromRoot == toRoot){
    return(0)
  }# if same: nothing
  else {
    if (weight[fromRoot] > weight[toRoot]){
      tmp <- fromRoot
      fromRoot <- toRoot
      toRoot <- tmp
    }
    # write connection "from root node" to "to root node"
    people[fromRoot] <<- toRoot 
    weight[toRoot] <<- weight[fromRoot] + weight[toRoot] # this line takes longer (50-200x) when checkLargestCluster = True. wat?
    #.subset(weight)[toRoot] <<- .subset(weight)[fromRoot] + .subset(weight)[fromRoot]
    weight[fromRoot] <<- 0
    merges <<- merges +1
    
  }
  return(weight[toRoot])
}


isPercolating <- function(){ # this does not work for lattices with boundary conditions
  leftRoots <- sapply(leftIndices, findRoot)
  rightRoots <- sapply(rightIndices, findRoot)
  topRoots <- sapply(topIndices, findRoot)
  botRoots <- sapply(botIndices, findRoot)
  percolatingLeftRight <- numeric()
  percolatingTopBot <- numeric()
  percolatingLeftRight <-Reduce(intersect, list(leftRoots,rightRoots))
  percolatingTopBot <-Reduce(intersect, list(topRoots,botRoots))
  if(length(percolatingLeftRight) != 0){
    return(TRUE)
  }
  if(length(percolatingTopBot) != 0){
    return(TRUE)
  }
  return(FALSE)
}

erf <- function(x,a,b) (pnorm(a*(x-b) * sqrt(2)))

canonical <- function(micObs){#find p from n
  # does this still work when shortcuts are introduced? question is wether it makes a difference that upper bound of n rises but p stays the same
  # I think this is okay, but need to talk about it
  canObs <- numeric(nProb)
  i <- 0
  for (p in c(1:nProb)){ 
    i <- i+1
    canObs[i] <- sum(binoms[,i]*micObs)
  }
  #percolProb[nProb] <- percolTest[nCon] #per def
  return(canObs)
}

erf <- function(x,a,b) (pnorm(a*(x-b) * sqrt(2)))

initializeBinoms <- function(){
  i <- 0
  print(nCon)
  for (p in seq(pFrom, pTo, (pTo-pFrom)/(nProb-1))){
    i <- i+1
    binoms[,i] <<- dbinom(x = c(1:nCon),size = nCon, prob = p)
  }
}


#-----------------------main--------------------------------------

if(sBool){
  sDistribution <- sample(c(0,1), replace=TRUE, size=N,prob = c(immunity,1-immunity))
}
lattice[,,3] <- sDistribution # set a uniform suceptibility 
findConn()#find all possible connections and their likelihoods
#remove the impossible connections
possConn <- possConn[-which(possConn[,3]==0),]
nCon <- nrow(possConn)
binoms <- array(0, dim=c(nCon,nProb))
initializeBinoms()
for (a in c(1:nTest)){
  people <- c(1:N)
  weight <- rep(1,N)
  maxWeight <- 1
  connections <- possConn[sample(nCon,nCon,replace=FALSE,prob = possConn[,3]),c(1,2)] #choose bonds which will be occupied gradually
  for(x in c(1:nCon)){ 
     if(checkLargestCluster){ # this if is really just for comfort, can remove it if it takes too long
       newWeight <- addConnection(connections[x,2], connections[x,1])
       if(newWeight > maxWeight){
         maxWeight <- newWeight
       }
       largestWeight[x,a] <- maxWeight
       #largestWeight[x,a] <- weight[which.max(weight)] # this is sligthly faster than max(weight)
     }else {
       addConnection(connections[x,2], connections[x,1])
     }
    
    #does not make sense to check for percolation if shprtcuts are introduced, od does it?
    #certainly not when closed boundary conditions are present
    # if boundary conditions are turned on this needs changing as well
    # if(isPercolating()){
    #   percolTest[c(x:nCon),a] <- percolTest[c(x:nCon),a] + 1
    #   break
    # }
    # if(checkLargestCluster){ # this if is really just for comfort, can remove it if it takes too long
    #   largestWeight[x,a] <- weight[which.max(weight)] # this is sligthly faster than max(weight)
    # }
  }
  
}

#----------- display data -----------------


# percolation threshold finden als checkup




if(checkLargestCluster){
  xData <- seq(pFrom, pTo, (pTo-pFrom)/(nProb-1))
  yData <- array(0,c(nProb,nTest))
  ySEM <- numeric(nProb)
  for(i in c(1:nTest)){
    yData[,i] <- canonical(largestWeight[,i]) 
  }
  for(i in c(1:nProb)){
    ySEM[i] <- sd(yData[i,])/(sqrt(nTest)*N)
  }
  yMean <- rowMeans(yData)/N
  yData <- yData/N
  largestCluster <- data.frame( x= xData, y=yMean)
  plot(largestCluster,xlab = "", ylab = "",xlim = c(0,1), ylim = c(0,1))
  arrows(xData, yMean-ySEM, xData, yMean+ySEM, length=0.05, angle=90, code=3)
  title(main="Largest Cluster Size",  xlab="Q(p)", ylab="Size") 
  text(x = pFrom, y = 1, labels = paste("N:", N,"  nShort:", nShort, "  nTest:", nTest,sep = " "),pos = 4)
  abline(h=1-immunity)
  
  #ourFit <- nls(y ~ erf(x,a,b), data = largestCluster, start=list(a=50, b=0.4))
  # this ignores the known errors
  #x <- seq(0.2,1,by=0.001)
  #lines(y=erf(x,summary(ourFit)$coefficients[1], summary(ourFit)$coefficients[2]), x = x)
}

#percolProb <- data.frame( x= seq(pFrom, pTo, (pTo-pFrom)/(nProb-1)), y=canonical(percolTest))
#plot(percolProb)

#endTime <- proc.time()
#runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")


#})
