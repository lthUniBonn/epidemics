# create 2d lattice 
# add shortcuts (density of these is important says newman)
# SIR model - suscetible, infectious, removed
# S: prob that an exposed individual contract
# transmissiblity: the same as S?

#binomial koeffizienten in file schreiben 


# visualization? 
library('plot.matrix')
library('gmp')
library('profvis')
library('Brobdingnag')
#startTime <- proc.time()
#---------- Parameters to be set ----------------
profile <- profvis({
N = 200**2 # number of people
nShort <- 0 # how many shortcuts are created
pFrom <- 0.2 # which canociacl Q(p) are created
pTo <- 0.8
nProb <- 100 # how many datapoints are calculated in the above range
nTest <- 1 # how many times is the same thing done
# we completely ignored this so far....
checkPercolation <- FALSE 
checkLargestCluster <- TRUE

#-----------global declarations---------------------
lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N),3),dimnames = list(c(),c() , c("id", "I", "S")))

#set.seed(1)

merges <- 0 # how many connections were already made
dof <- N*(N+1)/2 # degrees of freedom in symmetric matrix
possConn <- array(0,dim=c((2*N+nShort),2)) # these are the possible connections

topIndices <- seq(1, N-sqrt(N)+1, by = sqrt(N))
botIndices <- seq(sqrt(N), N, by = sqrt(N))
leftIndices <- c(1:sqrt(N))
rightIndices <- c((N-sqrt(N)+1):N)
#now the indices that need to be checked are all in useful indices, this reduces the array size by N^2-N*(N+1)/2, for 100*100 grid is 4950 sites smaller

#connections is the list of random connections that will be made

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
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j,1] 
      possConn[counter,2] <<- lattice [i,j+1,1]
      
    }
    #poss con add C(lattice[i,j], lattice [i+1,j]) to the bottom and to the right
    counter <- counter + 1
    possConn[counter,1] <<- lattice[sqrt(N),j,1]  
    possConn[counter,2] <<- lattice [sqrt(N),j+1,1] 
  }
  for (i in c(1:(sqrt(N)-1))){
    counter <- counter + 1
    possConn[counter,1] <<- lattice[i,sqrt(N),1]  
    possConn[counter,2] <<- lattice [i+1,sqrt(N),1]  
  }
  for(i in c(1:length(topIndices))){
    counter <- counter + 1
    possConn[counter,1] <<- topIndices[i]
    possConn[counter,2] <<- botIndices[i]
  }
  for(i in c(1:length(leftIndices))){
    counter <- counter + 1
    possConn[counter,1] <<- leftIndices[i]
    possConn[counter,2] <<- rightIndices[i]
  }
  if(nShort != 0){
    counter <- counter + 1
    possConn[c(counter:(counter+nShort-1)),1] <<- sample(c(1:N),size = nShort,replace = TRUE)
    possConn[c(counter:(counter+nShort-1)),2] <<- sample(c(1:N),size=nShort,replace = TRUE)
    # get duplicates and self connections
  }
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
    weight[toRoot] <<- weight[fromRoot] + weight[toRoot]
    weight[fromRoot] <<- 0
    merges <<- merges +1
  }
  
}


isPercolating <- function(){ # this does not work for lattices with boundary conditions
  leftRoots <- sapply(leftIndices, findRoot)
  rightRoots <- sapply(rightIndices, findRoot)
  #allIndices <- c(1:N)
  #topIndices <- allIndices[which((allIndices %% sqrt(N))== 1)]
  #botIndices <- allIndices[which((allIndices %% sqrt(N))== 0)]
  
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
  # again: I think this is okay, but need to talk about it
  canObs <- numeric(nProb)
  i <- 0
  for (p in seq(pFrom, pTo, (pTo-pFrom)/(nProb-1))){ 
    # think about next lines, it meight be fine after all
    # we need to change something about how p is calculated
    # if we want to compare different shortcuts
    # For example: 100*100 Lattice has 20000 possible connections, 
    # if we add 1000 shortcuts p should not be 100% if 20000 connections are filled but when 21000 connections were made
    i <- i+1
    binoms <- double(length=nCon)
    binoms <- dbinom(x = c(1:nCon),size = nCon, prob = p)*micObs
    canObs[i] <- sum(binoms)
  }
  #p <- 1 
  
  #percolProb[nProb] <- percolTest[nCon] #per def
  return(canObs)
  
}

erf <- function(x,a,b) (pnorm(a*(x-b) * sqrt(2)))

#main

findConn()#find all possible connections
nCon <- nrow(possConn)
percolTest <- array(0, dim = c(nCon,nTest))
largestWeight <- array(0, dim = c(nCon,nTest))
for (a in c(1:nTest)){
  people <- c(1:N)
  weight <- rep(1,N)
  connections <- possConn[sample(nCon,nCon,replace=FALSE),] #choose bonds which will be occupied gradually
  for(x in c(1:nCon)){ 
    addConnection(connections[x,2], connections[x,1])
    #does not make sense to check for percolation if shprtcuts are introduced, od does it?
    #certainly not when closed boundary conditions are present
    # if boundary conditions are turned on this needs changing as well
    # if(isPercolating()){
    #   percolTest[c(x:nCon),a] <- percolTest[c(x:nCon),a] + 1
    #   break
    # }
    
    largestWeight[x,a] <- max(weight)
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
  plot(largestCluster,xlab = "", ylab = "")
  arrows(xData, yMean-ySEM, xData, yMean+ySEM, length=0.05, angle=90, code=3)
  title(main="Largest Cluster Size",  xlab="Q(p)", ylab="Size") 
  text(x = 0.2, y = 1, labels = paste("N:", N,"  nShort:", nShort, "  nTest:", nTest,sep = " "),pos = 4)
  #mtext(paste("N:", N,"  nShort:", nShort, "  nTest:", nTest,sep = " "),side = 3)
  
  ourFit <- nls(y ~ erf(x,a,b), data = largestCluster, start=list(a=50, b=0.4))
  # this ignores the known errors
  x <- seq(0.2,1,by=0.001)
  lines(y=erf(x,summary(ourFit)$coefficients[1], summary(ourFit)$coefficients[2]), x = x)
}

#percolProb <- data.frame( x= seq(pFrom, pTo, (pTo-pFrom)/(nProb-1)), y=canonical(percolTest))
#plot(percolProb)

#endTime <- proc.time()
#runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")


})
