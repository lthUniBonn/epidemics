library('plot.matrix')
library('gmp')
library('profvis')
library('Brobdingnag')


patientZero <- function(){
  zeroId <- sample(x = which(sDistribution == 1),size = 1) #choose one of the susceptibles to be patient 0
  infectionTime[zeroId] <<- 1 
  #latticeID <- which(lattice[,,1] == zeroId, arr.ind = TRUE)
  #lattice[latticeID[1], latticeID[2],2] <<- 1 
}

# find possible connections in 2d case next neighbor
findConn <- function(){
  counter <- 0
  for(j in c(1:(sqrt(N)-1))){
    for(i in c(1:(sqrt(N)-1))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,j,1]  
      possConn[counter,2] <<- lattice [i+1,j,1]
      possConn[counter,3] <<- transProb(lattice[i,j,3],lattice [i+1,j,3])
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j,1] 
      possConn[counter,2] <<- lattice [i,j+1,1]
      possConn[counter,3] <<- transProb(lattice[i,j,3],lattice [i,j+1,3])
    }
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
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[1,i,1]
      possConn[counter,2] <<- lattice[sqrt(N),i,1]
      possConn[counter,3] <<- transProb(lattice[1,i,3],lattice[sqrt(N),i,3])
    }
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,1,1]
      possConn[counter,2] <<- lattice[i,sqrt(N),1]
      possConn[counter,3] <<- transProb(lattice[i,1,3],lattice[i,sqrt(N),3])
    }
  }
  while (TRUE) {
    if(nShort != 0){
      counter <- counter + 1
      condition <- TRUE
      while (condition == TRUE) {
        fromRow <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
        fromCol <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
        toRow   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
        toCol   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
        for (i in c(1:nShort)){
          if (!(fromRow[i] == toRow[i] && fromCol[i] == toCol[i])){
            condition <- FALSE
            break
          }
        }
      }
      for (i in c(0:(nShort-1))){
      possConn[counter + i,1] <<- lattice[fromRow[i+1],fromCol[i+1],1]
      possConn[counter + i,2] <<- lattice[toRow[i+1],toCol[i+1],1]
      possConn[counter + i,3] <<- transProb(lattice[fromRow[i+1],fromCol[i+1],3],lattice[toRow[i+1],toCol[i+1],3])
    # not left < right with shortcuts #!!
      }
    counter <- counter - 1
    if (!(TRUE %in% duplicated(possConn))){break}
    }
  }
  # possConn <<- possConn[-which(possConn[,3]==0),] only works if which retruns something! only needed if 0 possible
}

#!! change later // get bond occupation probability  
transProb <- function(x,y){ #this function determines the likelihood of transmission along the connection from x to y
  return(bondOccProb) 
}


    #cluster evaluation
#------
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
    #merges <<- merges +1
    
  }
  return(weight[toRoot])
}

#------

erf <- function(x,a,b) (pnorm(a*(x-b) * sqrt(2)))

# High efficient Newman from microcanonical to canonical (binomial interpret) 
canonical <- function(micObs){#find p from n
  canObs <- numeric(nProb)
  i <- 0
  for (p in c(1:nProb)){ 
    i <- i+1
    canObs[i] <- sum(binoms[,i]*micObs)
  }
  return(canObs)
}

initializeBinoms <- function(){
  i <- 0
  print(nCon)
  for (p in seq(pFrom, pTo, (pTo-pFrom)/(nProb-1))){
    i <- i+1
    binoms[,i] <<- dbinom(x = c(1:nCon),size = nCon, prob = p) 
  }
}


#----------- display data -----------------

plotLargestCLuster <- function(){
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
  text(x = pFrom, y = 1, labels = paste("N:", N,"  nShort:", nShort, "  nTest:", nTest, "do not trust errorbars yet",sep = " "),pos = 4)
  abline(h=1-immunity)
}

