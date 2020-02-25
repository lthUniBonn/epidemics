# create 2d lattice 
# add shortcuts (density of these is important says newman)
# SIR model - suscetible, infectious, removed
# S: prob that an exposed individual contract
# transmissiblity: the same as S?

#binomial koeffizienten in file schreiben 


# visualization? 
library('gmp')
library('profvis')
library('Brobdingnag')
#startTime <- proc.time()
profile <- profvis({
N = 100**2 # number of people

lattice <- array(data=c(1:N), dim = c(sqrt(N),sqrt(N),3),dimnames = list(c(),c() , c("id", "I", "S")))
nShort <- sqrt(N)
#set.seed(1)

merges <- 0 # how many connections were already made
dof <- N*(N+1)/2 # degrees of freedom in symmetric matrix
possConn <- array(0,dim=c((2*N-2*sqrt(N)+nShort),2)) # these are the possible connections

topIndices <- seq(1, N-sqrt(N)+1, by = sqrt(N))
botIndices <- seq(sqrt(N), N, by = sqrt(N))
leftIndices <- c(1:sqrt(N))
rightIndices <- c((N-sqrt(N)+1):N)
#now the indices that need to be checked are all in useful indices, this reduces the array size by N^2-N*(N+1)/2, for 100*100 grid is 4950 sites smaller

#connections is the list of random connections that will be made

# find possible connections in 2d case
findConn <- function(){
  counter <- 0
  for(j in c(1:(sqrt(N)-1))){
    for(i in c(1:(sqrt(N)-1))){
      #poss con add C(lattice[i,j], lattice [i+1,j]) to the bottom and to the right
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,j]  
      possConn[counter,2] <<- lattice [i+1,j] 
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j] 
      possConn[counter,2] <<- lattice [i,j+1]
      
    }
    #poss con add C(lattice[i,j], lattice [i+1,j]) to the bottom and to the right
    counter <- counter + 1
    possConn[counter,1] <<- lattice[sqrt(N),j]  
    possConn[counter,2] <<- lattice [sqrt(N),j+1] 
  }
  for (i in c(1:(sqrt(N)-1))){
    counter <- counter + 1
    possConn[counter,1] <<- lattice[i,sqrt(N)]  
    possConn[counter,2] <<- lattice [i+1,sqrt(N)]  
  }

  counter <- counter + 1
  possConn[c(counter:(counter+nShort-1)),1] <<- sample(c(1:N),size = 10,replace = TRUE)
  possConn[c(counter:(counter+nShort-1)),2] <<- sample(c(1:N),size=10,replace = TRUE)
  # get duplicates and self connections
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

addConnection <- function(from, to){
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


isPercolating <- function(){
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


#main
nTest <- 1
findConn()#find all possible connections
nCon <- nrow(possConn)
percolTest <- numeric(length=nCon)
for (a in c(1:nTest)){
  people <- c(1:N)
  weight <- rep(1,N)
  connections <- possConn[sample(nCon,nCon,replace=FALSE),] #choose bonds which will be occupied gradually
  print(a)
  for(x in c(1:nCon)){ 
    addConnection(connections[x,2], connections[x,1])
    #if(isPercolating()){
     # percolTest[c(x:nCon)] <- percolTest[c(x:nCon)] + 1
      #break
    #}
  }
}

# percolation threshold finden als checkup



#maybe calc once and write to file? takes looong
#nOverK <- as.brob(chooseZ(nCon, 1:nCon))

canonical <- function(){#find p from n
  for (i in seq(1,(nProb-1))){ #all p but p=1 makes problem in logarithmic expression
    print(i)
    p <- i / nProb
    percolBinom <- double(length=nCon)
    for (x in c(1:nCon)){#should start at 0 if observable is not 0 for n=0 
      percolBinom[x] <- as.numeric(nOverK[x]*(as.brob(p)**x)*(as.brob(1-p)**(nCon-x))*percolTest[x])
      if (x %% 1000 == 0){  print(x)}
    }
    percolProb[i] <- sum(percolBinom)
  }
#p <- 1 
percolProb[nProb] <- percolTest[nCon] #per def
}



#endTime <- proc.time()
#runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")


})
