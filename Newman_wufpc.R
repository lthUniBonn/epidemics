#fast checking of trees (?)
# slow finding possible connections : findCon
# visualization? 

library('profvis')
startTime <- proc.time()
profile <- profvis({
N = 10**2 # number of people
No = sqrt(N)*2 #the No. of connections that will be added
people <- c(1:N)
set.seed(1)
weight <- rep(1,N)
visual <- rep(0,N)
merges <- 0 # how many connections were already made
dof <- N*(N+1)/2 # degrees of freedom in symmetric matrix
possConn <- array(0,dim=c(dof,2)) # these are the possible connections

#now the indices that need to be checked are all in useful indices, this reduces the array size by N^2-N*(N+1)/2, for 100*100 grid is 4950 sites smaller

#connections is the list of random connections that will be made

# find possible connections in 2d case
calcDistance<- function(first, second){ # calculate distance between people in 2dim
  rowfirst <- first %% sqrt(N)
  columnfirst <- first%/% sqrt(N)
  rowsecond <- second %% sqrt(N)
  columnsecond <- second %/% sqrt(N)
  distance <- sqrt((rowfirst-rowsecond)**2+(columnfirst - columnsecond)**2)
  return(distance)
}

findConn <- function(){
  counter <- 0
  for(j in c(1:N)){
    for(i in c(j:N)){
      counter <- counter + 1
      possConn[counter,1] <<- j # filters out the lower triangle 
      possConn[counter,2] <<- i 
      
    }
  }
  
  for(s in c(1:nrow(possConn))){
    if(calcDistance(possConn[s,1],possConn[s,2]) != 1){
      possConn[s,1] <<- 0 #filter out non-next-neighbour i.a.
      possConn[s,2] <<- 0
    }
  }
  possConn <<- possConn[-which(possConn ==0),] #remove filtered from array
}

findRoot <- function(startIndex){
  #print(startIndex)
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
  visual[from] <<- 1
  visual[to] <<- 1
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


#main
findConn()                                                        #find all possible connections
connections <- possConn[sample(nrow(possConn),No,replace=FALSE),] #choose bonds which will be occupied gradually
#for(x in c(1:No)){
  #addConnection(connections[x,2], connections[x,1])
  #evaluate!
#}
addConnection(connections[1,2], connections[1,1])
endTime <- proc.time()
runTime <- endTime-startTime
#write(x = c(N, runTime[1]), file = "timesNewman.txt", append = TRUE,sep = "\t")

})
 # for (i in c(1:N)){
 #   if (people[i] != 25){
 #     people[i] <- 0
 #   }
 # }

#find percolation
leftRoots <- sapply(c(1:sqrt(N)), findRoot)

