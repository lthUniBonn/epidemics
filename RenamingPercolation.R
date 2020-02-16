library('plot.matrix')
library('data.tree')
library('DiagrammeR')
startTime <- proc.time()
#set.seed()

p <- 0.5


M <- 100
# M   t
# 10  0.04
# 10  0.06
# 10  0.08
# 10  0.07
# 20  0.04
# 20  0.05
# 20  0.07
# 20  0.03
# 30  0.09
# 30  0.07
# 30  0.07
# 30  0.09
# 40  0.16
# 40  0.15
# 40  0.15
# 40  0.08
# 50  0.19
#
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges
No <- 30 # number of lattices inspected per p 
NoPoints <- 20 #number of data points that should appear in the plot


lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
oldLattice <- lattice 
lattice[1,] <-0
lattice[L,] <-0
lattice[,1] <-0
lattice[,L] <-0
#plot(lattice) # plotten dauert EEEEEEEWIG
NoOccupied <- length(lattice[which(lattice==1)]) # useful quantity for cross checking

renamingVector <- data.frame() # rename second entry to first entry


clusterCounter <- 1
usedClusters <- numeric()
searchTime<- 0
renameTime <- 0



checkNeighbours <- function(i,j){
  if(lattice[i,j]==0){
    return(0)
  }
  top <- lattice[i-1,j]
  left <- lattice[i,j-1]
  if (left == 0 && top == 0){
    clusterCounter <<- clusterCounter +1
    lattice[i,j] <<- clusterCounter
    usedClusters <<- append(usedClusters,clusterCounter,length(usedClusters))
    return(1)
  }
  else if ((left == 0 && top != 0) || left == top){
    lattice[i,j] <<- top
  }
  else if (left != 0 && top == 0){
    lattice[i,j] <<- left
  }
  else{ 
    if(left<top){
      lattice[i,j] <<- left
      usedClusters[top]<<-0
      #-------this works fine, but scales quadratically with N------
      tmp <- proc.time()
      lattice[which(lattice == top)] <<- left # this takes about 90 % of the entire runTime
      renameTime <<- renameTime + proc.time() - tmp
      #-------------------------------------------------------------
      
      
    } else {
      lattice[i,j] <<- top
      usedClusters[left]<<-0
      #-------this works fine, but scales quadratically with N------
      tmp <- proc.time()
      lattice[which(lattice == left)] <<- top# this takes about 90 % of the entire runTime
      renameTime <<- renameTime + proc.time() - tmp
      #--------------------------------------------------------------
    }
  }
}

rows <- c(2:(L-1))
outer(rows,rows,Vectorize(checkNeighbours)) # runs checkNeighbour

endTime <- proc.time()
runTime <- endTime-startTime

# #-------------------------------------------------------------------------------
# #------------------------------ Calculate Observables---------------------------
# #-------------------------------------------------------------------------------

#--------------------------------Cluster Sizes------------------------------------

sizes <- table(as.vector(lattice))
sizes[names(sizes)==0] <- 0
hist(sizes)
largestCluster <- max(sizes)
# usedClusters <- unique(as.vector(lattice))
# for(x in usedClusters){
#   
# }
#write(x = c(M**2, runTime[1]), file = "timesRenaming.txt", append = TRUE,sep = "\t")
