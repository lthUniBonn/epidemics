library('plot.matrix')
library('data.tree')
library('DiagrammeR')
library('dplyr')
#library('tidyverse')
library('profvis')
#startTime <- proc.time()
profile <- profvis({
  set.seed(1)
  
  p <- 0.5
  
  
  M <- 500
  # M   t
  # 
  L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges
  No <- 30 # number of lattices inspected per p 
  NoPoints <- 20 #number of data points that should appear in the plot
  
  
  lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
  lattice[1,] <-0
  lattice[L,] <-0
  lattice[,1] <-0
  lattice[,L] <-0
  #plot(lattice) # plotten dauert EEEEEEEWIG
  #NoOccupied <- length(lattice[which(lattice==1)]) # useful quantity for cross checking
  
  
  clusterCounter <- 1
  usedClusters <- numeric()
  checkNeighbours <- function(i,j){
    if(lattice[i,j]==0){
      return(0)
    }
    top <- lattice[i-1,j]
    left <- lattice[i,j-1]
    if (left == 0 && top == 0){
      clusterCounter <<- clusterCounter +1
      lattice[i,j] <<- clusterCounter
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
        #-------this works fine, but scales quadratically with N------
        #tmp <- proc.time()
        lattice[which(lattice == top)] <<- left # this takes about 90 % of the entire runTime
        #renameTime <<- renameTime + proc.time() - tmp
        #-------------------------------------------------------------
        
        
      } else {
        lattice[i,j] <<- top
        usedClusters[left]<<-0
        #-------this works fine, but scales quadratically with N------
        #tmp <- proc.time()
        lattice[which(lattice == left)] <<- top# this takes about 90 % of the entire runTime
        #renameTime <<- renameTime + proc.time() - tmp
        #--------------------------------------------------------------
        
      }
    }
  }
  
  
  
  
  rows <- c(2:(L-1))
  
  outer(rows,rows,Vectorize(checkNeighbours)) # runs checkNeighbour
  
  
  
})
#write(x = c(M**2, runTime[1]), file = "timesRenamingAfter.txt", append = TRUE,sep = "\t")
profile
