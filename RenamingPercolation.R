library('plot.matrix')
library('data.tree')
library('DiagrammeR')
library('dplyr')
#library('tidyverse')
startTime <- proc.time()
set.seed(1)

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


clusterCounter <- 1
usedClusters <- numeric()
searchTime <- 0
renameTime <- 0
labelVector <- data.frame()
newLabelVector <- data.frame()


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
      #-------this works fine, but scales quadratically with N------
      #tmp <- proc.time()
      #lattice[which(lattice == top)] <<- left # this takes about 90 % of the entire runTime
      #renameTime <<- renameTime + proc.time() - tmp
      #-------------------------------------------------------------
      labelVector <<- rbind(labelVector, c(left, top))
      
    } else {
      lattice[i,j] <<- top
      usedClusters[left]<<-0
      #-------this works fine, but scales quadratically with N------
      #tmp <- proc.time()
      #lattice[which(lattice == left)] <<- top# this takes about 90 % of the entire runTime
      #renameTime <<- renameTime + proc.time() - tmp
      #--------------------------------------------------------------
      labelVector <<- rbind(labelVector, c(top, left))
    }
  }
}


findDuplicateInitials <- function(initial = 90){
  duplicateInitials <- labelVector[which(labelVector[,2] == initial), ]
  if (nrow(duplicateInitials) >1 ){
    # sort by target
    labelVector <<- labelVector[order(labelVector[,1], labelVector[,2]),]
    #remove duplicate initials 
    #labelVector <<- labelVector[-which(labelVector[,2] == initial), ]
    target <- duplicateInitials[1,1]
    for (i in c(2:length(duplicateInitials))){
      duplicateInitials[i,2] <- duplicateInitials[i,1]
      duplicateInitials[i,1] <- target
      }
  }
  newLabelVector <<- rbind(newLabelVector, duplicateInitials)
}

rows <- c(2:(L-1))
outer(rows,rows,Vectorize(checkNeighbours)) # runs checkNeighbour

endTime <- proc.time()
runTime <- endTime-startTime

# vector dupes raus

labelVector <- labelVector[which(!duplicated(labelVector)),]
# sort by target
labelVector <- labelVector[order(labelVector[,1], labelVector[,2]),]

#reduce chains of targeting

for (i in c(1:nrow(labelVector))){
  initial <- labelVector[i,2]
  target <- labelVector[i,1]
  labelVector[which(labelVector[,1] == initial),1] <- target
}
labelVector <- labelVector[which(!duplicated(labelVector)),]
labelVector <- labelVector[order(labelVector[,1], labelVector[,2]),]
labelVectorFull <- labelVector
# >= dupes in startClusters: start 
#for (i in c(clusterCounter:1)){findDuplicateInitials(i)}
existingInitialsVec <- labelVector[,2]
existingInitialsVec <- existingInitialsVec[which(!duplicated(existingInitialsVec))]
existingInitialsVec <- sort(existingInitialsVec, decreasing = TRUE)
lapply(existingInitialsVec,FUN=findDuplicateInitials)

newLabelVector <- newLabelVector[order(newLabelVector[,2], newLabelVector[,1], decreasing = TRUE),]

for(s in c(1:nrow(newLabelVector))){
  print(s)
  lattice[which(lattice == newLabelVector[s,2])] <- newLabelVector[s,1]
}
#findDuplicateInitials(90)
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
  