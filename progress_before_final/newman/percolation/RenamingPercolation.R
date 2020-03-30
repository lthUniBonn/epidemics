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


findDuplicateInitials <- function(initial){
  
  duplicateInitials <- labelVector[which(labelVector[,2] == initial), ]
  
  # print(duplicateInitials)
  if (nrow(duplicateInitials) >1 ){
    
    # sort by target
    #labelVector <<- labelVector[order(labelVector[,1], labelVector[,2]),]
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


# vector dupes raus

labelVector <- labelVector[which(!duplicated(labelVector)),]
# sort by target
labelVector <- labelVector[order(labelVector[,1], labelVector[,2]),]
#print("C")
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
#tmp <- proc.time()
while (TRUE) {
  
  lapply(existingInitialsVec,FUN=findDuplicateInitials)
  
  #print(which(newLabelVector[,1]==newLabelVector[,2]))
  
  y <- 1
  while(y <= nrow(newLabelVector)){
    #print(y)
    if(newLabelVector[y,1]==newLabelVector[y,2]){
      newLabelVector <- newLabelVector[-y,]
    }
    y <- y +1
  }
  existingInitialsVec <- newLabelVector[,2]
  existingInitialsVec <- existingInitialsVec[which(!duplicated(existingInitialsVec))]
  existingInitialsVec <- sort(existingInitialsVec, decreasing = TRUE)
  #print("loop")
  #print(newLabelVector[which(duplicated(newLabelVector[,2])),])
  #wenn links gleich rechts raus!! 
  
  
  
    
  if (nrow(newLabelVector[which(duplicated(newLabelVector[,2])),]) == 0){
    break
  }
  labelVector <- newLabelVector
  newLabelVector <- data.frame()
}
for (i in c(1:nrow(newLabelVector))){# vielleicht gute idee das auch in die fors
  initial <- newLabelVector[i,2]
  target <- newLabelVector[i,1]
  newLabelVector[which(newLabelVector[,1] == initial),1] <- target
}
#timeReorganizing <- proc.time() - tmp
newLabelVector <- newLabelVector[order(newLabelVector[,2], newLabelVector[,1], decreasing = TRUE),]

for(s in c(1:nrow(newLabelVector))){
  lattice[which(lattice == newLabelVector[s,2])] <- newLabelVector[s,1]
}

#for (i in c(1:(L**2))){
# for (i in c(1:L)){
#   for (j in c(1:L)){
#     if (lattice[i] != 0){
#       found <- newLabelVector[which(newLabelVector[,2] == lattice[i]),1]
#       if (length(found) != 0){
#         lattice[i] <- found
#    }
#   }
#   }
# }

#endTime <- proc.time()
#runTime <- endTime-startTime
# #-------------------------------------------------------------------------------
# #------------------------------ Calculate Observables---------------------------
# #-------------------------------------------------------------------------------

#--------------------------------Cluster Sizes------------------------------------

sizes <- table(as.vector(lattice))
sizes[names(sizes)==0] <- 0
#hist(sizes)
largestCluster <- max(sizes)
# usedClusters <- unique(as.vector(lattice))
# for(x in usedClusters){
#   
# }
})
#write(x = c(M**2, runTime[1]), file = "timesRenamingAfter.txt", append = TRUE,sep = "\t")
profile
