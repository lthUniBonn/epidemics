library('plot.matrix')
library('data.tree')
library('DiagrammeR')
library('dplyr')
library('data.table')
#library('tidyverse')
library('profvis')
#startTime <- proc.time()
#profile <- profvis({
set.seed(1)

p <- 0.5


M <- 100
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
newLabelVector <- data.table(x=rep(0,(M**2)),y=rep(0,M**2)) # is the maximum length of this really M**2?


usedLength <- function(dataTable){ # returns the index of the last row which is not 0
  
  if(length(dataTable)==0){
    return(0)
    break
    }
  for(i in c(1:length(dataTable))){
      if(dataTable[i,1]== 0){
        return(i-1)
        break
      }
  }
  return(length(dataTable))
}

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
  
  if (nrow(duplicateInitials) >1 ){
    
    # sort by target
    #labelVector <<- labelVector[order(labelVector[,1], labelVector[,2]),]
    target <- duplicateInitials[1,1]
    for (i in c(2:length(duplicateInitials))){
      duplicateInitials[i,2] <- duplicateInitials[i,1]
      duplicateInitials[i,1] <- target
    }
    
  }
  
  #newLabelVector <<- rbind(newLabelVector, duplicateInitials) # this is replace by next 4 lines 
  
  end <- usedLength(newLabelVector) # last entry which is not 0
  
  
  set(newLabelVector,i = c((end+1):(end+nrow(duplicateInitials))),1L,value = duplicateInitials[,1])#c(end:(end+nrow(duplicateInitials)) chooses the rows in which to append to
  set(newLabelVector,i = c((end+1):(end+nrow(duplicateInitials))),2L,value = duplicateInitials[,2])
  
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
stop()
labelVector <- labelVector[which(!duplicated(labelVector)),]
labelVector <- labelVector[order(labelVector[,1], labelVector[,2]),]
labelVectorFull <- labelVector
# >= dupes in startClusters: start 

existingInitialsVec <- labelVector[,2]
existingInitialsVec <- existingInitialsVec[which(!duplicated(existingInitialsVec))]
existingInitialsVec <- sort(existingInitialsVec, decreasing = TRUE)
print("D")

while (TRUE) {
  
  lapply(existingInitialsVec,FUN=findDuplicateInitials)
  print("jsda")
  
  
  
  y <- 1
  while(y <= usedLength(newLabelVector)){
    
    
    if(newLabelVector[y,1]==newLabelVector[y,2] && newLabelVector[y,1]!= 0){ # remove entries which have identical initial and final
      
      newLabelVector <- newLabelVector[-y,]
    }
    y <- y +1
    
  }
  
  existingInitialsVec <- unlist(newLabelVector[,2])
  existingInitialsVec <- existingInitialsVec[which(!duplicated(existingInitialsVec))]
  existingInitialsVec <- sort(existingInitialsVec, decreasing = TRUE)
  
  
  
  #wenn links gleich rechts raus!! 
  
  
  
  
  if (nrow(newLabelVector[which(duplicated(newLabelVector[c(1:usedLength(newLabelVector)),2])),]) == 0){
    #want: how many duplicates are there which are not 0
    print("C")
    break
  }
  labelVector <- newLabelVector
  newLabelVector <- data.table(x=rep(0,(M**2)),y=rep(0,M**2))
}

print(class(labelVector))
stop()
for (i in c(1:nrow(labelVector))){# vielleicht gute idee das auch in die fors
  
  initial <- labelVector[i,2]
  target <- labelVector[i,1]
  labelVector[which(labelVector[,1] == initial),1] <- target
  
}

labelVector <- labelVector[order(labelVector[,2], labelVector[,1],decreasing = TRUE),]
print("E")

#timeReorganizing <- proc.time() - tmp
#print("D")
#labelVector <- setorder(x = labelVector,-1)

#print("C")

for(s in c(1:nrow(labelVector))){
  lattice[which(lattice == labelVector[s,2])] <- labelVector[s,1]
  #print(s)
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
#})
#write(x = c(M**2, runTime[1]), file = "timesRenamingAfter.txt", append = TRUE,sep = "\t")
profile
