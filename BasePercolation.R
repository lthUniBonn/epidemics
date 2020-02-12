library('plot.matrix')
set.seed(3)
L <- 10 #Lattice side length
p <- 0.5

#initialise and plot lattice
lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
plot(lattice)

# what are <<-???
belongsToThisCluster <-array(0,dim = c(L,L))
neighbours <- array(0, dim=c(L,L,4))
# isItSpanning()

# isItSpanning <- function(){
#    for(i in c(1:L)){
#      if(lattice[i,1]== 1){
#        print("got here")
#       belongsToCluster[i,1] <<- 1
#     }
#   }
#   for(j in c(1:L)){
#       if(lattice[1,j]== 1){
#         belongsToCluster[1,j] <<- 1
#       }
#   }
# }


hasNeighbour <- function(i,j){
  isAtLeftEdge <<- FALSE
  isAtTopEdge <<- FALSE
  isAtBottomEdge <<- FALSE
  isAtRightEdge <<- FALSE
  
  if(i == L){# have to check this, otherwise produces out of bounds
    isAtBottomEdge <<- TRUE
    
  }
  if(j == L){
    isAtRightEdge <<- TRUE
    
  }
  if(i == 1){
    isAtTopEdge <<- TRUE
    
  }
  if(j == 1){
    isAtLeftEdge <<- TRUE
    
  }
  if(!isAtTopEdge){
  if(lattice[i-1,j]==1 ){
    neighbours[i,j,1] <<- 1 #has neighbour at top
    
  } 
  }
  if(!isAtRightEdge){
  if(lattice[i,j+1]==1){
    neighbours[i,j,2] <<- 1 #has neighbour at right
  } 
  }
  if(!isAtBottomEdge){
  if(lattice[i+1,j]==1){
    neighbours[i,j,3] <<- 1 #has neighbour at bottom
  } 
  }
  if(!isAtLeftEdge){
  if(lattice[i,j-1]==1){
    neighbours[i,j,4] <<- 1 #has neighbour at left
  } 
  }
}


# checkSurrounding <- function(i,j){
#   isAtLeftEdge <- FALSE
#   isAtTopEdge <- FALSE
#   isAtBottomEdge <- FALSE
#   isAtRightEdge <- FALSE
#   if(i == L || j == L){
#     print("percolates")
#     return(1)
#   }
#   if(i == 1){# have to check this, otherwise produces out of bounds
#     isAtTopEdge <<- TRUE
#   }
#   if(j == 1){
#     isAtLeftEdge <<- TRUE
#   }
#   if((lattice[i-1,j]==1) && (isAtTopEdge)){
#     belongsToThisCluster[i,j] <<- 1 #has neighbour at top
#   } 
#   if(lattice[i,j+1]==1){
#     belongsToThisCluster[i,j] <<- 1 #has neighbour at right
#   } 
#   if(lattice[i+1,j]==1){
#     belongsToThisCluster[i,j] <<- 1 #has neighbour at bottom
#   } 
#   if(lattice[i+1,j-1]==1 && isAtLeftEdge){
#     belongsToThisCluster[i,j] <<- 1 #has neighbour at left
#   } 
#   if(belongsToThisCluster[i,j] == 0){
#     return(0)
#   }
# }


move <-function(i,j, neighb){
  neiNew <- neighb
  print("-") 
  print(sum(neighb))
  print(sum(neighb[i,j,]))
  print(i)
  print(j)

  if(i == L){
    print("percolates")
    return(1)
  }
  if(j == L){
    print("percolates")
    return(1)
  }
  if(sum(neighb[i,j,])== 0){
    print("dead end")
    return(0)
  }
  else if(sum(neighb[i,j,])==1){  #if only one neighbour its a dangling end and should be removed before moving on
    print("one neigbour")
    if(neighb[i,j,1]==1){
      neiNew[i-1,j,3] <- 0 # remove dead end
      move(i-1,j, neiNew)
    }
    if(neighb[i,j,2]==1){
      neiNew[i,j+1,4] <- 0 # remove dead end
      move(i,j+1, neiNew)
    }
    if(neighb[i,j,3]==1){
      print("hier?")
      neiNew[i+1,j,1] <- 0 # remove dead end
      move(i+1,j, neiNew)
    }
    if(neighb[i,j,4]==1){
      neiNew[i,j-1,2] <- 0 # remove dead end
      move(i,j-1, neiNew)
    }
  }
  
  else if(sum(neighb[i,j,])>1){  #if more than one neighbour it can move anywhere
    print("multiple neigbour")
    if(neighb[i,j,1]==1){
      print("top")
      move(i-1,j, neighb)
    }
    else if(neighb[i,j,2]==1){
      print("right")
      move(i,j+1, neighb)
    }
    else if(neighb[i,j,4]==1){ #order changed, caught in top-bottom spiral --> check! nor left-right spiral!
      print("left")
      move(i,j-1, neighb)
    }
    else if(neighb[i,j,3]==1){
      print("bottom")
      move(i+1,j, neighb)
    }

  }
}

isAtLeftEdge <- FALSE
isAtTopEdge <- FALSE
isAtBottomEdge <- FALSE
isAtRightEdge <- FALSE
n0<-0
for(i in c(1:L)){
  for(j in c(1:L)){
    hasNeighbour(i,j)
    n0 <- n0+1
  }
}

#test with a certain start position
n <- 5
m <- 1
move(n,m, neighbours) # column i, row j
