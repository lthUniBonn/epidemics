library('plot.matrix')
set.seed(5)

p <- 0.5

#initialise and plot lattice
M <- 20 #size of usable array
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges


lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
lattice[1,] <-0
lattice[L,] <-0
lattice[,1] <-0
lattice[,L] <-0
plot(lattice)

#make list of possible starting nodes

startTop <- which(lattice[2,]==1)
startLeft <- which(lattice[,2]==1)+1 # +1 necessary to avoid starting at [2,2] twice



# strat parameters
i <- 2
j <- 3
recursionDepth <- 0
move <- function(i,j,left,top){
  recursionDepth <<- recursionDepth +1
  print(i)
  print(j)
  if(i== L-1 && top==TRUE){
    print("percolates from top to bottom")
    return(1)
  }
  if(j== L-1 && left==TRUE){
    print("percolates from left to right")
    return(1)
  }
  lattice[i,j] <<- lattice[i,j] +1 #increment this site so that it shows we were here bofore
  moveThatWay<-dirDecision(i,j)
  if(moveThatWay==0){
    return(0)
  }
  if(moveThatWay == 1){
    print("Top")
    return(move(i-1,j,left,top))
  }
  if(moveThatWay == 2){
    print("Left")
    return(move(i,j-1,left,top))
  }
  if(moveThatWay == 3){
    print("Right")
    return(move(i,j+1,left,top))
  }
  if(moveThatWay == 4){
    print("Bottom")
    return(move(i+1,j,left,top))
  }
  
}

dirDecision <- function(i,j){ # this evaluates where the next step should ge if multiple are possible
  visitedTop <- lattice[i-1,j]
  visitedLeft <- lattice[i,j-1]
  visitedRight <- lattice[i,j+1]
  visitedBot <- lattice[i+1,j]
  neigbs <- 4 # counts how many neibghours the node has 
  if(visitedTop == 0){
    visitedTop <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }
  if(visitedLeft == 0){
    visitedLeft <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }
  if(visitedRight == 0){
    visitedRight <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }
  if(visitedBot == 0){
    visitedBot <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }
  visitList <- c(1/visitedTop,1/visitedLeft,1/visitedRight,1/visitedBot)
  if(neigbs==1){
    lattice[i,j]<<- 0 # burn dangling end
  }
  if(neigbs == 0){
    print("no further connections")
    return(0)
  } else { return(which.max(x = visitList)) }
  
  
}



left <- FALSE
top <- FALSE
#while(percolation found or excluded){
#  for(i in startTop){
  
    if( i == 2){
      top <<- TRUE
    } else if(i!=1){
      top <<- FALSE
    }
    if( j == 2){
      left <<- TRUE
    } else if(j!=1){
      left <<- FALSE
    }
    if(lattice[i,j]==1){
      move(i,j, left, top) # column i, row j left, top are bool when start is at left,top
    } else {print("This node is not infected")}
  #}
  #for()
#}




