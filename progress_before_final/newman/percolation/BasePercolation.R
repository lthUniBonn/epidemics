library('plot.matrix')
set.seed(5)

p <- 0.9

#initialise and plot lattice
M <- 100 #size of usable array
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges

#0: not infected ; >=1 infected, (-1) visited
lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
lattice[1,] <-0
lattice[L,] <-0
lattice[,1] <-0
lattice[,L] <-0
plot(lattice)
plot(lattice[c(1:20), c(1:20)])

#make list of possible starting nodes

startTop <- which(lattice[2,]==1)
startLeft <- which(lattice[,2]==1)+1 # +1 necessary to avoid starting at [2,2] twice



# start parameters
i <- 2
j <- 11
recursionDepth <- 0

move <- function(i,j,left,top){
  recursionDepth <<- recursionDepth +1
  print("---")
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
  lattice[i,j] <<- 2 #increment this site so that it shows we were here bofore
  moveThatWay <- dirDecision(i,j)
  if(moveThatWay == 0){
    return(0)
  }
  if(moveThatWay == 1){
    print("Top")
    backMove <<- 4
    return(move(i-1,j,left,top))
  }
  if(moveThatWay == 2){
    print("Left")
    backMove <<- 3
    return(move(i,j-1,left,top))
  }
  if(moveThatWay == 3){
    print("Right")
    backMove <<- 2
    return(move(i,j+1,left,top))
  }
  if(moveThatWay == 4){
    print("Bottom")
    backMove <<- 1
    return(move(i+1,j,left,top))
  }
  
}
#direction:  top:1 left:2 right:3 bottom:4
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
  visitList <- c(visitedTop,visitedLeft,visitedRight,visitedBot)
  possibleList <- visitList[which(visitList != 100000000)]
  possibleIndex <- which(visitList != 100000000)

  if(neigbs==1){
    print("burn danglnig end")
    lattice[i,j] <<- 0 # burn dangling end
    return(which.min(visitList))
  }
  else if(neigbs == 0){
    print("no further connections")
    lattice[i,j] <<- 0 # burn single site
    return(0)
  } 
  else if (length(unique(possibleList)) == 1 && unique(possibleList) != 1){
    lattice[i,j] <<- 0 # burn loop end
    print("burn loop end")
    if (backMove %in% possibleIndex){
      return(backMove)
    }
    else {
      return(which.min(visitList))
    }
    }
  else {return(which.min(visitList)) }
}



left <- FALSE
top <- FALSE
backMove <- 0
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

# rekursion auflösen?? while loop stattdessen?
# moves speichern -> so weit zurück bis wieder 1? --> burn ausgang aus loop? 
# zurück bis start: ganzes cluster, keine percolation --> alle 2en removen
# 2 an beiden seiten: percolating



