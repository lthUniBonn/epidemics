library('plot.matrix')
startTime <- proc.time()
#set.seed(1)

p <- 0.5

#initialise and plot lattice
M <- 150 #size of usable array
# runTime rises faster than the sqaure of M
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges
No <- 30 # number of lattices inspected per p 
# No is roughly proportional to runTime
NoPoints <- 20 #number of data points that should appear in the plot
#0: not infected ; >=1 infected, (2) visited
 lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
# lattice[1,] <-0
# lattice[L,] <-0
# lattice[,1] <-0
# lattice[,L] <-0
#plot(lattice)
#plot(lattice[c(1:20), c(1:20)])

#make list of possible starting nodes

startTop <- which(lattice[2,]==1)
startLeft <- which(lattice[,2]==1)+1 # +1 necessary to avoid starting at [2,2] twice



# start parameters
n <- 2
m <- 2
moves <- 0
percolations <- 0
doesItPercolate <- FALSE
moveHist <- numeric(length = 1)
percResults <- numeric(length = NoPoints)
pResults <- numeric(length = NoPoints)
#while(percolation found or excluded){
#  for(i in startTop){




fullClusterChecked <- function(i,j,i0, j0){ 
  #return 0 if checking of cluster needs to continue
  # retrun 1 if percolation was found
  #return 2 if cluster is finished and not percolating
 # print("---")
  #print(i)
  #print(j)
  neighbours <- c(lattice[i-1,j], lattice[i,j-1], lattice[i+1,j],lattice[i,j+1])
  #
  if(i == (L-1) && i0==2){
    #print("percolates from top to bottom")
    doesItPercolate <<- TRUE
    return(1)
  }
  else if(j== (L-1) && j0==2){
    #print("percolates from left to right")
    doesItPercolate <<- TRUE
    return(1)
  }
  else if(i == i0 && j == j0 && !(1 %in% neighbours) ){
    #print("full cluster ckecked - no percolation found")
    return(2) 
  }
  else {return(0)}
}

saveHistory <- function (step){
  if (moveHist[length(moveHist)] == (step+2) %% 4){
    moveHist <<- moveHist[-length(moveHist)]
  }
  else{moveHist <<- c(moveHist, step)}
}

move <- function(i,j, moveThatWay){
  #print(i)
  #print(j)
  #print(moveThatWay)
  if(moveThatWay == 1){
    #print("Top")
    saveHistory(1)
    return(c(i-1,j))
  }
  if(moveThatWay == 2){
    #print("Left")
    saveHistory(2)
    return(c(i,j-1))
  }
  if(moveThatWay == 3){
    #print("Bottom")
    saveHistory(3)
    return(c(i+1,j))
  }
  if(moveThatWay == 0){
    #print("Right")
    saveHistory(0)
    return(c(i,j+1))
  }
}
#direction:  top:1 left:2 bottom:3 right:0
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
  if(visitedBot == 0){
    visitedBot <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }
  if(visitedRight == 0){
    visitedRight <- 100000000 # make sure that this is not visited if 0
    neigbs <- neigbs -1
  }

  visitList <- c(visitedRight, visitedTop,visitedLeft,visitedBot)
  possibleList <- visitList[which(visitList != 100000000)]
  possibleMove <- which(visitList != 100000000)-1

  if(neigbs==1){
    #print("burn dangling end")
    #lattice[i,j] <<- 0 # burn dangling end
    return(which.min(visitList)-1)
  }
  else if(neigbs == 0){
    return(-1)
  } 
  else if (length(unique(possibleList)) == 1 && unique(possibleList) == 2 && 
           ((moveHist[length(moveHist)]+2)%%4) %in% possibleMove){
    #lattice[i,j] <<- 3 # burn loop end
    #print("burn loop end")
  #  if {
      return((moveHist[length(moveHist)]+2) %% 4) # move backwards
  #  }
  #  else {
  #    return(which.min(visitList))
  #  }
    }
  else {
    return(which.min(visitList)-1) }
}



#------------------------------------------------------------------------------------------------------
# rekursion auflösen?? while loop stattdessen?
clusterParsing <- function(i0,j0){
  pos <- c(i0,j0)
  if(lattice[pos[1],pos[2]]==0){
    #print("This node is not infected")
    return(0)
  }
  else{
    while (TRUE) {
      #print(pos)
      lattice[pos[1],pos[2]] <<- 2 #increment this site so that it shows we were here bofore
      moves <<- moves +1
      #if(moves == 10){break} #only for debugging
      moveThatWay <- dirDecision(pos[1],pos[2])
      #print(moveThatWay)
      if (moveThatWay == -1){
        #print("no further connections")
        #return(1)
        ##lattice[pos[1],pos[2]] <<- 0 # burn single site
        return(FALSE)
        break
      }
      else{
        pos <- move(pos[1],pos[2],moveThatWay)
        #print(pos)
      }
      if (fullClusterChecked(pos[1],pos[2],i0,j0) == 1 ){ # this cluster percolates
        return(TRUE)  
        break
      }else if (fullClusterChecked(pos[1],pos[2],i0,j0) == 2 ){ # this cluster does not percolate
        return(FALSE)  
        break
      }  
    }
  }
}
# moves speichern -> so weit zurück bis wieder 1? --> burn ausgang aus loop? 
# zurück bis start: ganzes cluster, keine percolation --> alle 2en removen
# 2 an beiden seiten: percolating
for(p in seq(0,1,by = 1/NoPoints)){
  for(l in c(1:No)){
    lattice <<- array(rbinom(n = L*L,1,p), dim = c(L,L))
    lattice[1,] <-0
    lattice[L,] <-0
    lattice[,1] <-0
    lattice[,L] <-0
    startTop <- which(lattice[2,]==1)
    startLeft <- which(lattice[,2]==1)
    #print("checking lattice")
    #print(l)
    doesItPercolate <- FALSE
    while(length(startTop)!= 0 && !doesItPercolate){ # continue scanning this cluster until no further starting points are available or a percolation is found
      #print("got here")
      doesItPercolate <- clusterParsing(2,startTop[1])
      #print(doesItPercolate)
      if(doesItPercolate){
        #print("percolation counter up")
        percolations <<- percolations +1
      } else {
        #remove all sites previously visited
        #print("removing all sites with 2")
        lattice[which(lattice==2)] <-0
        startTop <- which(lattice[2,]==1)
        if(length(startTop)== 0){
          break
        }
      }
    }
    while(length(startLeft)!= 0 && !doesItPercolate){ # continue scanning this cluster until no further starting points are available or a percolation is found
      
      doesItPercolate <- clusterParsing(startLeft[1],2)
      #print(doesItPercolate)
      if(doesItPercolate){
        #print("percolation counter up")
        percolations <<- percolations +1
      } else {
        #remove all sites previously visited
        #print("removing all sites with 2")
        lattice[which(lattice==2)] <-0
        startLeft <- which(lattice[,2]==1)
        if(length(startLeft)== 0){
          #print("no percolation is this lattice at all")
        }
      }
    }
  }
  percResults[p*NoPoints] <- percolations
  pResults[p*NoPoints] <- p
  percolations <- 0
}
plot(x = pResults, y=percResults/No,xlab = "Infection Probability", ylab = "Percolation Probability")
endTime <- proc.time()
runTime <- endTime-startTime