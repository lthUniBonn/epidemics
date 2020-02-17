library('plot.matrix')
library('data.tree')
library('DiagrammeR')
startTime <- proc.time()
set.seed(1)

p <- 0.5

#initialise and plot lattice
M <- 100 #size of usable array
# M   t
# 10  0.1
# 20  0.25
# 30  0.37
# 40  0.63
# 50  0.94
# 60  1.38
# 80  2.09
# 100 3.44
# 150 7.08
# 200 14.97
# 500 85.5
# rougly quadratic in time
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges
No <- 30 # number of lattices inspected per p 
# No is roughly proportional to runTime
NoPoints <- 20 #number of data points that should appear in the plot
#0: not infected ; >=1 infected, (2) visited
nodeTime <- 0
lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
oldLattice <- lattice # use this to plot the initial lattice afterwards
NoOccupied <- length(lattice[which(lattice==1)])
#lattice <- array(c(0,0,0,0,0,1,1,0,0,0,1,0,0,1,1,1,0,0,0,0), dim=c(4,5))
 lattice[1,] <-0
 lattice[L,] <-0
 lattice[,1] <-0
 lattice[,L] <-0
plot(lattice)
#plot(lattice[c(17:40), c(39:45)])

# start parameters

treeList <- list()
clusterCounter <- 1
usedClusters <- numeric()


newTree <- function(i, j){
  treeList[[clusterCounter]] <<- Node$new(paste0(i,",",j))
  #treeList[[clusterCounter]]$AddChild(paste0(i,",",j))
  
  # print(treeList[[clusterCounter]])
  # parent <- Climb(treeList[[clusterCounter]], name = paste0(i,",",j)) # finds which node is the parent
  # print(parent)
  # parent$AddChild(paste0(2,",",3))
  # print(treeList[[clusterCounter]])
  return(0)
}

addNode <- function(i0, j0, i, j, clusterID){ # check runtime for findNode, might be long
  #create a new child for the node that already exists at i0, j0
  #name this parent node by its lattice i, j
  tempTime <- proc.time()
  parent <- FindNode(treeList[[clusterID]], name = paste0(i0,",",j0)) # finds which node is the parent
  nodeTime <<- nodeTime + proc.time() - tempTime
  #print(parent)
  parent$AddChild(paste0(i,",",j))
  
}

mergeTrees <- function(leftID, topID){
  
  if(leftID<topID){
    treeList[[leftID]]$AddChildNode(treeList[[topID]])
    
  } else {
    treeList[[topID]]$AddChildNode(treeList[[leftID]])
  }
  
}

checkNeighbours <- function(i,j){
  #print(i)
  #print(j)
  #print(lattice[i,j])
  if(lattice[i,j]==0){
    return(0)
  }
  if(lattice[i,j]==0){
    print("fail")
  }
  top <- lattice[i-1,j]
  left <- lattice[i,j-1]
  if (left == 0 && top == 0){
    clusterCounter <<- clusterCounter +1
    lattice[i,j] <<- clusterCounter
    usedClusters <<- append(usedClusters,clusterCounter,length(usedClusters))
    newTree(i,j) 
   # print(treeList[[clusterCounter]])
    #clusterCounter strats at 1
    # ??? add node via numbers not names? acme$children[[1]]$children[[2]]$name
    #wie genau finden wir die node zu der wir schrieben wollen <-> nur top node verbindungen? 
    # tree mit find node mit position --> dran schreiben 
    return(1)
  }
  else if ((left == 0 && top != 0) || left == top){
    addNode(i-1, j, i, j, top)
    lattice[i,j] <<- top
    #print(treeList[[top]])
    #read counter top node 
    #add node to top node with top counter
  }
  else if (left != 0 && top == 0){
    addNode(i, j-1, i, j, left)
    lattice[i,j] <<- left
    #print(treeList[[left]])
  }
  else{ #merge trees larger number goes to smaller number
    if(left<top){
      addNode(i,j-1,i,j,left)
      lattice[i,j] <<- left
      treeList[[left]]$AddChildNode(treeList[[top]])
      usedClusters[top]<<-0
      
    } else {
      addNode(i-1,j,i,j,top)
      lattice[i,j] <<- top
      treeList[[top]]$AddChildNode(treeList[[left]])
      usedClusters[left]<<-0
    }
    
    
    
  }
  
}

# for(n in c(2:(L-1))){
#   for(m in c(2:(L-1))){
#   checkNeighbours(m,n)
#   }
# }
rows <- c(2:(L-1))
# 
# indexMat <- array(c(2:(L-1)))
#apply(X = c(1:L),MARGIN = 1,FUN = checkNeighbours,j=2)
outer(rows,rows,Vectorize(checkNeighbours))

endTime <- proc.time()
runTime <- endTime-startTime

#-------------------------------------------------------------------------------
#------------------------------ Calculate Observables---------------------------
#-------------------------------------------------------------------------------

#------------------------largest Cluster and order Parameter--------------------
counts <- numeric()
for(s in (which(usedClusters !=0))){
  #print(s)
  if(s != length(treeList)){
    counts[s] <- treeList[[s+1]]$totalCount
  }
}
for(s in (which(usedClusters ==0))){
  counts[s] <- 0
}
largestCluster <- max(counts)
orderPara <- largestCluster/NoOccupied


largestIndex <- which.max(counts)
plot(treeList[[largestIndex+1]])
#debugging attempts: for M=100 the plotted maximal tree yields a weird display of the following two points:
#oldLattice[25,32] <-2
#oldLattice[18,40] <-2
#but these ponits are part of the large cluster, so it is only a problem in displaying the tree, not in the program
#plot(oldLattice[c(1:30), c(1:40)]) 

#plot(oldLattice)

#----------------------spanning Cluster------------------------------------------


