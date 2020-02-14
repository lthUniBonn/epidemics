library('plot.matrix')
library('data.tree')
startTime <- proc.time()
#set.seed(1)

p <- 0.5

#initialise and plot lattice
M <- 60 #size of usable array
# runTime rises roughly as a sqaure of M
L <- M+2 #expand array by 1 in each direction to make it uneccesary to inculde special cases for edges
No <- 30 # number of lattices inspected per p 
# No is roughly proportional to runTime
NoPoints <- 20 #number of data points that should appear in the plot
#0: not infected ; >=1 infected, (2) visited

lattice <- array(rbinom(n = L*L,1,p), dim = c(L,L))
 lattice[1,] <-0
 lattice[L,] <-0
 lattice[,1] <-0
 lattice[,L] <-0
plot(lattice)
#plot(lattice[c(1:20), c(1:20)])

# start parameters
n <- 2
m <- 2

checkNeigbours <- function(i,j, clusterCounter){
  top <- lattice[i-1][j]
  left <- lattice[i][j-1]
  
  if (left == 0 && top == 0){
    newTree(i,j, clusterCounter) #make new tree and increment counter, node gets counter value
    # ??? add node via numbers not names? acme$children[[1]]$children[[2]]$name
    #problem: cannot add/make new nodes without new name in function?? 
    
    #wie genau finden wir die node zu der wir schrieben wollen <-> nur top node verbindungen? 
    #schreiben wir die zahl in den Tree? 
    # tree list position entspricht counter zahl? 
    # new tree ansprechen über listenposition ( haben alle selben namen) 
    # tree mit find node mit position --> dran schreiben 
    return(0)
  }
  else if (left == 0 && top != 0){
    #read counter top node 
    #add node to top node with top counter
  }
  else if (left != 0 && top == 0){
    
  }
  else {
    
  }
  
}









endTime <- proc.time()
runTime <- endTime-startTime