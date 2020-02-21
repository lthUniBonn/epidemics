
L = 10 # grid side length
N = L^2 # number of people
No = 10 #the No. of connections that will be added
people <- array(data=0,dim = N)
edges <- array(0, dim = c(N,N)) # this is 0 if there is no connection between the people and 1 if there is one
# it could be made so that the adges are bidirectional e.g. there is a connection between person 5 
#and person 3   and vice versa (edges[5,3] = edges[3,5]) or that there is a directional connection only (edges[5,3] = 1 and edges[3,5] = 0)
# in the first case only the upper triangle of the connections will be needed
# the Newman algorithm is bidirectional
made <- 0 # how many connections were already made
dof <- N*(N+1)/2 # degrees of freedom in symmetric matrix
usefulIndices <- array(0,dim=c(0,2)) # these are the possible connections

#now the indices that need to be checked are all in useful indices, this reduces the array size by N^2-N*(N+1)/2, for 100*100 grid is 4950 sites smaller
#connections <- usefulIndices[sample(nrow(usefulIndices),No,replace=FALSE),]
#connections is the list of random connections that will be made


calcDistance<- function(first, second){ # calculate distance between people
  rowfirst <- first %% sqrt(N)
  columnfirst <- first%/% sqrt(N)
  rowsecond <- second %% sqrt(N)
  columnsecond <- second %/% sqrt(N)
  distance <- sqrt((rowfirst-rowsecond)**2+(columnfirst - columnsecond)**2)
  return(distance)
  
}

addConnection <- function(){
  edges[connections[made+1,1],connections[made+1,2]] <<- 1 
  made <<- made +1
}
for(j in c(1:N)){
  for(i in c(j:N)){
    usefulIndices <- rbind(usefulIndices, c(j,i)) # filters out the lower triangle 
  }
}

for(s in c(1:nrow(usefulIndices))){
  if(calcDistance(usefulIndices[s,1],usefulIndices[s,2]) != 1){
    print('here')
    usefulIndices[s,1] <- 0
    usefulIndices[s,2] <- 0
  }
}

usefulIndices <- usefulIndices[-which(usefulIndices ==0),]

for(x in c(1:No-1)){
  addConnection()
  #merge(connections[made-1])
}