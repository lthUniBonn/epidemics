library('plot.matrix')
library('gmp')
library('profvis')
library('Brobdingnag')


patientZero <- function(){
  zeroId <- sample(x = which(sDistribution != 0),size = 1) #choose one of the susceptibles to be patient 0
  infected[zeroId] <<- TRUE
  #latticeID <- which(lattice[,,1] == zeroId, arr.ind = TRUE)
  #lattice[latticeID[1], latticeID[2],2] <<- 1 
}

findRecBool <- function(t){
  return(t>=rnorm(length(t), mean = avgRecoveryTime, sd = sdRecoveryTime))
}

timesteps <- function(){
  # find infected people
  infPeople <- which(infected==TRUE)
  if (length(infPeople) == 0){
    #("no more infected")
    return(TRUE)
  }
  
  
  # find infected connections
  # only get bonds where exactly one site is infected
  a <- possConn[,1] %in% infPeople
  b <- possConn[,2] %in% infPeople
  infConn <- possConn[which(xor(a, b)),c(1,2)]
  
  # split infConn for convenience
  #infBondProb <- infConn[,3]
  possiblyInf1 <- infConn[,1]
  possiblyInf2 <- infConn[,2]
  
  # test for bond probability
  #bondCheck <- runif(n = nrow(infConn))
  #possiblyInf1 <- possiblyInf1[which(infBondProb >= bondCheck)]
  #possiblyInf2 <- possiblyInf2[which(infBondProb >= bondCheck)]
  
  susc1 <- sDistribution[possiblyInf1]/sDistFactor
  susc2 <- sDistribution[possiblyInf2]/sDistFactor
  
  # test for susceptibility
  newlyInf1 <- possiblyInf1[which(susc1 >= runif(n = length(possiblyInf1)))] 
  newlyInf2 <- possiblyInf2[which(susc2 >= runif(n = length(possiblyInf2)))]
  
  #recovery // remove people
  recoveryCheck <- findRecBool(infectionTime[infected])
  recPeople <- infPeople[recoveryCheck]
  infectionTime[recPeople] <<- 0
  
  # recovered people are now immune
  sDistribution[recPeople] <<- 0
  
  
  # check if people who are recovering are newly infected and dismiss these infections
  #newlyInf1 <- setdiff(newlyInf1,recPeople)
  #newlyInf2 <- setdiff(newlyInf2,recPeople)
  
  
  
  primaryInfected <- length(infPeople)
  
  infected[newlyInf1] <<- TRUE
  infected[newlyInf2] <<- TRUE
  totalInfected <- length(which(infected == TRUE)) #newly infected (including recovering)
  infected[recPeople] <<- FALSE
  
  #increase infection time 
  infectionTime[infected] <<- infectionTime[infected] + 1
  
  #R0 <- (totalInfected-primaryInfected)/primaryInfected
  #statistical participation in disease spread for this timestep
  
  recovered[recPeople] <<- 1 #only take R0 measurement of recovered (for a total disease cycle and not partial timesteps)
  #R0OverInfectiousPeriod[infPeople] <<- R0OverInfectiousPeriod[infPeople] + R0
 
  #if(anyNA(evalR0[statRun,])){
  #  if (sum(recovered)>=checkR0Here){
  #     R0Mean <<- mean(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])#!! check!! 
  #     R0Sd <<- sd(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])
  #   }
  #   
  # }
  # if(anyNA(evalR02[statRun,])){
  #   if (sum(recovered)>=2){
  #     R0Mean2 <<- mean(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])#!! check!! 
  #     R0Sd2 <<- sd(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])
  #   }
  #   
  # }
  # if(anyNA(evalR03[statRun,])){
  #   if (sum(recovered)>=3){
  #     R0Mean3 <<- mean(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])#!! check!! 
  #     R0Sd3 <<- sd(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])
  #   }
  #   
  # }
  # if(anyNA(evalR04[statRun,])){
  #   if (sum(recovered)>=4){
  #     R0Mean4 <<- mean(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])#!! check!! 
  #     R0Sd4 <<- sd(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])
  #   }
  #   
  # }
  # if(anyNA(evalR07[statRun,])){
  #   if (sum(recovered)>=7){
  #     R0Mean7 <<- mean(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])#!! check!! 
  #     R0Sd7 <<- sd(R0OverInfectiousPeriod[which((R0OverInfectiousPeriod!=0) & (recovered == 1))])
  #   }
  #   
  # }
  return(FALSE)
}

sAge <- function(age){
  if (age <= 3){return(sAgeDist[1])}
  else if (age <= 20){return(sAgeDist[2])}
  else if (age <= 40){return(sAgeDist[3])}
  else if (age <= 60){return(sAgeDist[4])}
  else if (age <= 80){return(sAgeDist[5])}
  else if (age <= 100){return(sAgeDist[6])}
}

# find possible connections in 2d case next neighbor
findConn <- function(){
  counter <- 0
  
  for(j in c(1:(sqrt(N)-1))){
    for(i in c(1:(sqrt(N)-1))){
      
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,j]  
      possConn[counter,2] <<- lattice [i+1,j]
      #possConn[counter,3] <<- transProb(1,0)
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j] 
      possConn[counter,2] <<- lattice [i,j+1]
      #possConn[counter,3] <<- transProb(0,1)
    }
    counter <- counter + 1
    possConn[counter,1] <<- lattice[sqrt(N),j]  
    possConn[counter,2] <<- lattice [sqrt(N),j+1] 
    #possConn[counter,3] <<- transProb(0,1)
    
  }
  for (i in c(1:(sqrt(N)-1))){
    counter <- counter + 1
    possConn[counter,1] <<- lattice[i,sqrt(N)]  
    possConn[counter,2] <<- lattice [i+1,sqrt(N)]
    #possConn[counter,3] <<- transProb(1,0)
  }
  
  if(periodicBoundaries){
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[1,i]
      possConn[counter,2] <<- lattice[sqrt(N),i]
      #possConn[counter,3] <<- transProb(1,0)
    }
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,1]
      possConn[counter,2] <<- lattice[i,sqrt(N)]
      #possConn[counter,3] <<- transProb(0,1)
    }
  }
  if(nShort != 0){   
    while (TRUE) {
      if(nShort != 0){
        counter <- counter + 1
        condition <- TRUE
        while (condition == TRUE) {
          fromRow <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
          fromCol <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
          toRow   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
          toCol   <- sample(c(1:sqrt(N)), size = nShort, replace = TRUE)
          for (i in c(1:nShort)){
            if (!(fromRow[i] == toRow[i] && fromCol[i] == toCol[i])){
              condition <- FALSE
              break
            }
          }
        }
        for (i in c(0:(nShort-1))){
          possConn[counter + i,1] <<- lattice[fromRow[i+1],fromCol[i+1]]
          possConn[counter + i,2] <<- lattice[toRow[i+1],toCol[i+1]]
          #possConn[counter + i,3] <<- transProb((fromRow[i+1]-toRow[i+1]), (fromCol[i+1]-toCol[i+1]))
          
          # not left < right with shortcuts #!!
        }
        counter <- counter - 1
        if (!(TRUE %in% duplicated(possConn))){break}
      }
    }
  }
  
}


transProb <- function(x,y){ #this function determines the likelihood of transmission along the connection from x to y
  if(abs(x) > sqrt(N)/2){
    x <- sqrt(N) - abs(x)
  }
  if(abs(y) > sqrt(N)/2){
    y <- sqrt(N) - abs(y) 
  }
  dist <- sqrt(x**2 + y**2)
  prob <- 1/dist**bondOccPow*0.5
  return(prob) 
}


#cluster evaluation
#------
findRoot <- function(startIndex){
  #find root plus path compression
  root <- startIndex
  while (people[root] != root){
    root <- people[root]
  }
  while (people[startIndex] != root){
    parent <- people[startIndex]
    people[startIndex] <<- root # path compression:  -> always point to root
    startIndex <- parent
  }
  return(root)
}

addConnection <- function(from, to){ # does this work with shortcuts? I think yes, but not 100% sure
  #definetly get a problem when people can recover, as this should interrupt the tree
  #this cannot be taken care of in hindsight, as the trees do not contain this information anymore
  # maybe sacrifice the path compression?
  #find roots for both parts of conn (A,B) [path compression]
  fromRoot <- findRoot(from)#index of root node 
  toRoot <- findRoot(to)
  if (fromRoot == toRoot){
    return(0)
  }# if same: nothing
  else {
    if (weight[fromRoot] > weight[toRoot]){
      tmp <- fromRoot
      fromRoot <- toRoot
      toRoot <- tmp
    }
    # write connection "from root node" to "to root node"
    people[fromRoot] <<- toRoot 
    weight[toRoot] <<- weight[fromRoot] + weight[toRoot] # this line takes longer (50-200x) when checkLargestCluster = True. wat?
    #.subset(weight)[toRoot] <<- .subset(weight)[fromRoot] + .subset(weight)[fromRoot]
    weight[fromRoot] <<- 0
    #merges <<- merges +1
    
  }
  return(weight[toRoot])
}


