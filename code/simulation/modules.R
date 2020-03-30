library('plot.matrix')

#----------- run functions (sim + write)
figCount <- 0

simulationRun <- function(statRun){
  sDistribution <<- initialsDistribution 
  #set Infection status
  infected <<- logical(length = N)
  infectionTime <<- numeric(length=N)
  recovered <<- logical(length=N)
  
  
  # infect Patient 0
  patientZero()
  
  
  x <<- 0
  count <- 0
  while (TRUE) {
    infPeople <- which(infected==TRUE)
    noMoreInfected <- timesteps()
    # returns NA if no more infected
    if(noMoreInfected){ break}
    # accumulated plot to png in figpath of infected | currently infected in red
    if((x %% plotEvery == 0) && (plotAccumulated == TRUE)){
      visibleLattice <- array(0, dim= c(sqrt(N),sqrt(N)))
      visibleLattice[which(sDistribution  != initialsDistribution)] <- 1
      plot(which(visibleLattice==1, arr.ind = TRUE)[,1], which(visibleLattice==1, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '.', pty = "s", xlab = paste(c("T:", x), sep = " ", collapse = ""), ylab ="")
      par(new=TRUE)
      visibleLattice[which(infected == TRUE)] <- 2
      plot(which(visibleLattice==2, arr.ind = TRUE)[,1], which(visibleLattice==2, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '*',col = 'red', pty = "s", xlab = "", ylab = "")
      figCount <- figCount +1
      
    }
    
    # evaluation of disease parameters
    if ((x %% clusterEvery == 0) && (checkCluster == TRUE)){
      count <- count + 1
      people <<- c(1:N)
      weight <<- numeric(N) 
      # now identifiy clusters use Newman Algo with infConn as used bonds
      infPeople <- which(infectionTime != 0)
      weight[infPeople] <<- 1
      infConn <- possConn[which((possConn[,1] %in% infPeople) & (possConn[,2] %in% infPeople)),c(1,2)]
      if (!is.array(infConn)){infConn <- array(data = infConn, dim = c(1,2))}
      if(nrow(infConn) != 0){
        for(i in c(1:nrow(infConn))){
          
          addConnection(infConn[i,1],infConn[i,2])
        }
      }
      #write to evaluationArray[x,statRun, observable]#!!
      evalArray[count,statRun, 'x'] <<- x
      evalArray[count,statRun, 'maxWeight'] <<- weight[which.max(weight)]
      evalArray[count,statRun, 'numberCluster'] <<- length(which(weight!=0))
      evalArray[count,statRun, 'numberInfected'] <<- sum(weight)
      evalArray[count,statRun, 'largeOverTotal'] <<- weight[which.max(weight)]/sum(weight)
      evalArray[count,statRun, 'accInfections'] <<- (length(which(sDistribution != initialsDistribution))+length(infPeople))/N
      
    }
    x <<- x + 1
  }
}

writeEval <- function(writeThis, params){
  #write observation array to files
  for(obs in obsNames[-1]){
    writeThisObs <- array(NA, dim= c(20000, nStatRun+1))
    writeThisObs[,1] <- seq(0,20000*clusterEvery-1,clusterEvery)
    writeThisObs[,c(2:(nStatRun+1))] <- writeThis[,c(1:(nStatRun)),obs] 
    
    isNa <- is.na(writeThisObs[,2])
    for(s in c(2:(nStatRun))){
      isNa <- isNa & is.na(writeThisObs[,s+1])
    }
    
    writeThisObs <- writeThisObs[c(1:(which(isNa)[1]-1)), c(1:(nStatRun+1))]
    
    if(is.vector(writeThisObs)){
      writeThisObs <-t(writeThisObs)
    }
    write.table(x = writeThisObs, file =paste(c(path, obs,"_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
  }
}


#--------------- simulation parts

patientZero <- function(){
  zeroId <- sample(x = which(sDistribution != 0),size = 1) #choose one of the susceptibles to be patient 0
  infected[zeroId] <<- TRUE
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
  
  # split infConn for significantly speedinf up the simulation
  possiblyInf1 <- infConn[,1]
  possiblyInf2 <- infConn[,2]
  
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
  
  primaryInfected <- length(infPeople)
  
  infected[newlyInf1] <<- TRUE
  infected[newlyInf2] <<- TRUE
  totalInfected <- length(which(infected == TRUE)) #newly infected (including recovering)
  infected[recPeople] <<- FALSE
  
  #increase infection time 
  infectionTime[infected] <<- infectionTime[infected] + 1
  
  
  recovered[recPeople] <<- 1
  return(FALSE)
}

sAge <- function(age){
  #point age to corresponding age distribution
  if (age <= 3){return(sAgeDist[1])}
  else if (age <= 20){return(sAgeDist[2])}
  else if (age <= 40){return(sAgeDist[3])}
  else if (age <= 60){return(sAgeDist[4])}
  else if (age <= 80){return(sAgeDist[5])}
  else if (age <= 100){return(sAgeDist[6])}
}

#--------------- # find possible connections in 2d case next neighbor
findConn <- function(){
  counter <- 0
  
  for(j in c(1:(sqrt(N)-1))){
    for(i in c(1:(sqrt(N)-1))){
      
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,j]  
      possConn[counter,2] <<- lattice [i+1,j]
      counter <- counter +1 
      possConn[counter,1] <<- lattice[i,j] 
      possConn[counter,2] <<- lattice [i,j+1]
    }
    counter <- counter + 1
    possConn[counter,1] <<- lattice[sqrt(N),j]  
    possConn[counter,2] <<- lattice [sqrt(N),j+1] 
    
  }
  for (i in c(1:(sqrt(N)-1))){
    counter <- counter + 1
    possConn[counter,1] <<- lattice[i,sqrt(N)]  
    possConn[counter,2] <<- lattice [i+1,sqrt(N)]
  }
  
  if(periodicBoundaries){
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[1,i]
      possConn[counter,2] <<- lattice[sqrt(N),i]
    }
    for(i in c(1:sqrt(N))){
      counter <- counter + 1
      possConn[counter,1] <<- lattice[i,1]
      possConn[counter,2] <<- lattice[i,sqrt(N)]
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
        }
        counter <- counter - 1
        if (!(TRUE %in% duplicated(possConn))){break}
      }
    }
  }
  
}


#--------------- cluster evaluation
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

addConnection <- function(from, to){ 
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
    weight[toRoot] <<- weight[fromRoot] + weight[toRoot] 
    weight[fromRoot] <<- 0
  }
  return(weight[toRoot])
}