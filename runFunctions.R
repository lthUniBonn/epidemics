simulationRun <- function(statRun){
  sDistribution <<- initialsDistribution 
  # R0OverInfectiousPeriod <<- numeric(length=N)
  # R0Mean <<- numeric()
  # R0Sd <<- numeric()
  # R0Mean2 <<- numeric()
  # R0Sd2 <<- numeric()
  # R0Mean3 <<- numeric()
  # R0Sd3 <<- numeric()
  # R0Mean4<<- numeric()
  # R0Sd4 <<- numeric()
  # R0Mean7 <<- numeric()
  # R0Sd7 <<- numeric()
  #set Infection status
  infected <<- logical(length = N)
  infectionTime <<- numeric(length=N)
  recovered <<- logical(length=N)
  
  
  # infect Patient 0
  patientZero()
  
  
  x <<- 0
  count <- 0
  while (TRUE) {
    noMoreInfected <- timesteps()
    # returns NA if no more infected
    if(noMoreInfected){ break}
    # accumulated plot of infected | currently infected in red
    if((x %% plotEvery == 0) && (plotAccumulated == TRUE)){
      
      png(paste(c(figPath, "/", figCount, ".png"), sep = "", collapse = ""), width = 500, height = 500)
      visibleLattice <- array(0, dim= c(sqrt(N),sqrt(N)))
      visibleLattice[which(sDistribution  != initialsDistribution)] <- 1
      plot(which(visibleLattice==1, arr.ind = TRUE)[,1], which(visibleLattice==1, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '.', pty = "s", xlab = paste(c("T:", x), sep = " ", collapse = ""), ylab ="")
      par(new=TRUE)
      visibleLattice[which(infected == TRUE)] <- 2
      plot(which(visibleLattice==2, arr.ind = TRUE)[,1], which(visibleLattice==2, arr.ind = TRUE)[,2],xlim = c(0,sqrt(N)), ylim = c(0,sqrt(N)),type ="p", pch = '*',col = 'red', pty = "s", xlab = "", ylab = "")
      dev.off()
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
    # if(anyNA(evalR0[statRun,])){
    #   if (sum(recovered)>=checkR0Here){
    #     evalR0[statRun, 1] <<- R0Mean
    #     evalR0[statRun, 2] <<- R0Sd
    #   }
    # }
    # if(anyNA(evalR02[statRun,])){
    #   if (sum(recovered)>=2){
    #     evalR02[statRun, 1] <<- R0Mean2
    #     evalR02[statRun, 2] <<- R0Sd2
    #   }
    # }
    # if(anyNA(evalR03[statRun,])){
    #   if (sum(recovered)>=3){
    #     evalR03[statRun, 1] <<- R0Mean3
    #     evalR03[statRun, 2] <<- R0Sd3
    #   }
    # }
    # if(anyNA(evalR04[statRun,])){
    #   if (sum(recovered)>=4){
    #     evalR04[statRun, 1] <<- R0Mean4
    #     evalR04[statRun, 2] <<- R0Sd4
    #   }
    # }
    # if(anyNA(evalR07[statRun,])){
    #   if (sum(recovered)>=7){
    #     evalR07[statRun, 1] <<- R0Mean7
    #     evalR07[statRun, 2] <<- R0Sd7
    #   }
    # }
    x <<- x + 1
  }
}

writeEval <- function(writeThis, params){#writeThisR0, params){
  
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
  
  
  #write.table(x = writeThisR0, file = paste(c(path, "R0Mean_", params, ".txt"),sep="", collapse=""), append = FALSE,sep = "\t",row.names = FALSE, col.names = FALSE)
}
