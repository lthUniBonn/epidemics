N <- 400**2
nShort <- N/100
immunity <- 0
avgRecoveryTime <- 4
sdRecoveryTime <- 1
i <- 1
sAgeDist0 <- c(1, 0.7, 0.5, 0.7, 0.8, 1)
sAgeDist <- sAgeDist0
sDistFactor <- 3
sChoice <- c(F,F,F,T)
sChoiceNames <- c("sBool", "sFixed", "sNot", "sReal")
sBool <- sChoice[1] # if True the susceptibility is 1 or 0 // other poss like gaussian with age etc
sFixed <- sChoice[2]
sNot <- sChoice[3]
sReal <- sChoice[4]

#params <- paste(c(sqrt(N), nShort, immunity, avgRecoveryTime, sdRecoveryTime, i, sDistFactor, sChoiceNames[sChoice]), sep="", collapse="_") 


meanPlot <- function(name,params,df, compare = FALSE, name2, params2, df2){
  thisObsMean <- rowMeans(df[,c(2:ncol(df))], na.rm = TRUE)
  thisObsSd <- apply(X = df[,c(2:ncol(df))],MARGIN = 1, FUN = sd, na.rm =TRUE)
  
  if(compare==TRUE){
    thisObsMean2 <- rowMeans(df2[,c(2:ncol(df2))], na.rm = TRUE)
    thisObsSd2 <- apply(X = df2[,c(2:ncol(df2))],MARGIN = 1, FUN = sd, na.rm =TRUE)
    xlim <- c(0,max(df[which.max(df[,1]),1], df2[which.max(df2[,1]),1]))
    print(max(thisObsSd))
    print(max(thisObsSd2))
    ylim <- c(0,max(max(thisObsMean, na.rm = T)+max(thisObsSd, na.rm = T),(max(thisObsMean2, na.rm = T)+max(thisObsSd2, na.rm = T))))
    
  } else {
    xlim <- c(0,df[which.max(df[,1]),1]) 
    ylim <- c(0,max(thisObsMean))
  }
  
  
  if(compare == TRUE){
    plot(df2[,1], thisObsMean2, xlab = "T",ylab = name, col = 'red', xlim = xlim, ylim = ylim)
    arrows(df2[,1], thisObsMean2-thisObsSd2, df2[,1], thisObsMean2+thisObsSd2, length=0.05, angle=90, code=3,col = 'red')
    par(new = TRUE)  
  }
  
  plot(df[,1], thisObsMean, xlab = "T",ylab = name, col = 'black',, xlim = xlim, ylim = ylim)
  arrows(df[,1], thisObsMean-thisObsSd, df[,1], thisObsMean+thisObsSd, length=0.05, angle=90, code=3)
  mtext(text = paste("N:", params[1]**2,"  nShort:", params[2], "  immunity:", params[3],"  recoveryTime:", params[4],"+-", params[5],"  suscDist:", params[6],"  suscFactor:", params[7], sep = " "),side = 3)
  if(compare == TRUE){
    mtext(text = paste("N:", params2[1]**2,"  nShort:", params2[2], "  immunity:", params2[3],"  recoveryTime:", params2[4],"+-", params2[5],"  suscDist:", params2[6],"  suscFactor:", params2[7], sep = " "),side = 3,col = 'red',line = 1)
  }
}
