bootstrap <- function(vec, nSample=50){
  means <- numeric(length=nSample)
  vec <- vec[which(is.na(vec)==FALSE)]
  for(i in c(1:nSample)){
    means[i] <- mean(sample(vec, length(vec), replace=TRUE))
  }
  err <- sd(means) 
  return(err)
}

meanPlot <- function(name,params,df, compare = FALSE, name2, params2, df2){
  thisObsMean <- rowMeans(df[,c(2:ncol(df))], na.rm = TRUE)
  thisObsErr <- apply(X = df[,c(2:ncol(df))],MARGIN = 1, FUN = bootstrap)
  
  if(compare==TRUE){
    thisObsMean2 <- rowMeans(df2[,c(2:ncol(df2))], na.rm = TRUE)
    thisObsErr2 <- apply(X = df2[,c(2:ncol(df2))],MARGIN = 1, FUN = bootstrap)
    xlim <- c(0,max(df[which.max(df[,1]),1], df2[which.max(df2[,1]),1]))
    ylim <- c(0,max(max(thisObsMean, na.rm = T)+max(thisObsErr, na.rm = T),(max(thisObsMean2, na.rm = T)+max(thisObsErr2, na.rm = T))))
    
  } else {
    xlim <- c(0,df[which.max(df[,1]),1]) 
    ylim <- c(0,max(thisObsMean)+max(thisObsErr, na.rm = T))
  }
  
  
  if(compare == TRUE){
    plot(df2[,1], thisObsMean2, xlab = "T",ylab = name, col = 'red', xlim = xlim, ylim = ylim)
    arrows(df2[,1], thisObsMean2-thisObsErr2, df2[,1], thisObsMean2+thisObsErr2, length=0.05, angle=90, code=3,col = 'red')
    par(new = TRUE)  
  }
  
  plot(df[,1], thisObsMean, xlab = "T",ylab = name, col = 'black', xlim = xlim, ylim = ylim)
  arrows(df[,1], thisObsMean-thisObsErr, df[,1], thisObsMean+thisObsErr, length=0.05, angle=90, code=3)
  mtext(text = paste("sqrt(N):", params[1],"  nShort:", params[2], "  immunity:", params[3],"  recoveryTime:", params[4],"+-", params[5],"  suscDist:", params[6],"  suscFactor:", params[7], sep = " "),side = 3)
  if(compare == TRUE){
    mtext(text = paste("sqrt(N):", params2[1],"  nShort:", params2[2], "  immunity:", params2[3],"  recoveryTime:", params2[4],"+-", params2[5],"  suscDist:", params2[6],"  suscFactor:", params2[7], sep = " "),side = 3,col = 'red',line = 1)
  }
}
