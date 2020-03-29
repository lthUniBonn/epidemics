parNames <- c('sqrt(N)',  'number of shortcuts', 'immunisation quota', 'avgerage recovery time', 'recovery time deviation', 'type of age distribution', 'social distancing factor', 'choice of susceptibility distribution')



bootstrap <- function(vec, nSample=1000, FUN=mean, qq = FALSE){
  popMean <- mean(vec, na.rm = T)
  means <- numeric(length=nSample)
  vec <- vec[which(is.na(vec)==FALSE)]
  for(i in c(1:nSample)){
    means[i] <- FUN(sample(vec, length(vec), replace=TRUE))
  }
  err <- sd(means) 
  if(qq == TRUE) {
    bias <- mean(means, na.rm = T)-popMean
    print(bias)
    qqnorm(y = means)
    qqline(y = means, distribution = qnorm)
    mtext(text = paste(c("sd: ", err, "  Bias: ", bias), sep = "", collapse = ""))
  }
  return(err)
}

outbreakProbability <- function(lastVal){
  return(length(which(lastVal>epidemicThreshold))/length(lastVal))
}


outbreakMeasure <- function(dfList, prob=T){
  p <- numeric(length(dfList))
  pErr <- numeric(length(dfList))
  lastValMean <- numeric(length(dfList))
  lastValErr <- numeric(length(dfList))
  for (dfIdx in c(1:length(dfList))){
    df <- dfList[[dfIdx]]
    x <- dfIdx
    lastVal <- apply(X = df[,c(2:ncol(df))],MARGIN = 2, FUN = max, na.rm =TRUE)
    lastValMean[x] <- mean(lastVal)
    lastValErr[x] <- bootstrap(lastVal)
    p[x] <- outbreakProbability(lastVal)
    pErr[x] <- bootstrap(lastVal, FUN = outbreakProbability)
  } 
  if(prob){
    return(list(p,pErr))
  } else { return(list(lastValMean, lastValErr))}
}


findObsvsParams <- function(obs='numberInfected', parIdx=3, params, checkSpecific=c()){
  #all files
  fileNames <- array(unlist(strsplit(list.files(path)[], "_")), dim=c(9,length(list.files(path))))
  paramList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files(path))))
  # find all files with obs 
  parList <- array(as.numeric(fileNames[c(2:8),]),dim=c(7,length(list.files(path))))
  obsList <- fileNames[1,]
  parList <- paramList[,which(obsList == obs)]
  if (length(checkSpecific) == 0){
    notParIdxVec <- c(1:7)
    notParIdxVec <- notParIdxVec[-which(notParIdxVec == parIdx)]
    
    #select param config, all but one fixed
    fixedParListCheck <- rep(TRUE, ncol(parList))
    for (idx in notParIdxVec){
      fixedParListCheck <- fixedParListCheck & (parList[idx,]==params[idx])
    }
    parList <- array(parList[,which(fixedParListCheck==TRUE)], dim = c(7, length(which(fixedParListCheck==TRUE))))
  } else {
    notParIdxVec <- c(1:7)
    notParIdxVec <- notParIdxVec[-which(notParIdxVec == parIdx)]
    #select param config, all but one fixed
    fixedParListCheck <- rep(TRUE, ncol(parList))
    for (idx in notParIdxVec){
      fixedParListCheck <- fixedParListCheck & (parList[idx,]==params[idx])
    }
    parList <- array(parList[,which(fixedParListCheck==TRUE)], dim = c(7, length(which(fixedParListCheck==TRUE))))
    specCheck <- rep(FALSE, ncol(parList))
    
    for (specIdx in c(1:length(checkSpecific))){
      specCheck[which(parList[parIdx,] == checkSpecific[specIdx])] <- TRUE
    }
    #fixedParListCheck <- fixedParListCheck & specCheck
    parList <- array(parList[,which(specCheck==TRUE)], dim = c(7, length(which(specCheck==TRUE))))
  }
  #read files into list of data frames 
  dfList <- list()
  for(idx in c(1:ncol(parList))){
    readParams <- paste(c(parList[1,idx], parList[2,idx], parList[3,idx], parList[4,idx], parList[5,idx], parList[6,idx], parList[7,idx], sChoice), sep="", collapse="_") 
    df <- read.table(file = paste(c(path,"/", obs, "_", readParams, ".txt"),sep="", collapse=""))
    dfList[[idx]] <- df
  }

  return(list(parList, dfList))
}



meanPlot <- function(name,params,df, compare = FALSE, name2="", params2=0, df2=0){
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
    plot(df2[,1], thisObsMean2, xlab = "",ylab = "", col = 'red', xlim = xlim, ylim = ylim, pch = 19, cex = 1, xaxt = 'n', yaxt = 'n')
    arrows(df2[,1], thisObsMean2-thisObsErr2, df2[,1], thisObsMean2+thisObsErr2, length=0.05, angle=90, code=3,col = 'red')
    
    par(new = TRUE)  
  }
  
  plot(df[,1], thisObsMean, xlab = "",ylab = "", col = 'black', xlim = xlim, ylim =ylim, pch = 19, cex = 1,xaxt = 'n', yaxt = 'n')
  title(ylab = name, line = 1.7, cex.lab = 1.5)
  title(xlab = "t", line = 1.7, cex.lab = 1.5)
  axis(2, mgp=c(3, .5, 0))
  axis(1, mgp=c(3, .5, 0))
  arrows(df[,1], thisObsMean-thisObsErr, df[,1], thisObsMean+thisObsErr, length=0.05, angle=90, code=3)
  mtext(text = paste( "r:", params2[3], sep = " "),side = 3, padj = 3.5, adj = 0.9, col = 'red', cex = 1.5)
  mtext(text = paste( "r:", params[3], sep = " "),side = 3, padj = 2, adj = 0.9, cex = 1.5)
  #mtext(text = paste("sqrt(N):", params[1],"  nShort:", params[2], "  immunity:", params[3],"  recoveryTime:", params[4],"+-", params[5],"  suscDist:", params[6],"  suscFactor:", params[7], sep = " "),side = 3)
   if(compare == TRUE){
     #mtext(text = paste("sqrt(N):", params2[1],"  nShort:", params2[2], "  immunity:", params2[3],"  recoveryTime:", params2[4],"+-", params2[5],"  suscDist:", params2[6],"  suscFactor:", params2[7], sep = " "),side = 3,col = 'red',line = 1)
   }
}
