library('plot.matrix')
startTime <- proc.time()
numberPeople = 30**2 #as the code is written right know this needs to be a square number
set.seed(4)
person <- array(0, dim = c(sqrt(numberPeople),sqrt(numberPeople)))
contagionProbability <- array(0, dim=c(numberPeople,3))
colnames(contagionProbability) <- c("proximityFactor","Immunity","Infected Time")
contagionProbability[,"Immunity"] <- 1 # 1 means not immune, 0 means fully immune, could in principle go for 0.5 as well
contagionProbability[,"Infected Time"] <- 0
ImmunisationQuote <- 0.5
for(i in c(1:numberPeople)){
  if(runif(1)<ImmunisationQuote){
    contagionProbability[i,"Immunity"] <- 0
  }
}
patientZero <- numeric(length = 1)
while(patientZero == 0){
  patientZero <- sample(1:numberPeople,1)
}
person[patientZero] <- 1
numberInfected <- 1




doTimeStep <- function(){#testet wie nah an den infizierten jede person ist und infiziert die neuen
  
  for(test in c(1:numberPeople)){
      contagionProbability[test,"proximityFactor"] <<- checkSurrounding(test=test) 
      if(person[test] == 1){
        contagionProbability[test,"Infected Time"] <<- contagionProbability[test,"Infected Time"] + 1
      }
  } 
  #return(contaminate(contagionProbability))
  contaminate(contagionProbability)
  numberInfected<<-sum(person)
  plot(person)
  return(person)
}

contaminate <- function(contagionProbability){
  #print(contagionProbability)
  for(i in c(1:numberPeople)){
    if(contagionProbability[i,"proximityFactor"]*contagionProbability[i,"Immunity"]> runif(1)){
      person[i] <<- 1
      #print("infected someone")
    }
    if(contagionProbability[i,"Infected Time"] == 2){ #should probably make this so only a chance to heal
      person[i] <<- 0
      contagionProbability[i,"Immunity"] <- 0
    }
  }
  #return(person)
}

checkSurrounding <- function(test){
  hereBePestilence <- which(person==TRUE)
  for(i in c(1:length(hereBePestilence))){
    if(test != hereBePestilence[i]){ #rules out that the infected person is tested for distance to itself
      distance <- calculateDistance(test=test,ill=hereBePestilence[i])
      contagionProbability[,"proximityFactor"] <- contagionProbability[,"proximityFactor"] + 1/distance**3 # 
    }
    
  }
  return(contagionProbability[test,"proximityFactor"])

}

calculateDistance <- function(test=1,ill=3){ #this assumes a square grid of people, could make it so that the opposing edges connect
  rowtest <- test %% sqrt(numberPeople)
  columntest <- test%/% sqrt(numberPeople)
  rowill <- ill %% sqrt(numberPeople)
  columnill <- ill %/% sqrt(numberPeople)
  distance <- sqrt((rowtest-rowill)**2+(columntest - columnill)**2)
  return(distance)
}
while(numberInfected !=0) {
  doTimeStep()
}
endTime <- proc.time()
runTime <- endTime-startTime
#person <- doTimeStep()


