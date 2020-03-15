library(data.tree)

a <- Node$new(paste0(1,",",2))
b<- a$AddChild(paste0(2,",",2))
c<- b$AddChild(paste0(2,",",3))
d <- b$AddChild(paste0(3,",",3))
d$value <- 29
b$i <- 2
b$j <- 2
d$i <- 3
d$j <- 3
print(a)
e <- Climb(a, name = paste0(3,",",3))

print(e$value)
e$AddChild("newChild")
print(a)
startTime <- proc.time()
print(FindNode(a, paste0(3,",",3)))
endTime <- proc.Time()
runTime <- endTime-startTime
