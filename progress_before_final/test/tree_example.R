library(data.tree)

acme <- Node$new(5)
accounting <- acme$AddChild("Accounting")
software <- accounting$AddChild("New Software")
standards <- accounting$AddChild("New Accounting Standards")
research <- acme$AddChild("Research")
newProductLine <- research$AddChild("New Product Line")
newLabs <- research$AddChild("New Labs")
it <- acme$AddChild("IT")
outsource <- it$AddChild("Outsource")
agile <- it$AddChild("Go agile")
goToR <- it$AddChild("Switch to R")
outsource$cost <- 25
it$cost <- 23

print(Climb(acme, name = "5", cost = 23))
testList <- list()
#testList[[2]] <- Node$new("new node")
#testNode$cost <- 2500
#trees <- c(acme, testNode)
#testNode <- Node$new("new node 2")
#trees[[1]]$cost <- 3000
print(acme)
print(acme$children[[1]]$children[[2]]$name)
print(acme$Climb(position = c(1, 2)))
# test <- trees[[3]]$AddChildNode(trees[[2]])
# test2 <- trees[[3]]$AddChildNode(testNode)
# check <- acme$Get("cost")
# print(check)
# 
clusterCounter <- 2
i <- 1
j <- 1
treeList <- list()
treeList[[clusterCounter]] <- Node$new(toString(clusterCounter))
treeList[[clusterCounter]]$i <- i
treeList[[clusterCounter]]$j <- j
bla <- treeList[[clusterCounter]]$AddChild("node")
bla$example <- 2
bla$j <- 4
fuu <- bla$AddChild("fuu")
fuu$example <- 2
fuu$j <- 4
print(treeList[[clusterCounter]])
print(Climb(treeList[[clusterCounter]], name = "node", example=2))
# print(treeList[[clusterCounter]]$i)
# print(treeList[[clusterCounter]]$j)
