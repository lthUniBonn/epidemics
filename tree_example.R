library(data.tree)

acme <- Node$new("Acme Inc.")
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

testNode <- Node$new("new node")
testNode$cost <- 2500
trees <- c(acme, testNode)
testNode <- Node$new("new node 2")
testNode$cost <- 3000


test <- trees[[1]]$AddChildNode(trees[[2]])
test2 <- trees[[1]]$AddChildNode(testNode)
check <- acme$Get("cost", filterFun = isLeaf)
print(check)
