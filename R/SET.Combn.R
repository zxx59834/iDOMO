SET.Combn <-
function(disSet, DrugSets, allGenes, orderedOut = TRUE) {
## To calculate p-value based drug combination score using SuperExactTest
# (disSet) Disease signature should be a list with up-regulated genes and down-regulated genes
# drugCombns, is a two-column drug pair matrix
# DrugSets, a list, drug signatures. DrugSets[[i]][[1]] is a up-regulated gene list of drug i, while 
# DrugSets[[i]][[2]] is the down-regulated gene list for the same drug i. The name of DrugSets[[i]] 
# should be drug i. 
# n is total number of genes in gene expression profiles. 

l <- length(DrugSets)
# n <- length(allSymbols)
drugCombns <- t(combn(names(DrugSets),2))
m <- dim(drugCombns)[1]
geneNums <- matrix(0, m, 18)
pvals    <- matrix(1, m, 18)
Score    <- rep(0, m)


k <- 1
for (i in 1:(l-1)) {
	for (j in (i+1):l) {
		NumPscore <- combinationSET(disSet, DrugSets[[i]], DrugSets[[j]], allGenes)
	
		geneNums[k, ] <- NumPscore[[1]]
		pvals[k, ]    <- NumPscore[[2]]
		Score[k]      <- NumPscore[[3]]
		k <- k+1
	}
}
colnames(geneNums) <- paste("n", c(1:3, 7:12, 16:21, 25:27), sep = "") ## 1:14, 16:27
colnames(pvals) <- paste("p", c(1:3, 7:12, 16:21, 25:27), sep = "")
output <- data.frame(drugCombns, geneNums, pvals, Score)
if (orderedOut) {
	return(output[order(-output$Score), ])
} else {
	return(output)
}
}
