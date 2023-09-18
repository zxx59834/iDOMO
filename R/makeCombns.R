makeCombns <-
function(rankedDrugs, topn) {
## make drug combinations for the rankedDrugs, with pairwise 
## combination among topn as well as all combinations between
## topn drugs and the rest drugs
topnDrugs <- rankedDrugs[1:topn]
otherDrugs <- rankedDrugs[(topn+1):length(rankedDrugs)]
topnPairwise <- combn(topnDrugs, 2)
topnPairwise <- t(topnPairwise)

otherCombns <- expand.grid(topnDrugs, otherDrugs)
colnames(otherCombns) <- c("drug1", "drug2")
colnames(topnPairwise) <- c("drug1", "drug2")
rbind(topnPairwise, otherCombns)
}
