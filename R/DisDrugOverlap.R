DisDrugOverlap <-
function(disSet, drugSets, RefGenes, 
							lower.tail=FALSE,log.p=FALSE) {
n <- length(RefGenes)
disSetUp <- disSet[[1]]
disSetDn <- disSet[[2]]
options(stringsAsFactors = FALSE)
disSetUp <- intersect(disSet[[1]],RefGenes)
disSetDn <- intersect(disSet[[2]],RefGenes)
disUpDn  <- intersect(disSetUp, disSetDn)

if (length(disUpDn) > 0) {
	warning(paste(length(disUpDn), "genes presenting in both up- AND down-regulated 
				disease signature will be removed from analysis. "))
	disSetUp <- setdiff(disSetUp, disUpDn)
	disSetDn <- setdiff(disSetDn, disUpDn)
}
## obtain union disease signature, drug1 signature and drug2 signature.
disSet  <- unique(c(disSetUp, disSetDn))

l <- length(drugSets)
Ol <- matrix(0,l,4)
Pv <- matrix(1,l,4)
SS <- rep(0,l)
for (i in 1:l) {
	drugSetUp <- intersect(drugSets[[i]][[1]],RefGenes)
	drugSetDn <- intersect(drugSets[[i]][[2]],RefGenes)

	DSup   <- length(disSetUp)
	DSdn   <- length(disSetDn)
	DRSup <- length(drugSetUp)
	DRSdn <- length(drugSetDn)

	Ol[i,1] <- length(intersect(disSetDn, drugSetUp))
	Ol[i,2] <- length(intersect(disSetUp, drugSetUp))
	Ol[i,3] <- length(intersect(disSetUp, drugSetDn))
	# x4 x5 x6 are Neutral 
	Ol[i,4] <- length(intersect(disSetDn,drugSetDn))
	# Ol[i,5] <- length(setdiff(drugSetUp,disSet))
	# Ol[i,6] <- length(setdiff(drugSetDn, disSet))

	Pv[i,1] <- phyper(Ol[i,1]-1,DSdn,n-DSdn,DRSup,lower.tail=lower.tail, log.p=log.p)
	Pv[i,2] <- phyper(Ol[i,2]-1,DSup,n-DSup,DRSup,lower.tail=lower.tail, log.p=log.p)
	Pv[i,3] <- phyper(Ol[i,3]-1,DSup,n-DSup,DRSdn,lower.tail=lower.tail, log.p=log.p)
	Pv[i,4] <- phyper(Ol[i,4]-1,DSdn,n-DSdn,DRSdn,lower.tail=lower.tail, log.p=log.p)

	# Pv[i,5] <- cpone(Ol[i,5],DRSdn,DSdn,DSup,n,lower.tail=lower.tail, log.p=log.p)
	# Pv[i,6] <- cpone(Ol[i,6],DRSup,DSdn,DSup,n,lower.tail=lower.tail, log.p=log.p)

	Bpvalue <- c(Pv[i,1],Pv[i,3])
	Hpvalue <- c(Pv[i,2],Pv[i,4]) ## ,p6,p8

	Bpvalue[Bpvalue == 0] <- 1e-323
	Hpvalue[Hpvalue == 0] <- 1e-323
	Bscore <- sum(-log10(Bpvalue))
	Hscore <- sum(-log10(Hpvalue))

	SS[i] <- Bscore - Hscore ## final score
}
Output <- data.frame(Drug=names(drugSets), Ol, Pv, SS)
colnames(Output) <- c("Drug","Ol1","Ol2","Ol3","Ol4","p1","p2","p3","p4","Score")
return(Output)
}
