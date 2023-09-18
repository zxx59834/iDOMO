combinationSET <-
function(disSet, drug1Set, drug2Set, allGenes, 
							lower.tail=FALSE,log.p=FALSE) {

## Calculate beneficial score using revised Super Exact Test

n <- length(allGenes)
disSetUp <- disSet[[1]]
disSetDn <- disSet[[2]]
options(stringsAsFactors = FALSE)
disSetUp <- intersect(disSet[[1]],allGenes)
disSetDn <- intersect(disSet[[2]],allGenes)
disUpDn  <- intersect(disSetUp, disSetDn)

if (length(disUpDn) > 0) {
	warning(paste(length(disUpDn), "genes presenting in both up- AND down-regulated 
				disease signature will be removed from analysis. "))
	disSetUp <- setdiff(disSetUp, disUpDn)
	disSetDn <- setdiff(disSetDn, disUpDn)
}
## obtain union disease signature, drug1 signature and drug2 signature.
disSet  <- unique(c(disSetUp, disSetDn))
drug1SetUp <- intersect(drug1Set[[1]],allGenes)
drug1SetDn <- intersect(drug1Set[[2]],allGenes)
drug1Set   <- c(drug1SetUp, drug1SetDn)
drug2SetUp <- intersect(drug2Set[[1]],allGenes)
drug2SetDn <- intersect(drug2Set[[2]],allGenes)
drug2Set   <- c(drug2SetUp, drug2SetDn)

DSup   <- length(disSetUp)
DSdn   <- length(disSetDn)
DRS1up <- length(drug1SetUp)
DRS1dn <- length(drug1SetDn)
DRS2up <- length(drug2SetUp)
DRS2dn <- length(drug2SetDn)

x1 <- intersect(intersect(disSetDn, drug1SetUp), drug2SetUp)
x2 <- setdiff(intersect(disSetDn, drug1SetUp), drug2Set)
x3 <- setdiff(intersect(disSetDn, drug2SetUp), drug1Set)
# x4 x5 x6 are Neutral 
# x4 <- intersect(intersect(drug1SetUp,drug2SetDn),disSetDn)
# x5 <- intersect(intersect(drug2SetUp,drug1SetDn),disSetDn)
# x6 <- setdiff(disSetDn, union(drug1Set,drug2Set))
## 
x7 <- setdiff(intersect(disSetDn, drug2SetDn), drug1Set)
x8 <- setdiff(intersect(disSetDn, drug1SetDn), drug2Set)
x9 <- intersect(intersect(disSetDn, drug1SetDn),drug2SetDn)
x10 <- setdiff(intersect(drug1SetUp, drug2SetUp), disSet)
x11 <- setdiff(drug1SetUp, union(disSet, drug2Set))
x12 <- setdiff(drug2SetUp, union(disSet, drug1Set))
## x13 x14 x15 are Neutral
# x13 <- setdiff(intersect(drug1SetUp,drug2SetDn),disSet)
# x14 <- setdiff(intersect(drug2SetUp,drug1SetDn),disSet)
##
x16 <- setdiff(drug2SetDn, union(disSet, drug1Set))
x17 <- setdiff(drug1SetDn, union(disSet, drug2Set))
x18 <- setdiff(intersect(drug1SetDn, drug2SetDn), disSet)
x19 <- intersect(intersect(disSetUp, drug1SetUp), drug2SetUp)
x20 <- setdiff(intersect(disSetUp, drug1SetUp), drug2Set)
x21 <- setdiff(intersect(disSetUp, drug2SetUp), drug1Set)
## x22 x23 x24 are Neutral
# x22 <- intersect(intersect(drug1SetUp,drug2SetDn),disSetUp)
# x23 <- intersect(intersect(drug2SetUp,drug1SetDn),disSetUp)
# x24 <- setdiff(disSetUp, union(drug1Set,drug2Set))
##
x25 <- setdiff(intersect(disSetUp, drug2SetDn), drug1Set)
x26 <- setdiff(intersect(disSetUp, drug1SetDn), drug2Set)
x27 <- intersect(intersect(disSetUp, drug1SetDn), drug2SetDn)

## Calculate P-value for each number
p1 <- cpsets(length(x1),c(DSdn, DRS1up, DRS2up),n,lower.tail=lower.tail, log.p=log.p)
p2 <- cpdiff(length(x2),DSdn,DRS1up,DRS2up+DRS2dn,n,lower.tail=lower.tail, log.p=log.p)
p3 <- cpdiff(length(x3),DSdn,DRS2up,DRS1up+DRS1dn,n,lower.tail=lower.tail, log.p=log.p)

p7 <- cpdiff(length(x7),DSdn,DRS2dn,DRS1up+DRS1dn,n,lower.tail=lower.tail, log.p=log.p)
p8 <- cpdiff(length(x8),DSdn,DRS1dn,DRS2up+DRS2dn,n,lower.tail=lower.tail, log.p=log.p)
p9 <- cpsets(length(x9),c(DSdn, DRS1dn, DRS2dn),n,lower.tail=lower.tail, log.p=log.p)

p10 <- cpdiff(length(x10),DRS1up,DRS2up,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)
p11 <- cpone(length(x11),DRS1up,DRS2up+DRS2dn,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)
p12 <- cpone(length(x12),DRS2up,DRS1up+DRS1dn,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)

p16 <- cpone(length(x16),DRS2dn,DRS1up+DRS1dn,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)
p17 <- cpone(length(x17),DRS1dn,DRS2up+DRS2dn,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)
p18 <- cpdiff(length(x18),DRS1dn,DRS2dn,DSup+DSdn,n,lower.tail=lower.tail, log.p=log.p)

p19 <- cpsets(length(x19),c(DSup, DRS1up, DRS2up),n,lower.tail=lower.tail, log.p=log.p)
p20 <- cpdiff(length(x20),DSup,DRS1up,DRS2up+DRS2dn,n,lower.tail=lower.tail, log.p=log.p)
p21 <- cpdiff(length(x21),DSup,DRS2up,DRS1up+DRS1dn,n,lower.tail=lower.tail, log.p=log.p)

p25 <- cpdiff(length(x25),DSup,DRS2dn,DRS1up+DRS1dn,n,lower.tail=lower.tail, log.p=log.p)
p26 <- cpdiff(length(x26),DSup,DRS1dn,DRS2up+DRS2dn,n,lower.tail=lower.tail, log.p=log.p)
p27 <- cpsets(length(x27),c(DSup, DRS1dn, DRS2dn),n,lower.tail=lower.tail, log.p=log.p)

Bpvalue <- c(p1,p2,p3,p25,p26,p27)
Hpvalue <- c(p7,p8,p9,p11,p12,p16,p17,p19,p20,p21)## p10,p18,

Bpvalue[Bpvalue == 0] <- 1e-323
Hpvalue[Hpvalue == 0] <- 1e-323
Bscore <- sum(-log10(Bpvalue))
Hscore <- sum(-log10(Hpvalue))

SS <- Bscore - Hscore ## final score
## Output
# return(c(x1,x2,x3,x7,x8,x9,x10,x11,x12,x16,x17,x18,x19,x20,x21,x25,x26,x27,
	# p1,p2,p3,p7,p8,p9,p10,p11,p12,p16,p17,p18,p19,p20,p21,p25,p26,p27,SS)) ## output vector
return(list(geneNums = c(length(x1),length(x2),length(x3),length(x7),length(x8), ## ,length(x4),length(x5),length(x6)
	length(x9),length(x10),length(x11),length(x12),length(x16),length(x17),length(x18), ## ,length(x13),length(x14)
	length(x19),length(x20),length(x21),length(x25),length(x26),length(x27)), ## ,length(x22),length(x23),length(x24)
	pvals = c(p1,p2,p3,p7,p8,p9,p10,p11,p12,p16,p17,p18,p19,p20,p21,p25,p26,p27), Score = SS))
}
