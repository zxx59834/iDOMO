\name{combinationSET}
\alias{combinationSET}

\title{
Calculate Drug Combination Score for Individual Drug Pair
}
\description{
This function calculates drug combination score for a disease signature and 2 drug signatures using the super exact test. 
}
\usage{
combinationSET(disSet, drug1Set, drug2Set, allGenes, 
			lower.tail = FALSE, log.p = FALSE)
}

\arguments{
  \item{disSet}{
Disease signature, a list of length 2 with the up-regulated genes in the 1st cell and down-regulated genes in the 2nd cell.  
}
  \item{drug1Set}{
Drug signature A, a list of length 2 with the up-regulated genes in the 1st cell and down-regulated genes in the 2nd cell.  
}
  \item{drug2Set}{
Drug signature B, a list of length 2 with the up-regulated genes in the 1st cell and down-regulated genes in the 2nd cell.  
}
  \item{allGenes}{
The gene set in which the disease signature and drug signature were identified. 
}
  \item{lower.tail}{
Logical. Whether the p-values will be calculated for the low tail. Default is FALSE. 
}
  \item{log.p}{
Logical. Whether the p-values will be log transformed. Default is FALSE. 
}
}

\value{
A list of "geneNums", "pvals" and "Score". geneNums contains the overlapping gene numbers in each category. 
pvals contains the p-values of the overlapping genes in each category. 
Score is the final combination score considering the beneficial and harmful categories. 
}
\references{
Xianxiao Zhou et al. iDOMO: Identification of Drug Combinations via Multi-Set Operations for Treating Complex Diseases. 
Minghui Wang, Yongzhong Zhao, and Bin Zhang (2015). Efficient Test and Visualization of MultiSet Intersections. Scientific Reports 5: 16923.
}
\author{
Xianxiao Zhou <xianxiao.zhou@mssm.edu>
}

\seealso{
\code{\link{SET.Combn}}
}
\examples{

## Run the TNBC example for the top 2 drugs
data("TNBCsig")
data("Top20DrugSig")
data("allGenes")
combinationSET(TNBCsig, Top20DrugSig[[1]], Top20DrugSig[[2]],allGenes)

}


\keyword{ combination }
\keyword{ p-value }
