\name{DisDrugOverlap}
\alias{DisDrugOverlap}

\title{
Overlapping of Disease Signature and Individual Drug Signatures
}
\description{
This function calculates the scores of a disease signature and a list of individual drug signatures. 
}
\usage{
DisDrugOverlap(disSet, drugSets, RefGenes, 
				lower.tail = FALSE, log.p = FALSE)
}

\arguments{
  \item{disSet}{
Disease signature, containing 2 cells with up- and down-regulated genes in the diease signature. 
}
  \item{drugSets}{
A list of drug signatures. Each cell contains 2 cells with up- and down-regulated genes induced by the drug. 
}
  \item{RefGenes}{
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
\item{ Output }{ A data.frame of overlapping gene number and p-values as well as scores of the input drug signatures. }
}
\references{
Xianxiao Zhou et al. iDOMO: Identification of Drug Combinations via Multi-Set Operations for Treating Complex Diseases. 
}
\author{
Xianxiao Zhou <xianxiao.zhou@mssm.edu>
}

\examples{

data("TNBCsig")
data("Top20DrugSig")
data("allGenes")
DisDrugOverlap(TNBCsig,Top20DrugSig,allGenes, lower.tail = FALSE, log.p = FALSE)

}


