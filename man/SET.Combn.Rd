\name{SET.Combn}
\alias{SET.Combn}

\title{
Calculate drug combination scores based on the signatures
}
\description{
Calculate the drug combination scores using the Super Exact Test based on disease signature gene sets and drug signature gene sets of drug pairs. 
}
\usage{
SET.Combn(disSet, DrugSets, allGenes, orderedOut = TRUE)
}

\arguments{
  \item{disSet}{
Disease signature. It should be a list with up-regulated genes in the first item and down-regulated genes in the 2nd item.
}
  \item{DrugSets}{
A list of drug signatures. It should be a list that each cell contains up-regulated genes in the first item and down-regulated genes in the 2nd item.
}
  \item{allGenes}{
The gene set in which the disease signature and drug signature were identified. 
}
  \item{orderedOut}{
Whether the output will be ordered by the combination score. Default is TRUE. 
}
}

\value{
\item{output}{A data frame of ranked drug combination list with final SET scores as well as beneficial scores and harmful scores. }

}
\references{
Xianxiao Zhou et al. iDOMO: Identification of Drug Combinations via Multi-Set Operations for Treating Complex Diseases. 

}
\author{
Xianxiao Zhou <xianxiao.zhou@mssm.edu>
}

\seealso{
\code{\link{combinationSET}}
}
\examples{

data("TNBCsig")
data("Top20DrugSig")
data("allGenes")
SET.Combn(TNBCsig,Top20DrugSig[1:3],allGenes, orderedOut = TRUE)

}


\keyword{ combination }
\keyword{ pair-wise }

