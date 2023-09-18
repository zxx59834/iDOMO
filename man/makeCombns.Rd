\name{makeCombns}
\alias{makeCombns}
\title{
Make pair-wise combination for a drug list. 
}
\description{
For the top n drug list, the function makes pair-wise combination. And for the rest drugs, it expand grid combinations with the top n drugs. 
}
\usage{
makeCombns(rankedDrugs, topn)
}

\arguments{
  \item{rankedDrugs}{
A string vector with ranked drugs. For example, drug list from the EMUDRA prediction. 
}
  \item{topn}{
Positive integer, should be no less than 2. Top n drugs that will be used to make pair-wise combinations. 
}
}

\value{
A list of length 2. 
  \item{topnPairwise}{The pair-wise combinations from the top n drugs. }
  \item{otherCombns}{The combinations between the top n drugs and the rest drugs. }
}
\references{
Xianxiao Zhou et al. iDOMO: Identification of Drug Combinations via Multi-Set Operations for Treating Complex Diseases. 
}
\author{
Xianxiao Zhou <xianxiao.zhou@mssm.edu>
}

\examples{

makeCombns(letters, topn = 10)

}

\keyword{combination}

