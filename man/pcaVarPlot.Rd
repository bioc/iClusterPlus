\name{pcaVarPlot}
\alias{pcaVarPlot}
\title{A function to make PCA variance plot of the combined matrix of
  multi-omics data matrices.}
\description{
  Multi-omics data matrices are column-combined and then the PCA
  variance plot of the combined matrix is made. 
}
\usage{
pcaVarPlot(xList,K=10)
}
\arguments{
  \item{xList}{A list object containing m data matrices representing m
    multi-omics data measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{K}{The number of principle components.}
}
\seealso{
  \code{\link{tune.iCluster2b}}, \code{\link{iClusterPlus2}}, \code{\link{iClusterBayes}}, \code{\link{plotHeatmap}}
}
\value{
  no value returned.
}
\examples{
# see tutorials 
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
# TBD 
}
\keyword{models}

