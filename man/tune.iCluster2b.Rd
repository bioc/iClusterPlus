\name{tune.iCluster2b}
\alias{tune.iCluster2b}
\title{Integrative clustering of multi-omics data}
\description{
This function finds optimal lasso regularization parameters for
iCluster2b.  
}
\usage{
tune.iCluster2b(xList,K=3,method=c("lasso","enet"),min.lambda=10,max.lambda=500,
    lambda.iter=25,EM.iter=25,min.shrinkage.rate=rep(0.05,length(xList)), eps=1e-4, eps2=1e-6)
}
\arguments{
  \item{xList}{A list object containing m data matrices representing m
    multi-omics data measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{K}{An positive integer representing the number of the latent factors used for modeling.}
  \item{method}{Method used for the regularization of model parameters.}
  \item{min.lambda}{The minimum value of the lasso regularization parameter.}
  \item{max.lambda}{The maximum value of the lasso regularization parameter.}
  \item{lambda.iter}{the number of iteration to find an optimal lambda.}
  \item{EM.iter}{The number of iteration for the EM algorithm.}
  \item{min.shrinkage.rate}{The minimum lasso shrinkage rates for
    multi-omics features.}
  \item{eps}{EM algorithm convergence criterion 1.}
  \item{eps2}{EM algorithm convergence criterion 2.}
}
\seealso{
   \code{\link{iCluster2b}}, \code{\link{iClusterPlus2}}, \code{\link{iClusterBayes}}
}
\value{
  A list with the following elements.
  \item{meanZ}{A n x k matrix; the rows represent samples and the columns represent the K factors.}
  \item{beta}{A list object containing m coefficient matrices
    corresponding to the m data matrices.}
  \item{lambda}{Final lasso regularization parameters used for iCluster modeling.}
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
  Ronglai Shen, Sijian Wang, Qianxing Mo. (2013).  Sparse Integrative Clustering of Multiple Omics Data Sets. \emph{Annals of Applied Statistics.} 7(1):269-294
}
\keyword{models}

