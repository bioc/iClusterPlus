\name{iClusterPlus2}
\alias{iClusterPlus2}
\title{Integrative clustering of multi-omics data}
\description{
Given multi-omics data (e.g., somatic mutation,copy number, gene
expression, DNA methylation) measured in the same set of samples, iClusterPlus2 fits a
regularized latent factor(variable) model that generates an latent factor matrix that can be used for integrative clustering of samples.
In addition, the driver features of sample clustering can be identified
by sparse coefficient matrices.
}
\usage{
iClusterPlus2(xList,type=c("gaussian","binomial","poisson","multinomial"),K=3,
                           n.burnin=100,n.draw=200,maxiter=25,sdev=0.05,lambda.scale=1/3,
                           BICrate.cutoff=rep(0.01,4),min.shrinkage.rate=rep(0.05,4))
}
\arguments{
  \item{xList}{A list object containing m data matrices representing m
    multi-omics data measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{type}{Data type, which can be gaussian, binomial, poisson, multinomial.}
  \item{K}{An integer representing the number of the latent factors used for modeling.}
  \item{n.burnin}{Number of MCMC burnin.}
  \item{n.draw}{Number of MCMC draw.}
  \item{maxiter}{Maximum iteration for the EM algorithm.}
  \item{sdev}{standard deviation of random walk proposal.}
  \item{lambda.scale}{scaling factor for the lasso regularization
    parameter when the data type is binomial.Empirically, lambda.scale
    is within (0.1, 1).}
  \item{BICrate.cutoff}{BIC rate between iteration i+1 and i.  If BIC
    rate < BICrate.cutoff, the search for optimal lambda will stop.}
  \item{min.shrinkage.rate}{The minimum lasso shrinkage rates for
    multi-omics features.}
}
\seealso{
 \code{\link{tune.iCluster2b}}, \code{\link{iClusterPlus}}, \code{\link{iClusterBayes}}
}
\value{
  A list with the following elements.
  \item{alpha}{A list of object for the intercept parameters corresponding to the m data matrices.}
  \item{beta}{A list object containing m coefficient matrices corresponding to the m data matrices.}
  \item{meanZ}{A n x k matrix; the rows represent samples and the columns represent the K factors.}
  \item{BIC}{Bayesian information criterion.}
  \item{dev.ratio}{see dev.ratio defined in glmnet package.}
  \item{lambda}{Final lasso regularization parameters used for iCluster modeling.}
}
\examples{
# see iManual.pdf
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2013).  
Pattern discovery and cancer gene identification in integrated cancer genomic data. 
\emph{Proc. Natl. Acad. Sci. USA.} 110(11):4245-50. 
}
\keyword{models}

