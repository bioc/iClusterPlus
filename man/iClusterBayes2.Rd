\name{iClusterBayes2}
\alias{iClusterBayes2}
\title{Integrative clustering of multi-omics data}
\description{
Given multi-omics data (e.g., Somatic mutation, copy number, gene
expression, DNA methylation) measured in the same set of samples, iClusterBayes2 fits a
Bayesian latent factor model that generates an latent factor matrix that
can be used for integrative clustering of samples. In addition, the
driver features of sample clustering can be identified by the posterior
probability of the model parameters.
}
\usage{
iClusterBayes2(xList,type = c("gaussian","binomial","poisson"),K=3,n.burnin=1000,n.draw=1200,
	prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5)
}
\arguments{
  \item{xList}{A list object containing m data matrices representing m
    multi-omics data measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.
  The allowed maximum number of data matrices is 6 (m < 7).}
  \item{type}{Data type corresponding to the data matrices in xList,
    which can be gaussian, binomial, or poisson.} 
  \item{K}{An integer representing the number of the latent factors used for modeling.}
  \item{n.burnin}{Number of MCMC burnin.}
  \item{n.draw}{Number of MCMC draw.}
  \item{prior.gamma}{Prior probability for the indicator variable gamma of each data set.}
  \item{sdev}{Standard deviation of random walk proposal for the latent variable.}
  \item{beta.var.scale}{A positive value to control the scale of covariance matrix of the proposed beta.}
  \item{thin}{A parameter to thin the MCMC chain in order to reduce autocorrelation.
  		Discard all but every 'thin'th sampling values. When thin=1, all sampling values are kept.}
  \item{pp.cutoff}{Posterior probability cutoff for the indicator variable gamma.  The BIC and deviance ratio will be
  	calculated by setting parameter beta to zero when the posterior probability of gamma <= cutoff.}
}
\seealso{
 \code{\link{iClusterPlus2}},\code{\link{plotHMBayes}},\code{\link{plotHeatmap}}
}
\value{
  A list with the following elements.
  \item{alpha}{A list of object for the intercept parameters corresponding to the m data matrices.}
  \item{beta}{A list object containing m coefficient matrices corresponding to the m data matrices.}
  \item{meanZ}{A n x k matrix; the rows represent samples and the columns represent the K factors.}
  \item{beta.pp}{A list of posterior probability for the parameter beta.
    The higher the beta.pp, the more likely the beta should be included
    in the model.} 
  \item{gamma.ar}{A list of acceptance ratio for the parameter gamma.}
  \item{beta.ar}{A list of acceptance ratio for the parameter beta.}
  \item{Z.ar}{Acceptance ratio for the latent variable.}
  \item{BIC}{Bayesian information criterion.}
  \item{dev.ratio}{ see dev.ratio defined in glmnet package.}
}
\examples{
# see iManual.pdf
}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
Mo Q, Shen R, Guo C, Vannucci M, Chan KS, Hilsenbeck SG. (2018).  
A fully Bayesian latent variable model for integrative clustering analysis of multi-type omics data.
Biostatistics 19(1):71-86. 
}
\keyword{models}
