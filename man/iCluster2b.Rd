\name{iCluster2b}
\alias{iCluster2b}
\title{Integrative clustering of multi-omics data}
\description{
Given multi-omics data (e.g., copy number, gene
expression, DNA methylation) measured in the same set of samples, iCluster2b fits a
regularized latent factor(variable) model that generates an latent factor matrix that can be used for integrative clustering of samples.
In addition, the driver features of sample clustering can be identified by sparse coefficient matrices. 
}
\usage{
iCluster2b(xList,K=3, lambda, method=c("lasso","enet","flasso","glasso","gflasso"),
  chr=NULL, EM.iter=25, eps=1e-4, eps2=1e-6)
}
\arguments{
  \item{xList}{A list object containing m data matrices representing m
    multi-omics data measured in a set of n samples. For each matrix, the
    rows represent samples, and the columns represent genomic features.}
  \item{K}{An integer representing the number of the latent factors used for modeling.}
  \item{lambda}{A list with m elements corresponding to the method
    argument. For lasso and glasso, each element is a positive value.
    For enet, flasso and gflasso, each element is a vector with two
    positive values.} 
  \item{method}{Method used for variable selection.}
  \item{chr}{Chromosome labels}
  \item{EM.iter}{Maximum iteration for the EM algorithm.}
  \item{eps}{EM algorithm convergence criterion 1.}
  \item{eps2}{EM algorithm convergence criterion 2.}  
}
\seealso{
 \code{\link{tune.iCluster2b}}, \code{\link{iClusterPlus2}}, \code{\link{iClusterBayes}}
}
\value{
  A list with the following elements.
  \item{meanZ}{A n x k matrix; the rows represent samples and the columns represent the K factors.}
  \item{beta}{A list object containing m coefficient matrices corresponding to the m data matrices.}	
  \item{iter}{EM iteration.}
}
\examples{

## clustering
n1 = 20
n2 = 20
n3 = 20
n = n1+n2+n3
p = 5
q = 100

x = NULL
x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[1]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[2]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[3]] = cbind(xa,xb)


x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[4]] = cbind(xa,xb)

x1a = matrix(rnorm(n1*p), ncol=p)
x2a = matrix(rnorm(n1*p, -1.5,1), ncol=p)
x3a = matrix(rnorm(n1*p, 1.5, 1), ncol=p)
xa = rbind(x1a,x2a,x3a)
xb = matrix(rnorm(n*q), ncol=q)
x[[5]] = cbind(xa,xb)

method = c('lasso', 'enet', 'flasso', 'glasso', 'gflasso')  
lambda=alist()
lambda[[1]] = 30
lambda[[2]] = c(20,1)
lambda[[3]] = c(20,20)
lambda[[4]] = 30
lambda[[5]] = c(30,20)

chr=c(rep(1,10),rep(2,(p+q)-10)) 
date()
fit2 = iCluster2b(x, K=3, lambda, method=method, chr=chr, EM.iter=20,eps=1e-4, eps2=1e-6)
date()

par(mfrow=c(5,1),mar=c(4,4,1,1))
for(i in 1:5){
barplot(fit2$beta[[i]][,1])
}

}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org},Ronglai Shen,Sijian Wang}
\references{
  Ronglai Shen, Sijian Wang, Qianxing Mo. (2013).  Sparse Integrative Clustering of Multiple Omics Data Sets. \emph{Annals of Applied Statistics.} 7(1):269-294
}
\keyword{models}

