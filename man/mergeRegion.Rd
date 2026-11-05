\name{mergeRegion}
\alias{mergeRegion}
\title{A function to combine neighboring genomic regions}
\description{
This function is used to combine neighboring genomic regions. 
}
\usage{
mergeRegion(peak.ct, rm.prob=0.75,rm.isolate=FALSE, min.gapwidth=10000)
}
\arguments{
  \item{peak.ct}{A sparse count Matrix of class "dgCMatrix". Rows are chromosome regions with format chrX:123-456. Columns are sample/cell IDs.}
  \item{rm.prob}{A value in (0, 1).  This value is used to control removing regions with low counts. For example, rm.prob = 0.75 means that 75\% of the least variable regions will be filtered out.}
  \item{rm.isolate}{logical; if TRUE, isolated regions will be removed. Isolated regions are regions without neighboring regions within min.gapwidth.}
  \item{min.gapwidth}{An integer. Neighboring regions less than min.gapwidth will be merged into a single regions.}
}
\seealso{
 \code{\link{CNregions}},\code{\link{iCluster2b}},\code{\link{tune.iCluster2b}},\code{\link{iClusterPlus2}}
}
\value{
A sparse matrix with merged genomic regions. The value in a merged region is the average counts of the merged regions.
}
\examples{

# see the example in the tutorial of iCluster analysis of scRNA-seq + scATAC-seq data
# atac <- mergeRegion(pbmc[["ATAC"]]$data,rm.prob=0.8,rm.isolate=FALSE,min.gapwidth=10000)

}
\author{Qianxing Mo \email{qianxing.mo@moffitt.org}}
\references{
Qianxing Mo, Sijian Wang, Venkatraman E. Seshan, Adam B. Olshen, Nikolaus Schultz, Chris Sander, R. Scott Powers, Marc Ladanyi, and Ronglai Shen.  (2013).  
Pattern discovery and cancer gene identification in integrated cancer genomic data. 
Proc. Natl. Acad. Sci. USA.
}
\keyword{models}


