
#xList is a matrix (rows for samples and columns for variables)
pcaVarPlot = function(xList,K=10){
  ns = length(xList)
  if(ns == 0){
    stop("xList must contain at least one data matrix! \n")
  }
  
  nrow1 = nrow(xList[[1]])
  if(ns >1){
    for(i in 2:ns){
      if(nrow(xList[[i]]) != nrow1){
        stop("The numbers of rows of the matrices in xList are not equal")
      }
    }
  }

  X = NULL
  for (i in 1:ns) {
    x = scale(xList[[i]],center=TRUE,scale=TRUE)
    X = cbind(X, x)	
  }
  
  if(K < 1 || K > min(dim(X))){
    stop("K must be an integer greater than 1 and less than the rank of the combined matrix. \n")
  }
  
  pcaX = prcomp_irlba(X, n=K,center=FALSE,scale.=FALSE)
  pca.var.prop = pcaX$sdev^2/pcaX$totalvar
  plot(pca.var.prop, type="b",xlab="PC",ylab="Proportion of total variance")
}
