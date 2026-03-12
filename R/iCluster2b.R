# iCluster version 2b for fast decomposition of large matrix to get initial values for iCluster2 modeling
# Last Updated: 11/01/2025 by Quincy Mo @ Moffitt Cancer center
# library(irlba)
iCluster2b = function(xList,K=3, lambda, method=c("lasso","enet","flasso","glasso","gflasso"),
  chr=NULL, EM.iter=25, eps=1e-4, eps2=1e-6){
  ## xList (small) is the array, xList[[k]]=...; rows are samples and columns are genes 
  ## X (big) is the pool-together matrix
  ## endpoints is the endpoints of each chromosome chains (for fused lasso)
  ## for example, if endpoints=c(10,20,30), it means there are three chains: x1:x10, x11:x20, x21:x30
  ## if there is only one chain, endpoints=p
  ## dim(X)=n*p
  ## X = (X1,...,Xn)^T, Xi is a p*1 column vector of the observation for subject i
  ## Z = (Z1,...,Zn)^T, Zi is a K*1 column vector of the latent variable for subject i
  ## B = (B1,...,Bp)^T, Bj is a K*1 column vector of coefficient for variable j
  
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
  
  if(any(unlist(lambda) < 0)){
    stop("Values in lambda must be greater than 0. \n")
  }
  
  endpoints = 0
  ID = 1
  lenID = 1
  lbda = list()
  allMethod = c("lasso","enet","flasso","glasso","gflasso")
  methodN = rep(NA,ns)
  pvec = rep(0,ns)
  for(i in 1:ns){
    pvec[i] = dim(xList[[i]])[2]
    methodN[i] = pmatch(method[i],allMethod)
    if(is.na(methodN[i])){
      stop("Methods are not correctly specified!\n")
    }

    if((methodN[i] == 3) || (methodN[i] == 5)){
      if(is.null(chr)){
        stop("User must supply chromosome indicator when using flasso or gflasso")
      }
      endpoints = cumsum(table(chr))
      ID = setdiff((1:pvec[i]), endpoints)
      lenID = length(ID)
    }

    if(methodN[i] == 1){
      if(length(lambda[[i]])==1){
        lbda[[i]] = lambda[[i]]
      }else{
        stop("Error in lambda: lambda for lasso must be single value !\n")
      }
    }else if(methodN[i] == 2){
      if(length(lambda[[i]])==2){
        lbda[[i]] = lambda[[i]]
      }else{
        stop("Error in lambda: lambda for enet must be a vector with two elements!\n")
      }
    }else if(methodN[i] == 3){
      if(length(lambda[[i]])==2){ 
        lbda[[i]] = lambda[[i]]
      }else{
        stop("Error in lambda: lambda for flasso must be a vector with two elements!\n")
      }
    }else if(methodN[i] == 4){
      if(length(lambda[[i]])==1){
        lbda[[i]] = lambda[[i]]
      }else{
        stop("Error in lambda: lambda for glasso must be a single value!\n")
      }
    }else if(methodN[i] == 5){
      if(length(lambda[[i]])==2){ 
        lbda[[i]] = lambda[[i]]
      }else{
        stop("Error in lambda: lambda for gflasso must be a vector with two elements!\n")
      }
    }
  }
  
  lbda = unlist(lbda)
  X = NULL
  for (i in 1:ns) {
    x = scale(xList[[i]],center=TRUE,scale=TRUE)
    X = cbind(X, x)		
  }
  n = dim(X)[1]
  p = dim(X)[2]
  XtX = t(X)%*%X
  XtXdiag = diag(XtX)
  
  if(K < 1 || K > min(dim(X))){
    stop("K must be an integer greater than 1 and less than the rank of the combined matrix. \n")
  }
  
  ## Initial value using K PCs
  pcaX = prcomp_irlba(X, n=K,center=FALSE,scale.=FALSE)
  #pcaX = prcomp(X, center=FALSE,scale.=FALSE)
  B = pcaX$rotation[,1:K] ## use PCA for B
  EZ = t(pcaX$x[,1:K])
  EZZt = EZ %*% t(EZ)
  #EZ = matrix(0,nrow=K,ncol=n)
  #EZZt = matrix(0,K,K)

  Xdiff = X - pcaX$x[,1:K] %*% t(B)
  Phivec = diag(t(X)%*%X)/(n-1)
  
  iter = 0
  dif = 1
  
  ires = .C("iClusterCore",as.integer(p),as.integer(K),as.integer(n), as.double(XtXdiag),as.double(X),
    B=as.double(B),EZ=as.double(EZ),EZZt=as.double(EZZt),Phivec=as.double(Phivec),dif=as.double(dif),
    iter = as.integer(iter),as.integer(pvec),as.double(lbda),as.double(eps),as.double(eps2),
    as.integer(EM.iter),as.integer(ns),as.integer(methodN),as.integer(ID),as.integer(lenID),
    PACKAGE="iClusterPlus")
        
  Bmat = matrix(ires$B,ncol=K)
  Bmat[abs(Bmat) <= eps2] = 0
  #EZZt = matrix(ires$EZZt,ncol=K)
  
  EZ = t(matrix(ires$EZ,ncol=n))
  rownames(EZ) = rownames(xList[[1]])
  colnames(EZ) = paste0("Factor",1:K)
  
  Bmat.list=split(data.frame(Bmat),f=rep(1:ns,pvec))
  for (i in 1:ns) {
    rownames(Bmat.list[[i]]) = colnames(xList[[i]])
    colnames(Bmat.list[[i]]) = paste0("Factor",1:K)
  }
  
  return(list(meanZ=EZ,beta=Bmat.list,iter=ires$iter))
}

