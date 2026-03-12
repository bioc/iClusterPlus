# tune.iCluster2b find the optimal lasso lambda value for iCluster2b function 
# Last Updated: 11/19/2025 by Quincy Mo @ Moffitt Cancer Center
# library(irlba)
tune.iCluster2b = function(xList,K=3, method=c("lasso","enet"),min.lambda=10, max.lambda=500,
                           lambda.iter=25,EM.iter=25,min.shrinkage.rate=rep(0.05,length(xList)), 
                           eps=1e-4, eps2=1e-6){
  ## xList is a list of data matrix, x[[i]] 
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
    stop("xList must contain at least one data matrix. \n")
  }

  nrow1 = nrow(xList[[1]])
  if(ns >1){
    for(i in 2:ns){
      if(nrow(xList[[i]]) != nrow1){
        stop("The numbers of rows of the matrices in xList are not equal")
      }
    }
  }
  
  if(min.lambda < 0 || max.lambda < 0 ){
    stop("min.lambda and max.lambda must be greater than 0. \n")
  }
  
  if(min.lambda >= max.lambda){
    stop("max.lambda must be greater than min.lambda. \n")
  }

  endpoints = 0
  ID = 1
  lenID = 1
  lbda = list()
  lbdIncr = (max.lambda - min.lambda)/lambda.iter
  #lbda = lbda + lbdIncr
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
    # set the initial lambda value 
    if(methodN[i] == 1){
      #lambda for lasso must be a single value
      lbda[[i]] = min.lambda
    }else if(methodN[i] == 2){
      #lambda for enet must be a vector with two elements
      lbda[[i]] = c(min.lambda, min.lambda)
    }else if(methodN[i] == 3){
      #lambda for flasso must be a vector with two elements
      lbda[[i]] = c(min.lambda, min.lambda)
    }else if(methodN[i] == 4){
     #lambda for glasso must be a single value
      lbda[[i]] = min.lambda
    }else if(methodN[i] == 5){
     #lambda for gflasso must be a vector with two elements
      lbda[[i]] = c(min.lambda, min.lambda)
    }
  }

  X = NULL
  x = list()
  for (i in 1:ns) {
    x[[i]] = scale(xList[[i]],center=TRUE,scale=TRUE)
    X = cbind(X, x[[i]])		
  }
  n = dim(X)[1]
  p = dim(X)[2]
  XtX = t(X)%*%X
  XtXdiag = diag(XtX)

  if(K < 1 || K > min(dim(X))){
    stop("K must be an integer greater than 1 and less than the rank of the combined matrix. \n")
  }
  
  pcaX = prcomp_irlba(X, n=K,center=FALSE,scale.=FALSE)
  #pcaX = prcomp(X, center=FALSE,scale.=FALSE)
  B = pcaX$rotation[,1:K] ## use PCA for B
  EZ0 = t(pcaX$x[,1:K])
  EZZt0 = EZ0 %*% t(EZ0)
  #EZ0 = matrix(0,nrow=K,ncol=n)
  #EZZt0 = matrix(0,K,K)
  
  Bmat.list0=split(data.frame(B),f=rep(1:ns,pvec))
  Xdiff = X - pcaX$x[,1:K] %*% t(B)
  Phivec0 = diag(t(Xdiff)%*%Xdiff)/(n-1)
  Phi.list0 = split(Phivec0,f=rep(1:ns,pvec))
  
  minBIC = rep(NA,ns)
  for(i in 1:ns){
    fit1 = .C("logNormAll",loglike = double(1),as.double(pcaX$x[,1:K]),as.double(rep(0,pvec[i])),
                as.double(unlist(Bmat.list0[[i]])),as.double(unlist(Phi.list0[[i]])),as.double(x[[i]]),
                as.integer(n),as.integer(pvec[i]),as.integer(K),PACKAGE="iClusterPlus")
    minBIC[i] = -2*fit1$loglike + sum(abs(Bmat.list0[[i]]) > eps2)*log(n * pvec[i])
  }
  #cat(minBIC,"\n")
  
  bestEZ = EZ0
  bestBeta = Bmat.list0
  
  outIter = 0
  BICrate = rep(1,ns) #
  while(outIter < lambda.iter && any(BICrate != 0)){
    outIter = outIter+1
    #cat(outIter," : ", unlist(lbda), "\n")
    cat(outIter," ")
    
    iter = 0  #initial value for iClusterCore 
    dif = 1
    ires = .C("iClusterCore",as.integer(p),as.integer(K),as.integer(n), as.double(XtXdiag),as.double(X),
      B=as.double(B),EZ=as.double(EZ0),EZZt=as.double(EZZt0),Phivec=as.double(Phivec0),dif=as.double(dif),
      iter = as.integer(iter),as.integer(pvec),as.double(unlist(lbda)),as.double(eps),as.double(eps2),
      as.integer(EM.iter),as.integer(ns),as.integer(methodN),as.integer(ID),as.integer(lenID),
      PACKAGE="iClusterPlus")
  
    Bmat = matrix(ires$B,ncol=K)
    Bmat[abs(Bmat) <= eps2] = 0
    EZ = matrix(ires$EZ,ncol=n)
    EZZt = matrix(ires$EZZt,ncol=K)
    Bmat.list=split(data.frame(Bmat),f=rep(1:ns,pvec))
    Phi.list = split(ires$Phivec,f=rep(1:ns,pvec)) 
    # some values can become < 0 for glasso, cause loglike = NaN
    # use Phi.list0 for calculating loglike
    for(i in 1:ns){
      fit1 = .C("logNormAll",loglike = double(1),as.double(t(EZ)),as.double(rep(0,pvec[i])),
                as.double(as.matrix(Bmat.list[[i]])),as.double(as.matrix(Phi.list0[[i]])),as.double(x[[i]]),
                as.integer(n),as.integer(pvec[i]),as.integer(K),PACKAGE="iClusterPlus")
      tempBIC = -2*fit1$loglike + sum(Bmat.list[[i]] != 0)*log(n * pvec[i])
      BICrate[i] = (tempBIC - minBIC[i])/minBIC[i]
      beta0rate = sum(apply(Bmat.list[[i]],1,function(x){all(x==0)}))/pvec[i]
      #cat("BICrate ",i, ":", BICrate[i], " labmbda: ", lbda[[i]], " beta0rate ", beta0rate,"\n")
      if(beta0rate < min.shrinkage.rate[i]){
        lbda[[i]] = lbda[[i]] + lbdIncr
      }else if(tempBIC > minBIC[i]){
       lbda[[i]] = lbda[[i]] + lbdIncr
      }else if(tempBIC < minBIC[i]){
        lbda[[i]] = lbda[[i]] - lbdIncr   
        minBIC[i] = tempBIC
        bestEZ = EZ
        bestBeta = Bmat.list
      }
      
      if(max(lbda[[i]]) < min.lambda){
        if(methodN[i] == 1 || methodN[i] == 3){
          lbda[[i]] = min.lambda
        }else{
          lbda[[i]] = c(min.lambda,min.lambda)
        }
      }else if(min(lbda[[i]]) > max.lambda){
        if(methodN[i] == 1 || methodN[i] == 3){
          lbda[[i]] = max.lambda
        }else{
          lbda[[i]] = c(max.lambda,max.lambda)
        }
      }
    }
  }
  cat("\n")
  
  EZ = t(EZ)
  rownames(EZ) = rownames(xList[[1]])
  colnames(EZ) = paste0("Factor",1:K)
  
  for (i in 1:ns) {
    rownames(bestBeta[[i]]) = colnames(xList[[i]])
    colnames(bestBeta[[i]]) = paste0("Factor",1:K)
  }
  
  return(list(meanZ=EZ, beta=bestBeta,lambda=lbda))
}

