# iClusterBayes v2, updated by Quincy Mo @ Moffitt Cancer Center

iClusterBayes2 <- function(xList,type = c("gaussian","binomial","poisson"),K=3,
                          n.burnin=1000,n.draw=1200,prior.gamma=rep(0.1,6),sdev=0.5,beta.var.scale=1,thin=1,pp.cutoff=0.5){

  ndt = length(xList)
  dttype = c("gaussian","binomial","poisson")
  n = nrow(xList[[1]]) # sample size
  
  if(ndt == 0){
    stop("xList must contain at least one data matrix! \n")
  }
  pvec = rep(NA,ndt)
  
  if(ndt >1){
    for(i in 2:ndt){
      if(nrow(xList[[i]]) != n){
        stop("The numbers of rows of the matrices in xList are not equal")
      }
    }
  }
  
  if(ndt > 6){
    stop("xList is only allowed to contain <=6 data matrices. \n")
  }
  
  if(!all(type %in% dttype)){
      cat("Error: ",type[!all(type %in% dttype)],"\n")
      stop("Allowed data types are gaussian, binomial and poisson. \n")
  }

  if(ndt != length(type)){
    stop("The numbers of data matrices and data types are inconsistent. \n
         Data type must be specified correctly. \n")
  } 
  
  x = list()
  X = NULL
  for (i in 1:ndt) {
    if((type[i] == "gaussian")){
      x[[i]] = scale(xList[[i]],center=TRUE,scale=TRUE)
    }else if((type[i] == "poisson")){
      x[[i]] = scale(xList[[i]],center=TRUE,scale=TRUE)
      type[i] = "gaussian" ## convert to Gaussian for computational efficiency
    }else{
      x[[i]] = xList[[i]]
    }
    X = cbind(X, x[[i]])
    pvec[i] = ncol(x[[i]])
  }

  if(K < 1 | K > min(dim(X))){
    stop("K must be an integer greater than 1 and less than the rank of combined data matrices. \n")
  }
  
  #pcaX = prcomp_irlba(X,n=K,center=FALSE,scale.=FALSE)
  #initialZ = pcaX$x ## use K PCA for initial latent variable 
  #Bmat.list = split(data.frame(pcaX$rotation), f = rep(1:ndt, pvec))
  
  ### burnin and draw for latent variable Z ##
  zBurnin = 0
  zDraw = 1
  ### the above settings make the draws for Z are the same as the draws for beta
  ### if zDraw > 1, the outcome Z from mcmcMix will be the mean Z of zDraw
  
  ### n.burnin and n.draw control the overall draw for Z and parameter beta
  betaBurnin = n.burnin
  betaDraw = n.draw
  
  Data = list()
  for(i in 1:ndt){
    Data[[i]] = dataType(x[[i]],type[i],K) #normalized data produce better results
    #Data[[i]]$Beta = as.matrix(Bmat.list[[i]])
  }
  
  ## if number of data matrix < 4, set Data[[ndt+1]] ... to Data[[4]] be NULL
  if(ndt < 6){
    for(i in (ndt+1):6){
      Data[[i]] = NULL
    }
  }

  ###### priors for Bayesian variable selection #######
  invSigma0 = diag(rep(1,K+1))
  beta0 = rep(0,K+1)
  invSigmaBeta0 = invSigma0 %*% beta0
  invga0=1
  invgb0=1

  res = mcmcBayes(dt1=Data[[1]],dt2=Data[[2]],dt3=Data[[3]],dt4=Data[[4]],dt5=Data[[5]],dt6=Data[[6]],ndt,sdev,n,K,zBurnin,zDraw,
    betaBurnin, betaDraw,thin,prior.gamma,beta0,invSigma0,invSigmaBeta0,invga0,invgb0,beta.var.scale)

  for(i in 1:ndt){
    names(res$Alpha[[i]]) = colnames(xList[[i]])
    rownames(res$Beta[[i]]) = colnames(xList[[i]])
    colnames(res$Beta[[i]]) = paste0("Factor",1:K)   
    names(res$Ratio[[i]]) = colnames(xList[[i]])
    names(res$acsGamma[[i]]) = colnames(xList[[i]])
    names(res$acsBeta[[i]]) = colnames(xList[[i]])
    Data[[i]]$Alpha = res$Alpha[[i]]
    Data[[i]]$Beta = res$Beta[[i]]
    Data[[i]]$sigma2 = res$Sigma2[[i]]
    Data[[i]]$Ratio = res$Ratio[[i]]
    Data[[i]]$acsGamma = res$acsGamma[[i]]
    Data[[i]]$acsBeta = res$acsBeta[[i]]
  }
  BIC = totalBICbayes(Data, res$meanZ, ndt, K, pp.cutoff)
  devRatio = dev.ratio.bayes(Data, res$meanZ,ndt, K, pp.cutoff)
  rownames(res$meanZ) = rownames(xList[[1]])
  colnames(res$meanZ) = paste0("Factor",1:K)
  list(alpha=res$Alpha, beta=res$Beta, meanZ=res$meanZ,beta.pp=res$Ratio,BIC=BIC, dev.ratio=devRatio,
       gamma.ar=res$acsGamma,beta.ar=res$acsBeta,Z.ar=res$acsZ)
}
