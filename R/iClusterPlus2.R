# iClusterPlus v2, updated by Quincy Mo @ Moffitt Cancer Center

iClusterPlus2 = function(xList,type=c("gaussian","binomial","poisson","multinomial"),K=3,
                           n.burnin=100,n.draw=200,maxiter=25,sdev=0.05,lambda.scale=1/3,
                         BICrate.cutoff=rep(0.01,4),min.shrinkage.rate=rep(0.05,4)){
  alpha=c(1,1,1,1) #Vector of elasticnet penalty terms; set to 1 because elasticnet has not been implemented.
  ndt = length(xList)
  dttype = c("gaussian","binomial","poisson","multinomial")
  n = nrow(xList[[1]]) # sample size
  
  if(ndt == 0){
    stop("xList must contain at least one data matrix! \n")
  }

  if(ndt >1){
    for(i in 2:ndt){
      if(nrow(xList[[i]]) != n){
        stop("The numbers of rows of the matrices in xList are not equal")
      }
    }
  }
  
  if(ndt > 4){
    stop("xList is only allowed to contain <=4 data matrices. \n")
  }
  
  if(!all(type %in% dttype)){
    cat("Error: ",type[!all(type %in% dttype)],"\n")
    stop("Allowed data types are gaussian, binomial, poisson and multinomial. \n")
  }
  
  if(ndt != length(type)){
    stop("The numbers of data matrices and data types are inconsistent. \n
         Data type must be specified correctly. \n")
  } 

  if(any(lambda.scale < 0.1) | (any(lambda.scale > 1))){
    stop("lambda.scale is expected to be a small value between 0.1 and 1 \n")   
  }
  
  if(any(BICrate.cutoff < 0.005) | (any(BICrate.cutoff > 0.05))){
    stop("BICrate.cutoff is expected to be a small value between 0.005 and 0.05 \n")   
  }
  
  x = list()
  for (i in 1:ndt) {
    if((type[i] == "gaussian")){
      x[[i]] = scale(xList[[i]],center=TRUE,scale=TRUE)
    }else if((type[i] == "poisson")){
      x[[i]] = scale(xList[[i]],center=TRUE,scale=TRUE)
      type[i] = "gaussian" ## convert to Gaussian for computational efficiency
    }else{
      x[[i]] = xList[[i]]
    }
  }

  X = NULL
  for (i in 1:ndt) {
    if((type[i] != "multinomial")){
      X = cbind(X, x[[i]])		
    }
  }

  if(K < 1 || K > min(dim(X))){
    stop("K must be an integer greater than 1 and less than the rank of the combined matrix. \n")
  }
  
  ## initial value of latent factor matrix
  lastZ = NULL
  if((ndt == 1) & (type[1] == "multinomial")){
    lastZ = matrix(rnorm(n*K,0,1),ncol=K)  
  }else{
    pcaX = prcomp_irlba(X,n=K,center=FALSE,scale.=FALSE)
    #pcaX = prcomp(X, center=FALSE,scale.=FALSE)
    lastZ = pcaX$x[,1:K] ## use K PCA for initial latent variable 
  }
  
  Alpha = list()
  Beta = list()
  Data = list()

  for(i in 1:ndt){
    Data[[i]] = dataType(x[[i]],type[i],K) #normalized data produce better results
    Alpha[[i]] = Data[[i]]$Alpha
    Beta[[i]] = Data[[i]]$Beta   
  }
  
  ## if number of data matrix < 4, set Data[[ndt+1]] ... to Data[[4]] be NULL
  if(ndt < 4){
    for(i in (ndt+1):4){
      Data[[i]] = NULL
    }
  }
  
  ### calculate the initial BIC values  ###
  minBIC = rep(NA,ndt)
  for(i in 1:ndt){
    if(Data[[i]]$type == 1){  # normal #
      fit1 = .C("logNormAll",loglike = double(1),as.double(lastZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.double(Data[[i]]$sigma2),as.double(Data[[i]]$con),
                as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      minBIC[i] = -2*fit1$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
    }else if(Data[[i]]$type == 2){ # binomial #
      fit2 = .C("logBinomAll",loglike = double(1),as.double(lastZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      minBIC[i] = -2*fit2$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
    }else if(Data[[i]]$type == 3){ # Poisson #
      fit3 = .C("logPoissonAll",loglike = double(1),as.double(lastZ),as.double(Data[[i]]$Alpha),
                as.double(Data[[i]]$Beta),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
      minBIC[i] = -2*fit3$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
    }else { # Multinomial #
      fit4 = .C("logMultAll",loglike = double(1),as.double(lastZ),as.double(Data[[i]]$Alpha),
                as.double(t(matrix(Data[[i]]$Beta,ncol=K*Data[[i]]$C))),as.integer(Data[[i]]$cat),
                as.integer(Data[[i]]$class),as.integer(Data[[i]]$nclass),as.integer(Data[[i]]$n),
                as.integer(Data[[i]]$p),as.integer(Data[[i]]$C),as.integer(K),PACKAGE="iClusterPlus")
      minBIC[i] = -2*fit4$loglike + sum(Data[[i]]$Beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
    }
  }
  
  ###### estimate the latent variable ########
  iter = 0
  lambda = rep(0.01,ndt)
  BICrate = rep(1,ndt)
  lambda.inc = 1/maxiter
  newZ = list(meanZ=lastZ)
  while(iter<maxiter && any(abs(BICrate) > 0.01)){
    iter = iter + 1
    cat(iter, " ") 
    
    for(i in 1:ndt){
      if(Data[[i]]$type == 1){  # normal #
        fit = .C("elnetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),sigma2 = double(Data[[i]]$p),
                 as.double(newZ$meanZ),as.double(Data[[i]]$con),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
                 as.double(alpha[i]),as.double(lambda[i]),PACKAGE="iClusterPlus")
        fit1 = .C("logNormAll",loglike = double(1),as.double(newZ$meanZ),as.double(fit$a0),
                  as.double(fit$beta),as.double(fit$sigma2),as.double(Data[[i]]$con),
                  as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
        tempBIC = -2*fit1$loglike + sum(fit$beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
        BICrate[i] = (tempBIC - minBIC[i])/minBIC[i]
        beta0rate = sum(apply(matrix(fit$beta,ncol=K),1,function(x){all(x==0)}))/Data[[i]]$p
        #cat("BICrate ",i, ":", BICrate[i], " labmbda: ", lambda[i], " beta0rate ", beta0rate,"\n")
        if(beta0rate < min.shrinkage.rate[i]){
          lambda[i] = lambda[i] + lambda.inc 
        }else if(BICrate[i] > BICrate.cutoff[i]){ 
          lambda[i] = lambda[i] + lambda.inc 
        }else if(BICrate[i] < -BICrate.cutoff[i]){
          minBIC[i] = tempBIC
          lambda[i] = lambda[i] - lambda.inc           
        }
        if(lambda[i] < 0.01){
          lambda[i] = 0.01
        }else if(lambda[i] > 0.99){
          lambda[i] = 0.99
        }
        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
        Data[[i]]$sigma2 = fit$sigma2
      }else if(Data[[i]]$type == 2){ # binomial #
        fit = .C("lognetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(newZ$meanZ),
                 as.integer(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
                 as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),as.integer(2),
                 as.integer(0),PACKAGE="iClusterPlus") #family=0 is binomial
        fit2 = .C("logBinomAll",loglike = double(1),as.double(newZ$meanZ),as.double(fit$a0),
                  as.double(fit$beta),as.integer(Data[[i]]$cat),
                  as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
        tempBIC = -2*fit2$loglike + sum(fit$beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
        BICrate[i] = (tempBIC - minBIC[i])/minBIC[i]
        beta0rate = sum(apply(matrix(fit$beta,ncol=K),1,function(x){all(x==0)}))/Data[[i]]$p
        #cat("BICrate ",i, ":", BICrate[i], " labmbda: ", lambda[i], " beta0rate ", beta0rate,"\n")
        if(beta0rate < min.shrinkage.rate[i]){
          lambda[i] = lambda[i] + lambda.inc*lambda.scale
        }else if(BICrate[i] > BICrate.cutoff[i]){
          lambda[i] = lambda[i] + lambda.inc*lambda.scale
        }else if(BICrate[i] < -BICrate.cutoff[i]){
          minBIC[i] = tempBIC
          lambda[i] = lambda[i] - lambda.inc*lambda.scale          
        }
        if(lambda[i] < 0.01){
          lambda[i] = 0.01
        }else if(lambda[i] > 0.99){
          lambda[i] = 0.99
        }
        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
      }else if(Data[[i]]$type == 3){ # Poisson #
        fit = .C("fishnetBatch",a0=double(Data[[i]]$p),beta=double(Data[[i]]$p*K),as.double(newZ$meanZ),
                 as.double(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
                 as.double(alpha[i]),as.double(lambda[i]),PACKAGE="iClusterPlus")
        fit3 = .C("logPoissonAll",loglike = double(1),as.double(newZ$meanZ),as.double(fit$a0),
                  as.double(fit$beta),as.integer(Data[[i]]$cat),
                  as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),as.integer(K),PACKAGE="iClusterPlus")
        tempBIC = -2*fit3$loglike + sum(fit$beta != 0)*log(Data[[i]]$n * Data[[i]]$p) 
        BICrate[i] = (tempBIC - minBIC[i])/minBIC[i]
        beta0rate = sum(apply(matrix(fit$beta,ncol=K),1,function(x){all(x==0)}))/Data[[i]]$p
        #cat("BICrate ",i, ":", BICrate[i], " labmbda: ", lambda[i], " ncoef= ", sum(fit$beta != 0),"\n")
        if(beta0rate < min.shrinkage.rate[i]){
          lambda[i] = lambda[i] + lambda.inc 
        }else if(BICrate[i] > BICrate.cutoff[i]){
          lambda[i] = lambda[i] + lambda.inc
        }else if(BICrate[i] < -BICrate.cutoff[i]){
          minBIC[i] = tempBIC
          lambda[i] = lambda[i] - lambda.inc          
        }
        if(lambda[i] < 0.01){
          lambda[i] = 0.01
        }else if(lambda[i] > 0.99){
          lambda[i] = 0.99
        }
        Data[[i]]$Alpha = fit$a0
        Data[[i]]$Beta = matrix(fit$beta,ncol=K)
      }else { # Multinomial #
        fit = .C("lognetBatch",a0=double(Data[[i]]$p * Data[[i]]$C),beta=double(Data[[i]]$p*K*Data[[i]]$C),
                 as.double(newZ$meanZ),as.integer(Data[[i]]$cat),as.integer(n),as.integer(K),as.integer(Data[[i]]$p),
                 as.double(alpha[i]),as.double(lambda[i]),as.integer(Data[[i]]$nclass),
                 as.integer(Data[[i]]$C),as.integer(1),PACKAGE="iClusterPlus") #family=1 is multinomial
        fit4 = .C("logMultAll",loglike = double(1),as.double(newZ$meanZ),as.double(fit$a0),
                  as.double(t(matrix(fit$beta,ncol=K*Data[[i]]$C))),as.integer(Data[[i]]$cat),as.integer(Data[[i]]$class),
                  as.integer(Data[[i]]$nclass),as.integer(Data[[i]]$n),as.integer(Data[[i]]$p),
                  as.integer(Data[[i]]$C),as.integer(K),PACKAGE="iClusterPlus")
        tempBIC = -2*fit4$loglike + sum(fit$beta != 0)*log(Data[[i]]$n * Data[[i]]$p)
        BICrate[i] = (tempBIC - minBIC[i] )/minBIC[i]
        beta0rate = sum(apply(matrix(fit$beta,ncol=K*Data[[i]]$C),1,function(x){all(x==0)}))/Data[[i]]$p #need verification
        #cat("BICrate ",i, ":", BICrate[i], " labmbda: ", lambda[i], " ncoef= ", sum(fit$beta != 0),"\n")
        if(beta0rate < min.shrinkage.rate[i]){
          lambda[i] = lambda[i] + lambda.inc*lambda.scale 
        }else if(BICrate[i] > BICrate.cutoff[i]){
          lambda[i] = lambda[i] + lambda.inc*lambda.scale
        }else if(BICrate[i] < -BICrate.cutoff[i]){
          minBIC[i] = tempBIC
          lambda[i] = lambda[i] - lambda.inc*lambda.scale          
        }  
        if(lambda[i] < 0.01){
          lambda[i] = 0.01
        }else if(lambda[i] > 0.99){
          lambda[i] = 0.99
        }
        Data[[i]]$Alpha = matrix(fit$a0,ncol=Data[[i]]$C)
        Data[[i]]$Beta = matrix(fit$beta,ncol=K*Data[[i]]$C)       
      }
    }
    
    newZ = mcmcMix(Data[[1]],Data[[2]],Data[[3]],Data[[4]],ndt,sdev,lastZ,n,K,n.burnin,n.draw)
    lastZ = newZ$lastZ
    
    for(i in 1:ndt){
      Alpha[[i]] = Data[[i]]$Alpha
      Beta[[i]] = Data[[i]]$Beta
    }
  }
  cat("\n")
  BIC = totalBIC(Data,newZ$meanZ,ndt,K)
  devRatio = dev.ratio(Data,newZ$meanZ,alpha,lambda,ndt,K)
  meanZ=newZ$meanZ
  rownames(meanZ) = rownames(xList[[1]])
  colnames(meanZ) = paste0("Factor",1:K)
  for(i in 1:ndt){
    rownames(Beta[[i]]) = colnames(xList[[i]])
    if(type[i] == "multinomial"){
      rownames(Alpha[[i]]) = colnames(xList[[i]])
      colnames(Alpha[[i]]) = paste0("C",1:(Data[[i]]$C))
      colnames(Beta[[i]]) = paste(rep(paste0("Factor",1:K), Data[[i]]$C),rep(paste0("C",1:K), each=Data[[i]]$C),sep="_")
    }else{
      names(Alpha[[i]]) = colnames(xList[[i]])
      colnames(Beta[[i]]) = paste0("Factor",1:K)      
    }
  }
  list(alpha=Alpha,beta=Beta,meanZ=meanZ,BIC=BIC,dev.ratio=devRatio,lambda=lambda)
}

