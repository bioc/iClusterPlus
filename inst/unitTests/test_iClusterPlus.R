###### Code used to test the iClusterPlus functions through simulations ########

library(iClusterPlus)
library(RUnit)

n1=n2=n3=n4=20  # assume there are 4 clusters, each cluster has 20 samples
n=n1+n2+n3+n4   # total number of samples 
p=200           # each sample has 200 genomic features (variables) 
p1=p2=p3=10     # For cluster 1,2,3, among the 200 genomic features, 10 of them are informative genomic features that contribute to clustering the samples
                # and the rest of them are noises, which do not contribute to the clustering.

true.clusters=c(rep("C1",n1),rep("C2",n2),rep("C3",n3),rep("C4",n4))

set.seed(339)

### For each of the following data type, the rows and columns represent samples and features, respectively.

################# normal distribution ####################

### Assume the normal data with the following back ground distribution ###
normData = matrix(rnorm(n*p, mean=0, sd=0.5),nrow=n, ncol=p) 

### For the samples in each cluster, assume the informative features with the following distributions ###
### note: if the informative features are similar to background, their coefficients may be shrunken to zero

c1=matrix(rnorm(n1*p1,mean=3,sd=0.5),ncol=p1,nrow=n1)
c2=matrix(rnorm(n2*p2,mean=2,sd=0.5),ncol=p2,nrow=n2)
c3=matrix(rnorm(n3*p3,mean=1,sd=0.5),ncol=p3,nrow=n3)

normData[1:n1,1:p1]= c1                                   #normData[1:n1,] are the data for cluster 1, 1:p1 are informative features
normData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2               #normData[(n1+1):(n1+n2),] are the data for cluster 2, (p1+1):(p1+p2) are informative features
normData[(n1+n2+1):(n1+n2+n3),(p1+p2+1):(p1+p2+p3)]= c3   #normData[(n1+n2+1):(n1+n2+n3),] are the data for cluster 3, (p1+p2+1):(p1+p2+p3) are informative features
                                                          #Note: by default, the data for cluster 4 are from the background distribution.

### For the following data types, we use similar simulation approach as the normal case. ###

######## binomial data #########################

### Assume the binomial data with the following background distribution
binomData = matrix(rbinom(n*p, size=1, prob=0.05),nrow=n, ncol=p)

### For the samples in each cluster, assume the informative features with the following distributions ###
c1=matrix(rbinom(n1*p1,size=1, prob=0.8), ncol=p1,nrow=n1)
c2=matrix(rbinom(n2*p2,size=1, prob=0.5), ncol=p2,nrow=n2)
c3=matrix(rbinom(n3*p3,size=1, prob=0.3), ncol=p3,nrow=n3)

binomData[1:n1,1:p1]= c1
binomData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
binomData[(n1+n2+1):(n1+n2+n3),(p1+p2+1):(p1+p2+p3)]= c3
#Note: by default, the data for cluster 4 are from the background distribution. 

### some features (variables) may result in only one category.     ####
### in that case, we need to remove these non-informative variables   ####
ncat = function(x){length(unique(x))}
len = apply(binomData,2,ncat)
unicat = which(len==1)
if(length(unicat)>0){
  binomData = binomData[,-unicat]
}
dim(binomData)

######### poisson data ############################
poisData = matrix(rpois(n*p, lambda=0.5),nrow=n, ncol=p)

c1=matrix(rpois(n1*p1,lambda=5), ncol=p1,nrow=n1)
c2=matrix(rpois(n2*p2,lambda=3), ncol=p2,nrow=n2)
c3=matrix(rpois(n3*p3,lambda=2), ncol=p3,nrow=n3)

poisData[1:n1,1:p1]= c1
poisData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
poisData[(n1+n2+1):(n1+n2+n3),(p1+p2+1):(p1+p2+p3)]= c3
#Note: by default, the data for cluster 4 are from the background distribution.

################## multinomial distribution ##############
### we assume the class 2 is the normal case and class 1 or 3 are rare. #########
multData = matrix(sample(1:3,size=n*p, replace=TRUE, prob=c(0.05,0.9,0.05)),nrow=n, ncol=p)

c1=matrix(sample(1:3,size=n1*p1, replace=TRUE, prob=c(0.8,0.15,0.05)), ncol=p1,nrow=n1)
c2=matrix(sample(1:3,size=n2*p2, replace=TRUE, prob=c(0.3,0.4,0.3)), ncol=p2,nrow=n2)
c3=matrix(sample(1:3,size=n3*p3, replace=TRUE, prob=c(0.05,0.15,0.8)), ncol=p3,nrow=n3)

multData[1:n1,1:p1]= c1
multData[(n1+1):(n1+n2),(p1+1):(p1+p2)]= c2
multData[(n1+n2+1):(n1+n2+n3),(p1+p2+1):(p1+p2+p3)]= c3
#Note: by default, the data for cluster 4 are from the background distribution. 

### function to plot the coefficient of the genomic features ###
plotBeta = function(res){
  nz = ncol(res)
  nf = nrow(res)
  maxbeta = max(res)
  minbeta = min(res)
  par(mfrow=c(nz,1),mar=c(4.5,4.5,1,1))
  for(i in 1:nz){
    plot(res[,i],ylim=range(res),xlab="Genomic feature",ylab=paste("beta",i),col=c(rep("red",30),rep("blue",nf-30)))
    abline(h=0,col="grey")
  }
}

### the lambda needs to be tuned in order to get the 'best' result.  The tune.iClusterPlus function can be used for this purpose ###
### here, the parameter lambda used is not the best.  We just want to get a sense of what the results look like ###

#norm.res = iClusterPlus(dt1=normData,type=c("gaussian"),K=3,alpha=1,lambda=0.03)
#table(norm.res$clusters,true.clusters)
norm.res = iClusterPlus2(xList=list(normData),type=c("gaussian"),K=3,BICrate.cutoff=0.005)
km4 = kmeans(norm.res$meanZ,4, nstart = 100)
table(km4$cluster,true.clusters)
plotBeta(norm.res$beta[[1]])

#binom.res = iClusterPlus(dt1=binomData,type=c("binomial"),K=3,alpha=1,lambda=0.05)
#table(binom.res$clusters,true.clusters)
binom.res = iClusterPlus2(xList=list(binomData),type=c("binomial"),K=3,BICrate.cutoff=0.005)
km4 = kmeans(binom.res$meanZ,4, nstart = 100)
plotBeta(binom.res$beta[[1]])

pois.res = iClusterPlus2(xList=list(poisData),type=c("poisson"),K=3,BICrate.cutoff=0.005)
km4 = kmeans(pois.res$meanZ,4, nstart = 100)
table(km4$cluster,true.clusters)
plotBeta(pois.res$beta[[1]])

#mult.res = iClusterPlus(dt1=multData,type=c("multinomial"),K=3,alpha=1,lambda=0.02)
#table(mult.res$clusters,true.clusters)
mult.res = iClusterPlus2(xList=list(multData),type=c("multinomial"),K=3)
km4 = kmeans(mult.res$meanZ,4, nstart = 100)
table(km4$cluster,true.clusters)

plotBeta(mult.res$beta[[1]][,1:3])   # coefficients for the genomic features of class 1
plotBeta(mult.res$beta[[1]][,4:6])   # coefficients for the genomic features of class 2 
plotBeta(mult.res$beta[[1]][,7:9])   # coefficients for the genomic features of class 3 

icfit = iClusterPlus2(xList=list(normData,binomData,poisData,multData), 
                        type=c("gaussian","binomial","poisson","multinomial"),K=3)

### The following are the plots for the coefficients of the genomic features for each data type
### Each coefficient (row) is a vector with k elements;
### For each coefficient vector, if any of the elements is not close to zero, it's corresponding feature is informative.
### For each coefficient vector, if all the elements are close to zero, it's corresponding feature is non-informative.

plotBeta(icfit$beta[[1]])   # coefficients for the genomic features for the normal data  
plotBeta(icfit$beta[[2]])   # coefficients for the genomic features for the binomial data 
plotBeta(icfit$beta[[3]])   # coefficients for the genomic features for the Poisson data 

plotBeta(icfit$beta[[4]][,1:3])  # coefficients for the genomic features of class 1  
plotBeta(icfit$beta[[4]][,4:6])  # coefficients for the genomic features of class 2
plotBeta(icfit$beta[[4]][,7:9])  # coefficients for the genomic features of class 3

### The following code is used to get the top 30 most significant coefficients.
### By design, the first 30 features are informative.
### As expected, for most of the data type, the first 30 coefficients are the most significant.
maxAbsVal.1 = apply(abs(icfit$beta[[1]]),1,max) #get the max absolute value for each coefficient (row)
sort(order(maxAbsVal.1,decreasing=TRUE)[1:30])
      
maxAbsVal.2 = apply(abs(icfit$beta[[2]]),1,max)
sort(order(maxAbsVal.2,decreasing=TRUE)[1:30])

maxAbsVal.3 = apply(abs(icfit$beta[[3]]),1,max)
sort(order(maxAbsVal.3,decreasing=TRUE)[1:30])

maxAbsVal.4 = apply(abs(icfit$beta[[4]]),1,max)
sort(order(maxAbsVal.4,decreasing=TRUE)[1:30])

### check if all the samples are correctly clustered ###
km4 = kmeans(icfit$meanZ,4, nstart = 100)
#um = umap(icfit$meanZ)
#plot(um$layout,col=km4$cluster)

clusterTab = table(km4$cluster,true.clusters) # all samples are correctly assigned in different groups 
clusterTab
clusterRowSum = apply(clusterTab,1,sum)
trueNum = rep(20,4)  #For the clusters 1-4, each of them should have 20 samples
names(trueNum) = 1:4
#checkEquals(clusterRowSum,trueNum)
