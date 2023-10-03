# A Gibbs Sampler for a Bayesian Mixture Model
# C=X'X, rhs=X'y, my=mean(y), vy=var(y),n=length(y)
# Note: we don't include an intercept, so, before computing X'y, X'X, you should center X and y arount their means.
# B0, a matrix with prior estiamtes of marker effects, as many columns as prior estimates. We also suggest including one column
#     full of zeroes, to allow for just plain shrinkage towards zero, if needed (this is the default value)
# nIter, burnIn: the number of iterations and burnIn used for the Gibbs Sampler
# R2: the prior proportion of variance of y explained by the model (this is used to derive hyper-parameters for the prior variances)
# priorProb a vector of prior probabilities for the mixture compoentns (as many entries as columsn in B0), we internally scale it to add up to one
# priorCounts, the number of counts associated to priorProb, using priorCounts very large (e.g., 1e8) fixes the prior probabilities. The default value is 2*ncol(B0)
#  ...
##

BMM=function(C,rhs,my,vy,n,B0=matrix(nrow=ncol(C),ncol=1,0),nIter=150,burnIn=50,R2=.5,nComp=matrix(ncol(B0)),
                df0.E=5,S0.E=vy*R2*df0.E,df0.b=rep(5,nComp), priorProb=rep(1/nComp,nComp),priorCounts=rep(2*nComp,nComp),verbose=TRUE){

 # nIter=150;burnIn=50;R2=.5;nComp=matrix(ncol(B0));df0.E=5;S0.E=vy*R2*df0.E;df0.b=rep(5,nComp);alpha=.1;my=mean(y); vy=var(y); B0=cbind(rep(0,p),-1,1)
 p=ncol(C) 
 b=rep(0,p)
 d=rep(1,p) # indicator variable for the group
 POST.PROB=matrix(nrow=p,ncol=nComp,0)
	
 S0.b=df0.b*vy*R2/(sum(diag(C))/n)
 varB=S0.b/df0.b

 priorProb=priorProb/sum(priorProb)
 compProb=priorProb

 postMeanCompProb=rep(0,nComp)
	
 postMeanB=rep(0,p)
 postMeanVarB=rep(0,nComp)
 postProb=rep(0,nComp)

 varE=S0.E/df0.E 
 counts=priorCounts/nComp
	
 PROBS=matrix(nrow=p,ncol=nComp)
 	
 	
for(i in 1:nIter){
 
	 ## Future C code
	 for(j in 1:p){
	
		offset=sum(C[j,-j]*b[-j])
		rhs_j=(rhs[j]-offset)/varE
		rhs_j=rhs_j+B0[j,d[j]]/varB[d[j]]
 
		lhs_j=C[j,j]/varE+1/varB[d[j]]	
		sol=rhs_j/lhs_j
		b[j]=rnorm(n=1,mean=sol,sd=sqrt(1/lhs_j))
	 }
	 ## End of C-code
 
	 postMeanB=postMeanB+b/nIter
 
	 ## Sampling mixture components 
	 for(k in 1:nComp){
		 PROBS[,k]=dnorm(b,mean=B0[,k],sd=sqrt(varB[k]))#*compProb[k]	
	 }
 
	 d=apply(FUN=sample,x=1:nComp,X=PROBS,size=1,MARGIN=1,replace=TRUE)

	 ## Sampling the variance and the prior probabilities of the mixture components
	 for(k in 1:nComp){
		 tmp=(d==k)
	 
		 DF=sum(tmp)
		 SS=S0.b[k]
		 if(DF>0){
			 bStar=b[tmp]-B0[tmp,k]
			 SS=SS+sum(bStar^2)
		 }
		 
		 varB[k]=SS/rchisq(df=DF+df0.b[k],n=1) 
		 postMeanVarB[k]= postMeanVarB[k]+varB[k]/nIter

		 counts[k]=DF
		 
	 }

	compProb=rDirichlet(counts+priorCounts)
	postProb=postProb+compProb/nIter
 
	 ## computing posterior mean of mixture probabilities
	 for(k in 1:nComp){
		 tmp=(d==k)
		 POST.PROB[tmp,k]=POST.PROB[tmp,k]+1/nIter
	  }

	 # Sampling error variances
	 # we also need to add prior probabilities for the mixtures...
	 if(verbose){ print(i) }
  } 

   return(list(b=postMeanB,POST.PROB=POST.PROB,postMeanVarB=postMeanVarB,postProb=postProb))
}
## A function to sample from a Dirichlet

rDirichlet=function(counts){
    x=rgamma(n=length(counts),shape=counts)
	return(x/sum(x))
}



## Example

if(FALSE){

  rm(list=ls())

 library(BGLR)
 data(wheat)

 X=scale(wheat.X,center=TRUE,scale=FALSE)

 X=X[,1:100]

 p=ncol(X)
 n=nrow(X)

 QTL=seq(from=5,to=95,by=10)
 b0=rep(0,p)
 b0[QTL]=rep(c(-1,1),each=5)

 
 signal=X%*%b0
 error=rnorm(n=n,sd=sd(signal))
 y=signal+error

 varE=var(error)

 C=crossprod(X)
 rhs=crossprod(X,y)

 B0=cbind(rep(0,p),-1,1)
 
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/misc/mixturesWithNonZeroPriorMeans.R')
 tmp=BMM(C=C,rhs=rhs,my=mean(y),vy=var(y),B0=cbind(rep(0,p),-1,1),nIter=1500,burnIn=500)

 colQTL=ifelse(b0==0,8,ifelse(b0==1,2,4))
 pointQTL=ifelse(colQTL==8,1,19)
  par(mfrow=c(3,1))
  plot(tmp$POST.PROB[,1],col=colQTL,pch=pointQTL,ylim=c(0,1))
  plot(tmp$POST.PROB[,2],col=colQTL,pch=pointQTL,ylim=c(0,1))
  plot(tmp$POST.PROB[,3],col=colQTL,pch=pointQTL,ylim=c(0,1))

}

