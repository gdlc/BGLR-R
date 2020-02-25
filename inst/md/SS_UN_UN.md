#### Spike slab

Fitting selection variable models with multi trait models with simulated data.

```R

library(BGLR)
data(simulated3t)

y<-as.matrix(simulated3t.pheno[,1:3])
g<-as.matrix(simulated3t.pheno[,4:6])
cov(g)
y<-scale(y,center=TRUE,scale=FALSE)
y.orig<-y
	
X<-simulated3t.X
X<-scale(X)/sqrt(ncol(X))

ETA1<-list(list(X=X,model="SpikeSlab",
		        inclusionProb=list(probIn=rep(1/100,ncol(y)),
		        counts=rep(1E6,ncol(y)))))

#Fit the model

fm1<-Multitrait(y=y,ETA=ETA1,nIter=1000,burnIn=500)

#Residual covariance, UN
fm1$resCov

#Compare against the TRUE residual covariance matrix
#6.0 6.0 1.0
#6.0 8.0 2.0  
#1.0 2.0 1.0

#Genetic co-variance
crossprod(fm1$ETA[[1]]$beta)/fm1$ETA[[1]]$p

#Compare against the TRUE genetic co-variance matrix
#1.00   0.34   0.07
#0.34   1.00   0.21 
#0.07   0.21   1.00 
	
#Covariance matrix for b, UN
fm1$ETA[[1]]$Cov
	           
```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
