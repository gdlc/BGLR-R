#### Fixed effects

Example adapted from  [MTM](http://quantgen.github.io/MTM/vignette.html) package.
Fixed effects can be added in any of the examples presented simply adding 
another term to the linear predictor and specifying the model as 
"FIXED". The parameter 'common' can be used to specify if all the effects are assumed
for all the traits.

```R

library(BGLR)
data(wheat)
y<-wheat.Y[,1:3]
K<-wheat.A

#Case 1: Simulating a fixed effect, assuming the same effects for all the traits
y1<-matrix(NA,nrow=nrow(y),ncol=ncol(y))
	
XF1 <- matrix(rnorm(2*599), ncol = 2)
B <- rbind(c(1, 2, 3), c(-2, 1, 0))
	
for (i in 1:3) {
  		y1[, i] <- y[, i] + XF1 %*% B[, i]
}

ETA1<-list(list(X=XF1,model="FIXED"),
		   list(K=K,model="RKHS"))

fm1<-Multitrait(y=y1,ETA=ETA1,nIter=1000,burnIn=500)

#Residual covariance matrix
fm1$resCov

#Fixed effects
fm1$ETA[[1]]

#Genetic covariance matrix
fm1$ETA[[2]]$Cov

#Random effects
fm1$ETA[[2]]$u

#Case 2: Simulating a fixed effect, with different incidence matrix for each trait
y2<-matrix(NA,nrow=nrow(y),ncol=ncol(y))
XF2 <- matrix(rnorm(2*599), ncol = 2)
XF3 <- matrix(rexp(2*599), ncol = 2)
XF3 <- scale(XF3,center=TRUE,scale=FALSE)
XF4 <- matrix(rnorm(2*599), ncol = 2)
XF4 <- scale(XF4,center=TRUE,scale=FALSE)
	
B <- rbind(c(1, 2, 3), c(-2, 1, 0))
y2[,1]<-y[, 1] + XF2 %*% B[, 1]
y2[,2]<-y[, 2] + XF3 %*% B[, 2]
y2[,3]<-y[, 3] + XF4 %*% B[, 3]
	           
ETA2<-list(list(X=cbind(XF2,XF3,XF4),model="FIXED",common=FALSE),
		   list(K=K,model="RKHS"))

fm2<-Multitrait(y=y2,ETA=ETA2,nIter=1000,burnIn=500)

#Fixed effects
fm2$ETA[[1]]

#Residual covariance matrix
fm2$resCov

#Genetic covariance matrix
fm2$ETA[[2]]$Cov

#Random effects
fm2$ETA[[2]]$u

```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
