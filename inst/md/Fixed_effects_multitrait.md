#### Fixed effects in Multitrait models

```R

#Fixed effects
	
rm(list=ls())
library(BGLR)
data(wheat)
y<-wheat.Y[,1:3]
X<-wheat.X
X<-scale(X)/sqrt(ncol(X))

# Simulating a fixed effect
B <- rbind(c(1, 2, 3), c(-2, 1, 0))

y2<-matrix(NA,nrow=nrow(y),ncol=ncol(y))
y3<-matrix(NA,nrow=nrow(y),ncol=ncol(y))

XF4 <- matrix(rnorm(2*599), ncol = 2)
XF5 <- matrix(rexp(2*599), ncol = 2)
XF5 <- scale(XF5,center=TRUE,scale=FALSE)
XF6 <- matrix(rnorm(2*599), ncol = 2)
XF6 <- scale(XF6,center=TRUE,scale=FALSE)

B <- rbind(c(1, 2, 3), c(-2, 1, 0))
y2[,1]<-y[, 1] + XF4 %*% B[, 1]
y2[,2]<-y[, 2] + XF5 %*% B[, 2]
y2[,3]<-y[, 3] + XF6 %*% B[, 3]


#Covariates in XF4 affects trait 1
#Covariates in XF5 affects trait 2
#Covariates in FX6 affects trait 3 
ETA<-list(list(X=X,model="BRR"),
		   list(X=cbind(XF4,XF5,XF6),
		        model="FIXED",common=FALSE,
		        idColumns=c(1,1,2,2,3,3),
		        saveEffects=TRUE))

fm<-Multitrait(y=y2,ETA=ETA,nIter=10000)


str(fm$ETA[[2]])

#Compare against B
fm$ETA[[2]]$beta
B

Bsamples<-readBinMat('ETA_2_beta.bin')
colMeans(Bsamples)

#Another example, covariates only for traits 1 and 3
y3[,1]<-y[, 1] + XF4 %*% B[, 1]
y3[,2]<-y[, 2] 
y3[,3]<-y[, 3] + XF6 %*% B[, 3]

ETA2<-list(list(X=X,model="BRR"),
		   list(X=cbind(XF4,XF6),
		        model="FIXED",common=FALSE,
		        idColumns=c(1,1,3,3),
		        saveEffects=TRUE))

fm2<-Multitrait(y=y3,ETA=ETA2,nIter=10000)


str(fm2$ETA[[2]])

#Compare against B
fm2$ETA[[2]]$beta
B[,-2]

Bsamples<-readBinMat('ETA_2_beta.bin')
colMeans(Bsamples)


#X=[XF4,XF5,XF6]
#common effects, all covariates in X affects all the traits
#common=TRUE by default
ETA3<-list(list(X=X,model="BRR"),
		  list(X=cbind(XF4,XF5,XF6),
		       model="FIXED",common=TRUE))
		       
fm3<-Multitrait(y=y2,ETA=ETA3,nIter=10000)
str(fm3$ETA[[2]])
fm3$ETA[[2]]$beta
B

```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
