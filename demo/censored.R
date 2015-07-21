rm(list=ls())
setwd(tempdir())

#loading libraries
library(survival)

#loading data included in BGLR
data(wheat) 

#simulation of data

X=wheat.X[,1:4]
n=nrow(X)
b=c(-2,2,-1,1)
error=rnorm(n)
y=X%*%b+ error
 
cen=sample(1:n,size=200)
yCen=y
yCen[cen]=NA
a=rep(NA,n)
b=rep(NA,n)
a[cen]=y[cen]-runif(min=0,max=1,n=200)
b[cen]=Inf

nIter=6000;
burnIn=1000;
thin=10;
saveAt='';
df0=5
S0=var(y)/2*(df0-2)
weights=NULL;
ETA=list(list(X=X,model='FIXED'))

fm1=BGLR(y=y,a=a,b=b,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,
         df0=df0,S0=S0,weights=weights)

#fits the model using survreg
event=ifelse(is.na(yCen),0,1)
time=ifelse(is.na(yCen),a,yCen)

surv.object=Surv(time=time,event=event,type='right')
fm2=survreg(surv.object~X, dist="gaussian")

plot(fm1$ETA[[1]]$b~fm2$coeff[-1],pch=19,col=2,cex=1.5, 
     xlab="survreg()", ylab="BGLR()")
abline(a=0,b=1,lty=2)
