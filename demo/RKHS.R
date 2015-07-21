#This example will run a RKHS model

rm(list=ls())
setwd(tempdir())

data(wheat)
set.seed(12345)
varB=0.5*(1/sum(apply(X=wheat.X,MARGIN=2,FUN=var)))
b0=rnorm(n=1279,sd=sqrt(varB))
signal=wheat.X%*%b0
error=rnorm(599,sd=sqrt(0.5))
y=100+signal+error
 	
nIter=500;
burnIn=100;
thin=3;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;

K=wheat.X%*%t(wheat.X)

ETA=list(list(K=K,model='RKHS'))
  
fit_RKHS=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
plot(fit_RKHS$yHat,y)
