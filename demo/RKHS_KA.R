#This example will run a RKHS with kernel averaging

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
weights=NULL;
R2=0.5;

#Euclidean distance matrix
D=as.matrix(dist(wheat.X,method="euclidean"))
h=quantile(as.vector(D)^2,probs=.05)
K1=exp(-5/h*(as.matrix(D)^2))
K2=exp(-1/h*(as.matrix(D)^2))
K3=exp(-1/5/h*(as.matrix(D)^2))
df=5
S=as.numeric(var(y))/2*(df-2)

ETA=list(list(K=K1,model='RKHS',df0=df,S0=S),list(K=K2,model='RKHS',df0=df,S0=S),list(K=K3,model='RKHS',df0=df,S0=S))

  
fit_RKHS=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,
              saveAt=saveAt,df0=5,S0=NULL,weights=weights,R2=R2)
plot(fit_RKHS$yHat,y)
