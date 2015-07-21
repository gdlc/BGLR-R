rm(list=ls())

setwd(tempdir())
 
data(wheat) 

n<-599   # should be <= 599
p<-300   # should be <= than 1279=ncol(X)
nQTL<-30 # should be <= than p
X<-wheat.X[1:n,1:p]

## Centering and standarization
for(i in 1:p)
{ 
  	X[,i]<-(X[,i]-mean(X[,i]))/sd(X[,i]) 
}
  
# Simulation
b0<-rep(0,p)
whichQTL<-sample(1:p,size=nQTL,replace=FALSE)
b0[whichQTL]<-rnorm(length(whichQTL),
                    sd=sqrt(1/length(whichQTL)))
signal<-as.vector(X%*%b0)
error<-rnorm(n=n,sd=sqrt(0.5))
y<-signal +error 


nIter=500;
burnIn=100;
thin=3;
saveAt='';
S0=NULL;
weights=NULL;
R2=0.5;
ETA<-list(list(X=X,model='BayesA'))
  
fit_BA=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn,thin=thin,saveAt=saveAt,df0=5,S0=S0,weights=weights,R2=R2)
plot(fit_BA$yHat,y)
