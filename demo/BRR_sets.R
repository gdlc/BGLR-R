#This example will run a standard Bayesian LASSO
rm(list=ls())
setwd(tempdir())

library(BGLR)

data(mice)
X=scale(mice.X,center=T,scale=F)
QTLs=seq(from=100,to=10000,length=10)
signal=rowSums(X[,QTLs])
signal=signal/sd(signal)
y=signal+rnorm(nrow(X))

mrkGroups=rep(1:2000,each=10)[1:ncol(X)]
fm=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=mrkGroups))) # Note: method and sets

plot(fm$ETA[[1]]$varB,cex=.1,col=4,type='o')
abline(v=QTLs,col=2,lty=2,lwd=.5)
