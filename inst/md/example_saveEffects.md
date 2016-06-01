The following example shows how to tell BGLR to save the samples of effects in binary files. The function ```readBinMat()``` can be used to read the samples.

```R
  library(BGLR); data(mice); X=scale(mice.X[,1:2000])
  h2=.5
  QTL=seq(from=50,to=ncol(X),length=20)
  nQTL=length(QTL); n=nrow(X)
  b=rep(1,nQTL)*sqrt(h2/nQTL)
  signal=X[,QTL]%*%b
  error=rnorm(n,sd=sqrt(1-h2))
  y=signal+error
  
  dir.create('test_saveEffects')
  setwd('test_saveEffects')
   fm=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=3000,thin=2,burnIn=1000)
   B=readBinMat('ETA_1_b.bin')
   plot(B[,1],type='o',col=4)
   plot(B[,QTL[1]],type='o',col=4)
   plot(B[,QTL[2]],type='o',col=4)
   plot(B[,QTL[3]],type='o',col=4)
  
```
[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
