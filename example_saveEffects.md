The following example shows how to tell BGLR to save the samples of effects in binary files. The function ```readBinMat()``` can be used to read the samples.

 
```R
 library(BGLR)
 data(wheat)
 y=wheat.Y[,1] ; X=scale(wheat.X)
 dir.create('test_saveEffects')
 setwd('test_saveEffects')
 fm=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=1200,thin=2,burnIn=200)
 B=readBinMat('ETA_1_b.bin')
 dim(B)
 plot(B[,1],type='o',col=4)
 
 setwd('../')
```
