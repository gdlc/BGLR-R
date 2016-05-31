The following example shows how to tell BGLR to save the samples of effects in binary files. 
  - The function ```readBinMat()``` can be used to read the samples.
  - The function ```getVariances()``` computes the sample-variance (```var()```) for sets of markers as well as the total variance.

 
```R
 library(BGLR)
 data(wheat)
  y=wheat.Y[,1] ; X=scale(wheat.X)
  dir.create('test_saveEffects')
  setwd('test_saveEffects')
  fm=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=12000,thin=2,burnIn=2000)
  B=readBinMat('ETA_1_b.bin')
  plot(B[,1],type='o',col=4)
  VAR=getVariances(B=B,X=X,sets=sample(1:20,size=1279,replace=T))
  head(VAR)
  plot(VAR[,"total"])
```
[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/inst/md/EXAMPLES.md)
