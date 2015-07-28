### Comparison of shrinkage and variable selection methods

** Download the following genotype [data](https://www.dropbox.com/s/tkrnzipro28gah2/X_3K_30K.RData?dl=0)

```R
 load('X_3K_30K.RData')
 n=nrow(X)
 p=ncol(X)

## Simulation
 nQTL=10
 QTLs=sample(1:ncol(X),size=nQTL)
 effects=sample(min=.8,max=1.2,n=nQTL)
 signal=X[,QTLs]%*%effects
 signal=scale(signal)
 error=rnorm(n)
 y=signal+error
##

 
 library(BGLR)
 fm=BGLR(y=error+signal,ETA=list(list(X=X,model='BayesB')))

```

