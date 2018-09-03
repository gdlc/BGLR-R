### Testing individual SNP sets

Toy simualtion

```r
  library(BGLR)
  data(mice)
  X=scale(mice.X[,1:1000])
  QTL=seq(from=10,to=990,by=100)
  p=ncol(X);q=length(QTL)
  b=rep(0,p);b[QTL]=runif(min=.5,max=1,n=q)
  signal=X%*%b
  error=rnorm(sd(signal))
  y=error+signal

```

## BayesA

This model has no point of mass at zero, but we can caluclate the posterior probability that the effect is positive, and
the posterior probability that is negative. If the posterior distribution is away from zero, one of the two will be close to one and the other will be close to zero. 
A posterior distribution providing no evidence that there is an effect will be symmetric around zero, in this case p(bj>0|data)=p(bj<0|data)=0.5.
```r
 fmBA=BGLR(y=y,ETA=list(list(X=X,model='BayesA',saveEffects=TRUE)),nIter=12000,burnIn=2000,saveAt='BA_')
 BA=readBinMat('BA_ETA_1_b.bin')
 
 ## One sided posterior probability of effects being different than zero
  pNeg=colMeans(BA<0)
  pPos=colMeans(BA>0)
  pBA=ifelse(pNeg>pPos,pNeg,pPos)
  plot(pBA);abline(v=QTL,col=2,lty=2) 
```
## BayesB

```r
 fmBB=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=12000,burnIn=2000,saveAt='BB_')
 BB=readBinMat('BB_ETA_1_b.bin')
 
 ## One sided posterior probability of effects being different than zero
  pNeg=colMeans(BB<0)
  pPos=colMeans(BB>0)
  pBB=ifelse(pNeg>pPos,pNeg,pPos)
  plot(pBB);abline(v=QTL,col=2,lty=2) 
```
