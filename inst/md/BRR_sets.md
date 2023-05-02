#### Gausian Regression Grouping Markers into Sets: `BRR_sets`

The method `BRR_sets` allows users to group markers and their effects into sets, each of which has its own variance. Markers can be 
grouped in many different ways including LD-blocks, pathway information, or simply results from GWAS. In the following example we illustrate
how to use `BRR_sets`  with a two-step procedure wherein step 1 a GWAS is conducted and used to group markers and in a second step this grouping
is used to inform `BRR_sets`.

```R
 library(BGLR)
 data(mice); X=scale(mice.X[,1:3000]); pheno=mice.pheno
 attach(pheno)
  h2=.5
  QTL=floor(seq(from=50,to=ncol(X),length=20))
  nQTL=length(QTL); n=nrow(X)
  b=rep(1,nQTL)*sqrt(h2/nQTL)
  signal=X[,QTL]%*%b
  error=rnorm(n,sd=sqrt(1-h2))
  y=signal+error
  n=ncol(X)
  nIter=6000
  burnIn=1000
 detach(pheno)
```
 
**Standard GWAS**
```R
 library(BGData)
 PC=svd(x=X,nu=5,nv=0)$u
 TMP=GWAS(y~PC,data=new('BGData',geno=X,pheno=data.frame(y=y),map=data.frame()),method='lm')
 thresholds=c(.99,.95,.9,.75,.66,.5,.33)
 
 par(mfrow=c(2,1))
 x=-log10(TMP[,4])
 
 plot(x,cex=.5,col=8)
 points(x=QTL,col=2,cex=.7,pch=19,y=x[QTL])
```

**BGLR Using BRR_sets**
```R
 x= -log10(TMP[,4])
 tmp=sort(c(max(x)+.1,min(x)-.1,quantile(x,prob=thresholds)))

 groups=cut(x,breaks=tmp)
 fm=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=as.integer(groups))),
           nIter=nIter,burnIn=burnIn)
 x=abs(fm$ETA[[1]]$b)
 plot(x,col=8,cex=.7)
 points(x=QTL,y=x[QTL],col=2,cex=.7,pch=19)

```


**Comparison with BRR**
```R
 fmBRR=BGLR(y=y,ETA=list(list(X=X,model='BRR')),
           nIter=nIter,burnIn=burnIn)

 par(mfrow=c(2,1))
 x=abs(fm$ETA[[1]]$b)
 plot(x,col=8,cex=.7)
 points(x=QTL,y=x[QTL],col=2,cex=.7,pch=19)
 
 x=abs(fmBRR$ETA[[1]]$b)
 plot(x,col=8,cex=.7)
 points(x=QTL,y=x[QTL],col=2,cex=.7,pch=19)

```

**Comparison with BayesB**
```R
 fmBB=BGLR(y=y,ETA=list(list(X=X,model='BayesB')),
           nIter=nIter,burnIn=burnIn)

 par(mfrow=c(2,1))
 x=abs(fm$ETA[[1]]$b)
 plot(x,col=8,cex=.7)
 points(x=QTL,y=x[QTL],col=2,cex=.7,pch=19)
 
 x=abs(fmBB$ETA[[1]]$b)
 plot(x,col=8,cex=.7)
 points(x=QTL,y=x[QTL],col=2,cex=.7,pch=19)

```

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
