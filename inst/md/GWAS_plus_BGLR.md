#### GWAS followed by BGLR

In this example we show how to implement a two-step method consisting of GWAS followed by BGLR with variances of effects defined base on results obtained from GWAS.


**Simulation**
```R
 library(BGLR)
 data(mice)
 X=scale(mice.X)
 h2=.5 ; nQTL=20; n=nrow(X);p=ncol(X)
 QTL=seq(from=5,by=floor(p/nQTL),length=nQTL)
 b=rep(1,nQTL)*sqrt(h2/nQTL)
 signal=X[,QTL]%*%b
 error=rnorm(n,sd=sqrt(1-h2))
 y=signal+error
```
**GWAS-step**

```R
 library(BGData)
 TMP=GWAS(y~1,data=new('BGData',geno=X,pheno=data.frame(y=y),map=data.frame()),method='lm')
 head(TMP)
```

**BGLR-1**
Here we cluster markers based on estimated effects

```R
 absB=abs(TMP[,1])
 tmp=quantile(absB,prob=c(.25,.5,.85,.95,.99))
 groups=rep(1,p)
 for(i in 1:length(tmp)){
   groups=ifelse(absB>tmp[i],i+1,groups)
 }
 fm1=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=groups)),nIter=6000,burnIn=1000)
  plot(abs(fm1$ETA[[1]]$b),col=groups); abline(v=QTL,col='grey',lty=2)
  plot(-log10(TMP[,1]),col=groups); abline(v=QTL,col='grey',lty=2)  
```

**BGLR-2**

Clustering based on z-statistic.

```R
 absB=abs(TMP[,3])
 tmp=quantile(absB,prob=c(.25,.5,.85,.95,.99))
 groups=rep(1,p)
 for(i in 1:length(tmp)){
   groups=ifelse(absB>tmp[i],i+1,groups)
 }
 fm2=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=groups)),nIter=6000,burnIn=1000)
  plot(abs(fm2$ETA[[1]]$b),col=groups); abline(v=QTL,col='grey',lty=2)
  plot(-log10(TMP[,1]),col=groups); abline(v=QTL,col='grey',lty=2)  
```
*BGLR-3**
 logPValue=-log10(TMP[,4])
 tmp=quantile(logPValue,prob=c(.25,.5,.85,.95,.99))
 groups=rep(1,p)
 for(i in 1:length(tmp)){
   groups=ifelse(logPValue>tmp[i],i+1,groups)
 }
 fm3=BGLR(y=y,ETA=list(list(X=X,model='BRR_sets',sets=groups)),nIter=6000,burnIn=1000)
  plot(abs(fm3$ETA[[1]]$b),col=groups); abline(v=QTL,col='grey',lty=2)
  plot(-log10(TMP[,1]),col=groups); abline(v=QTL,col='grey',lty=2)
```


