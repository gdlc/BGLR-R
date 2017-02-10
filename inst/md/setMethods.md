### Set Methods

In genomic regressions often one wants to consider to group SNPs (or other predictors) into sets and then assign set-specific priors.
This treatment allows set-specific rates of variable selection as well as set-specific shrinkage of estimates. For example, one may wish to
group SNPs based on bioinformatic categories, or based on allele frequencies or some other type of information. 

In `BGLR` models are specified using a list-interface, each element of the list is assigned a user-specified prior, this interface is very 
flexible and particularly convinient for implementing set methods. The following examples illustrate this. In the example results from
single-marker regressions are used to classify SNPs into sets.

**A simple simulation**
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

**Single Marker Regressions**

The following code performs an association analyses, the resulting z-scores are used to classify SNPs into sets.

```R
library(BGData)
 PC=svd(x=X,nu=5,nv=0)$u
 TMP=GWAS(y~PC,data=new('BGData',geno=X,pheno=data.frame(y=y),map=data.frame()),method='lm')
 thresholds=c(.99,.95,.9,.75,.66,.5,.33)

 par(mfrow=c(2,1))
 x=-log10(TMP[,4])

 plot(x,cex=.5,col=8)
 points(x=QTL,col=2,cex=.7,pch=19,y=x[QTL])
 
 x= abs(TMP[,3])
 tmp=sort(c(max(x)+.1,min(x)-.1,quantile(x,prob=thresholds)))

 groups=cut(x,breaks=tmp)
```

```R


```
