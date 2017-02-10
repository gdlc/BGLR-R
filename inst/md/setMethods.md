### Set Methods

In genomic regressions one may want to group SNPs (or other predictors) into sets each of which is then assigned set-specific prior distributions of effects. This treatment allows set-specific rates of variable selection as well as set-specific shrinkage of estimates. Groupings can be based on prior information (e.g., bioinformatics categories), allele frequencies or some other type of information. 
In `BGLR` models are specified using a list-interface: each element of the list is assigned a user-specified prior. This interface is particularly convenient for implementing set methods. The following examples illustrate this. In the example results from single-marker regressions are used to group SNPs based on their z-scores.


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

 zSq=TMP[,3]^2
 tmp=sort(c(max(zSq)+.1,min(zSq)-.1,quantile(zSq,prob=thresholds)))
 groups=cut(zSq,breaks=tmp)
 plot(zSq,cex=.5,col=groups)
 
```
**Bayesian Regression with SNP sets (using `BayesB` as an example)**

```R
ETA=list()
groups=as.integer(groups)

for(i in 1:length(unique(groups))){
   ETA[[i]]=list(X=X[,groups==i],model='BayesB')
}

fm=BGLR(y=y,ETA=ETA,nIter=nIter,burnIn=burnIn)

# collecting estimates
bHat=numeric()
for(i in 1:length(ETA)){ bHat=c(bHat,fm$ETA[[i]]$b) }
bHat=bHat[match(colnames(X),names(bHat)]
plot(abs(bHat),col=4,cex=.5)
abline(v=QTL,lty=2,lwd=.5)

```
