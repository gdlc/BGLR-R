#### Etimation of Genomic Heritability: two approaches

The following examples implement ideas described in *Genomic Heritabilit: What is it?* [(de los Campos et al., PLoS Genetics, 2015)](http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1005048) and references therein.

**Simple Simulation**
```R
 library(BGLR); data(wheat)
 X=scale(wheat.X)/sqrt(ncol(wheat.X))
 
 h20=.5
 nQTL=10; n=nrow(X)
 QTL=seq(from=100, to=1200,length=nQTL)
 
 b=rep(1,nQTL)
 signal=X[,QTL]%*%b
 signal=signal/sd(signal)*sqrt(h20)
 error=rnorm(n=n)
 error=error/sd(error)*sqrt(1-h20)
 y=signal+error
```

**Estimation using variance components**

```R
 fm=BGLR(y=y,ETA=list(list(X=X,model='BRR',saveEffects=T)),nIter=6000,burnIn=1000)
 varU=scan('ETA_1_varB.dat')
 varE=scan('varE.dat')
 h2=varU/(varU+varE)
 plot(h2,type='o',cex=.5,col=4);abline(h=c(h20,mean(h2[-c(1:200)])),lty=2,col=c(1,2),lwd=2)
```

**Estimation using the sample variance of genomic values at each iteration of the sampler.**

```R
 B=readBinMat('ETA_1_b.bin')
 h2_new=rep(NA,nrow(B))
 for(i in 1:length(h2_new)){
   u=X%*%B[i,]	
   h2_new[i]=var(u)
 }
  plot(h2_new,type='o',cex=.5,col=4);abline(h=c(h20,mean(h2_new)),lty=2,col=c(1,2),lwd=2)
```


```R
  d_h2=density(h2[-c(1:200)])
  d_h2_new=density(h2_new)
  plot(numeric()~numeric(),ylim=range(c(d_h2$y,d_h2_new$y)),xlim=range(c(d_h2$x,d_h2_new$x)) )
  lines(x=d_h2$x,y=d_h2$y,col=2)
  lines(x=d_h2_new$x,y=d_h2_new$y,col=4)
  h20=.5
  abline(v=h20) 
```

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
