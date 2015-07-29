### Comparison of shrinkage and variable selection methods

In the following example we use human genotypes (SNPs from distantly related individuals) to simulate an additive trait and subsequently analyze the simulated data using four Bayesian models that differ on the prior distirbution of effects. 

The genotypes used for the simulation can be downloaded from the following[link](https://www.dropbox.com/s/tkrnzipro28gah2/X_3K_30K.RData?dl=0).

#### Simulation
```R
 load('~/Dropbox/shortCourseUAB/X_3x10.RData')
 n=nrow(X) # number of individuals
 p=ncol(X) # number of markers
 nQTL=10   # number of loci with non-null effects
 h2=.5     # trait heritablity
 QTLs=sample(1:ncol(X),size=nQTL) # position of the QTL
 effects=runif(min=.8,max=1.2,n=nQTL) # QTL effects
 signal=X[,QTLs]%*%effects # genetic signal
 signal=scale(signal)*sqrt(h2)
 error=rnorm(n)*sqrt(1-h2)
 y=signal+error # simulated phenotype.
 
 # Generating a testing set
 tst=sample(1:n,size=500)
 yNA=y 
 yNA[tst]=NA # masking phenotypes in the testing set
##
```

#### Fitting four different models with BGLR
```R
 G=tcrossprod(scale(X))/ncol(X)
 EVD=eigen(G)
 plot(y=c(0,cumsum(EVD$values)/sum(EVD$values)),x=0:nrow(G),type='o',col=2,cex=.5)
 plot(EVD$vectors[,1:2],cex=.5,col=2)
 library(BGLR)
 
 fmGBLUP=BGLR(y=y,ETA=list(list(V=EVD$vectors,d=EVD$values,model='RKHS')))
 fmBRR=BGLR(y=y,ETA=list(list(X=X,model='BRR')), saveAt='brr_')
 

 
 fmGBLUP=BGLR(y=error+signal,ETA=list(list(X=X,model='BRR')))
 fmBB=BGLR(y=error+signal,ETA=list(list(X=X,model='BayesB')))

```

```

