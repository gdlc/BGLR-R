### Imlementing shrinkage and variable selection methods with BGLR

In the following example we use human genotypes (SNPs from distantly related individuals) to simulate an additive trait and subsequently analyze the simulated data using four Bayesian models that differ on the prior distirbution of effects. The following figure display three prior densities commonly used in genomic models: the Gaussian density (black, used in the G-BLUP or Bayesian Ridge Regression model), the double-exponential density (used in the Bayesian Lasso) and a point-of-mass-plus-slab prior (used in models BayesB and BayesC).


![Priors](https://github.com/gdlc/BGLR/blob/master/priors.jpg)

The simulation uses real SNP genotypes, these can be downloaded from the following [link](https://www.dropbox.com/s/tkrnzipro28gah2/X_3K_30K.RData?dl=0).

#### Simulation
The following code simulates a simple trait with heritability 0.5 and 10 QTL.

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
 nIter=12000
 burnIn=2000
 # Gaussian prior (equivalent to Genomic BLUP)
  fmBRR=BGLR(y=y,ETA=list(list(X=X,model='BRR')), saveAt='brr_',nIeter=nIter,burnIn=burnIn)
 
 # t-prior (BayesA)
  fmBA=BGLR(y=y,ETA=list(list(X=X,model='BayesA')), saveAt='ba_',nIeter=nIter,burnIn=burnIn)
 
 # double-exponential prior (Bayesian Lasso)
 
  fmBL=BGLR(y=y,ETA=list(list(X=X,model='BL')), saveAt='bl_',nIeter=nIter,burnIn=burnIn)

 
 # Point of mass at zero plus a slab (BayesB)
  fmBB=BGLR(y=y,ETA=list(list(X=X,model='BayesB')), saveAt='bb_',nIeter=nIter,burnIn=burnIn)

 # Gausian pior with window-specific variances
 fmBRRW=BGLR(y=y,ETA=list(list(X=X,model='BRR_windows')),nIeter=nIter,burnIn=burnIn)

```

```

