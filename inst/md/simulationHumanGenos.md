###(1) Implementing shrinkage and variable selection methods with BGLR

In the following example we use human genotypes (SNPs from distantly related individuals) to simulate an additive trait and subsequently analyze the simulated data using four Bayesian models that differ on the prior distirbution of effects. The following figure display three prior densities commonly used in genomic models: the Gaussian density (black, used in the G-BLUP or Bayesian Ridge Regression model), the double-exponential density (red, used in the Bayesian Lasso), the scaled-t prior (green, used in model BayesA) and a point-of-mass-plus-slab prior (blue, used in models BayesB and BayesC).

<img src="https://github.com/gdlc/BGLR/blob/master/inst/md/priors.jpg" width="500">


#### Simulation
The following code simulates a simple trait with heritability 0.5 and 10 QTL.

```R
 load('Z.RData')
 n=nrow(Z) # number of individuals
 p=ncol(Z) # number of markers
 nQTL=10   # number of loci with non-null effects
 h2=.5     # trait heritablity
 QTLs=sample(1:p,size=nQTL) # position of the QTL
 effects=runif(min=.5,max=1.5,n=nQTL) # QTL effects
 signal=Z[,QTLs]%*%effects # genetic signal
 signal=scale(signal)*sqrt(h2)
 error=rnorm(n)*sqrt(1-h2)
 y=signal+error # simulated phenotype.
 
 # Generating a testing set
 tst=sample(1:n,size=300)
 yNA=y 
 yNA[tst]=NA # masking phenotypes in the testing set
##
```

#### Fitting linear regressions with BGLR
```R 
 library(BGLR)
 nIter=6000 # note, we use a limited number of iterations to illustrate; for formal analyses longer chains are needed.
 burnIn=1000

 # Gaussian prior (equivalent to Genomic BLUP, a shrinkage estimation method)
  fmBRR=BGLR(y=yNA,ETA=list(list(X=Z,model='BRR')), saveAt='brr_',nIter=nIter,burnIn=burnIn)

 # Point of mass at zero plus a slab (BayesB, variable selection and shrinkage)
  fmBB=BGLR(y=yNA,ETA=list(list(X=Z,model='BayesB')), saveAt='bb_',nIter=nIter,burnIn=burnIn)

 # Estimated effects
  plot(abs(fmBRR$ETA[[1]]$b),main='Model=BRR',ylab='|estimated effect|',cex=.5,col=2,type='o')
  abline(v=QTLs,lty=2,col=4)
  
  plot(abs(fmBB$ETA[[1]]$b),main='Model=BRR', ylab='|estimated effect|',cex=.5,col=2,type='o')
  abline(v=QTLs,lty=2,col=4)
 
 # Probability of inclusion
  plot(fmBB$ETA[[1]]$d,col=2,cex=.5);abline(v=QTLs,lty=2,col=4)
 
 # Prediction correlation in a testing set
  cor(y[tst],fmBRR$yHat[tst])
  
  cor(y[tst],fmBB$yHat[tst])

```

###(2) GBLUP model 

```R
 G=tcrossprod(scale(Z,scale=F))
 G=G/mean(diag(G))
 fmGBLUP=BGLR(y=y,ETA=list(list(K=G,model='RKHS')),saveAt='gblup_',nIter=nIter,burnIn=burnIn)
 fmGBLUP$ETA[[1]]$varU ; fmGBLUP$varE
```

[More Examples](https://doi.org/10.1534/genetics.114.164442)
