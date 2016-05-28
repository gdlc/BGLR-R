

#### Parametric Random Regression with BGLR

The following examples illustrate the priors implemented so far for linear regression in BGLR. In thes examples we illustrate the use of these priors
one at a time.

**1. Flat Prior (FIXED)**

```R
 library(BGLR)
 data(mice); X=scale(mice.X); pheno=mice.pheno
 attach(pheno)
 fm=BGLR(y=Obesity.BMI,ETA=list( list(~GENDER+CoatColour+CageDensity,model='FIXED')), nIter=6000,burnIn=1000)
 fm2=lm(Obesity.BMI~GENDER+CoatColour+CageDensity)
 plot(cbind(c(fm$mu,fm$ETA[[1]]$b),coef(fm2))); abline(a=0,b=1)
```

**2. A simple simulation**

```R
 h2=.8
 QTL=seq(from=500,to=10000,by=500)
 nQTL=length(QTL); n=nrow(X)
 b=rep(1,nQTL)*sqrt(h2/nQTL)
 signal=X[,QTL]%*%b
 error=rnorm(n,sd=sqrt(1-h2))
 y=signal+error
```

**3. Gaussian Prior (BRR, RR-BLUP, BLUP)**

```R
 fmBRR=BGLR(y=y,ETA=list( list(X=X,model='BRR')), 
            nIter=nIter,burnIn=burnIn,saveAt='brr_')
 plot(fmBRR$ETA[[1]]$b,col=4,cex=.5, type='o')
```
**4. Scaled-t (BayesA)**


**5. Double-Exponential (Bayesian Lasso)**

**6. Point of mass at zero + Gaussian Slab (BayesC)**

**7. Point of mass at zero + t-Slab (BayesB)**

**8. Gaussian prior with set-specific variance (BRR_sets)**
