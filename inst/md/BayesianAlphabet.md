

#### Parametric Random Regression with BGLR

The following examples illustrate the priors implemented so far for linear regression in BGLR. In thes examples we illustrate the use of these priors
one at a time.

**0. Flat Prior (FIXED)**

```R
 library(BGLR)
 data(mice); X=scale(mice.X); pheno=mice.pheno
 attach(pheno)
 fm=BGLR(y=Obesity.BMI,ETA=list( list(~GENDER+CoatColour+CageDensity,model='FIXED')), nIter=6000,burnIn=1000)
 fm2=lm(Obesity.BMI~GENDER+CoatColour+CageDensity)
 plot(cbind(c(fm$mu,fm$ETA[[1]]$b),coef(fm2))); abline(a=0,b=1)
```
**1. Gaussian Prior (BRR, RR-BLUP, BLUP)**


**2. Scaled-t (BayesA)**


**3. Double-Exponential (Bayesian Lasso)**

**4. Point of mass at zero + Gaussian Slab (BayesC)**

**5. Point of mass at zero + t-Slab (BayesB)**

**6. Gaussian prior with set-specific variance (BRR_sets)**
