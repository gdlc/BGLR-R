

### Estimation of hyper-parameters

In BGLR the hyper-parameters entering in the prior distributions assigned to effects are treated as random; in all cases proper priors are
used. By default BGLR sets relatively flat priors; however the user can modify this. The following example illustrates this for two cases: BRR and 
BayesB. Further details about the prior used can be found in [Perez and de los Campos, Genetics (2014)](http://www.ncbi.nlm.nih.gov/pubmed/25009151).


**(1) Fitting GBLUP with fixed variance components**

BGLR assigns Scaled-inverse Chi-square priors to the variance parameters. In the parameterization used E[var]=S/(df-2) and Mode(var)=S/(df+2) where S>0 is a scale parameter and df>0 
is a degree of freedom parameter. To fix the variances, one can simply set `df` to a very large value and solve for the scale that sets
the prior mode equal to the desired value. The following example illustrates this.

```R
 library(BGLR)
 data(wheat)
 y=wheat.Y[,1]
 G=tcrossprod(scale(wheat.X))/ncol(wheat.X)
 
 # Fitting the model with  default priors (5df)
  fm1=BGLR(y=y,ETA=list(list(K=G,model='RKHS')),nIter=6000,burnIn=1000,saveAt='default_')
 
 # Fixing the variances
    h2=.8; DF=1e8
    Vy=var(y)
    Ve=Vy*(1-h2) 
    Vu=Vy*h2
    
    fm2=BGLR(y=y,ETA=list(list(K=G,model='RKHS',df0=DF,S0=Vu*(DF+2))),
             S0=Ve*(DF+2),df0=DF,
            nIter=6000,burnIn=1000,saveAt='h208_')
     plot(scan('h208_varE.dat'))
     plot(scan('h208_ETA_1_varU.dat'))
    
    c(fm1$varE,fm2$varE)
    c(fm1$ETA[[1]]$varU,fm2$ETA[[1]]$varU)
```

**(2) BayesB: estimating versus fixing the proportion of non-zero effects**

Models BayesB and BayesC use priors for effects which are finite-mixtures with a point of mass
at zero and a slab (Gaussian in case of BayesC, scaled-t in case of BayesB). One of the prior
hyper-parameters is the proportion of non-zero effects. BGLR assigns a beta prior to this parameter; the
Beta prior is index by two shape parameters (shape1 and shape2 in R) that can be thought as the number of prior
failures and successes in Bernoulli trials, we label them as counts0 and counts1, respectively. Alternatively you
can provide counts and probIn (success probability) and BGLR sets counts1=counts*probIn and counts0=counts*(1-probIn).
The mean of the prior distribution is given by counts1/(counts1+counts0)=probIn. Setting counts0=counts1=1 gives a uniform
prior.  As the number of counts approaches infinity the prior collapses to a point of mass at the prior mean. 
By default, BGLR sets counts0=counts1=5. The following example illustrates how to fit BayesB with flat prior and by fixing the
proportion of non-zero effects. 

```R
 rm(list=ls())
 library(BGLR)
  data(mice);X=scale(mice.X[,1:2000])
  h2=.5
  QTL=seq(from=50,to=ncol(X),length=20)
  nQTL=length(QTL); n=nrow(X)
  b=rep(1,nQTL)*sqrt(h2/nQTL)
  signal=X[,QTL]%*%b
  error=rnorm(n,sd=sqrt(1-h2))
  y=signal+error
 
 ## Using a flat prior
  fm1=BGLR(y=y,ETA=list(list(X=X,model='BayesB', counts0=1,counts1=1)),nIter=12000,burnIn=2000) 
   # alternatively, use above counts=2, probIn=.5
   # estimated effects (absolute value)
    plot(abs(fm1$ETA[[1]]$b));abline(v=QTL,col=4,lty=2) 

   # Estimated probablity of non-zero effects (average)
    fm1$ETA[[1]]$probIn

   # Estimated probablity of non-zero effects (per marker)    
    plot(fm1$ETA[[1]]$d);abline(v=QTL,col=4,lty=2)
   
   # Samples...
    TMP=read.table('ETA_1_parBayesB.dat',header=T)
    plot(TMP[,1],type='o',col=4)
    
 ## Using a extremely informative prior
   fm2=BGLR(y=y,ETA=list(list(X=X,model='BayesB', counts=1e5,probIn=1/1000)),nIter=12000,burnIn=2000)
  
   # estimated effects (absolute value)
    plot(abs(fm2$ETA[[1]]$b));abline(v=QTL,col=4,lty=2) 

   # Estimated probablity of non-zero effects (average)
    fm2$ETA[[1]]$probIn

   # Estimated probablity of non-zero effects (per marker)    
    plot(fm2$ETA[[1]]$d);abline(v=QTL,col=4,lty=2)
    
    #Suggestion: refit fm2 using probIn=.999
 
```


[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
