

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
 
```
