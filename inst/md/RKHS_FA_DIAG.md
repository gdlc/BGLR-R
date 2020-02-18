#### Factor analysis

Example adapted from  [MTM](http://quantgen.github.io/MTM/vignette.html) package.
In Factor Analysis (FA), a covariance matrix is decomposed into common and 
specific factors according to O = BB' + PSI where B is a matrix of loadings 
(regressions of the original random effects into common factors) and PSI
is a diagonal matrix whose non-null entries give the variances of factors 
that are trait-specific. The loadings are assigned flat priors 
(normal priors with null mean and large variance) and the variances 
of the specific factors are assigned scaled-invers chi-squared with df 
and scale given by parameters df0 and S0 
(which if not given are assigned default values). The following example 
specifies a 1-common factor model. 
The matrix M in the example is a logical matrix of the same dimensions as B
with TRUE for loadings that the user want to estimate, 
and FALSE for those that should be zeroed out.

```R

library(BGLR)
data(wheat)
K<-wheat.A
y<-wheat.Y

M <- matrix(nrow = 4, ncol = 1, TRUE)
ETA<-list(list(K=K,model="RKHS",Cov=list(type="FA",M=M)))

fm<-Multitrait(y=y,ETA=ETA,resCov=list(type="DIAG"), nIter=1000,burnIn=500)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Random effects
fm$ETA[[1]]$u

```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
