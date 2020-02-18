#### Recursive models

Example adapted from  [MTM](http://quantgen.github.io/MTM/vignette.html) package.
In recursive models covariances are modeled as regressions. For instance, 
a recursive model for the covariance matrix of the vector x can be 
x=Bx + delta where B is matrix of regression coefficients 
(typically lower-triangular) with zeros in the diagonal and delta
is a random vector (independent normal random variables in BGLR).
The covariance matrix of x is Cov(x,x')=(I-B)^{-1} * PSI * ((I-B)^{-1})', 
where PSI is the covariance matrix for delta, and we assume it is diagonal.
The parameters of the model include the recursive effects and the variance 
of delta. The first ones are assigned IID normal priors with null 
mean and variance var=100 and the variances are assigned IID 
scaled-inverse chi-squares priors.

```R

library(BGLR)
data(wheat)
K<-wheat.A
y<-wheat.Y

M <- matrix(nrow = 4, ncol = 4, FALSE)
M[3, 2] <- M[4, 2] <- TRUE # Adding recursion from trait 2 onto traits 3 and 4
M[4, 3] <- TRUE # Adding recursion from trait 3 on trait 4
	
ETA<-list(list(K=K,model="RKHS",Cov=list(type="REC",M=M)))

fm<-Multitrait(y=y,ETA=ETA,resCov=list(type="DIAG"), nIter=1000,burnIn=500)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Random effects
fm$ETA[[1]]$u

```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
