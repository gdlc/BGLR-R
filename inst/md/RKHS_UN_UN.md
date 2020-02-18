#### Random effect model with unstructured covariance matrices

Example adapted from  [MTM](http://quantgen.github.io/MTM/vignette.html) package.
The first random effect will be unstructured and will be assigned a Scaled-Inverse 
Chi-square with degree of freedom (scalar) df0, and scale (matrix, t x t) S0. 
The following example illustrates how to fit a multiple-trait model using 
the wheat dataset included in the package for 599 wheat lines and 4 traits. 
In this example the covariance matrix of the random effect and that of model 
residuals are UNstructured by default.

```R

library(BGLR)
data(wheat)
K<-wheat.A
y<-wheat.Y

ETA<-list(list(K=K,model="RKHS"))
fm<-Multitrait(y=y,ETA=ETA,nIter=1000,burnIn=500)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Random effects
fm$ETA[[1]]$u

```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
