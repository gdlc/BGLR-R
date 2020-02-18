#### Ridge regression

In this example we use the wheat 599 dataset included in the package.
Here we include directly the markers as predictors.
We assign a gaussian prior to marker effects (model="BRR").
The within subject covariance matrix (t x t) is modeled as UNstructured, 
and we assign by default a Scaled-Inverse 
Chi-square with degree of freedom (scalar) df0, and scale (matrix, t x t) S0. 
The variance covariance for the residuals is modelled also as UNstructured by default.
Other options that can be used for modelling the within subject 
covariance matrix are "DIAG", "FA" and "REC", likewise for the error.

```R

library(BGLR)
data(wheat)
y<-wheat.Y
X<-wheat.X
X<-scale(X)/sqrt(ncol(X))

ETA<-list(list(X=X,model="BRR"))
fm<-Multitrait(y=y,ETA=ETA,nIter=1000,burnIn=500)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix
fm$ETA[[1]]$Cov

#Marker effects
fm$ETA[[1]]$beta


```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
