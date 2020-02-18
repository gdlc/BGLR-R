#### Ridge regression + Additive relationship matrix

In this example we use the wheat 599 dataset included in the package.
Here we include directly the markers as predictors.
We assign a gaussian prior to marker effects (model="BRR").
The within subject covariance matrix (t x t) is modeled as UNstructured, 
and we assign by default a Scaled-Inverse 
Chi-square with degree of freedom (scalar) df0, and scale (matrix, t x t) S0.
We also include a random effect to include information about the 
relationship between individuals using an additive relationship matrix
derived from pedigree. The within subject covariance matrix for this random 
effect is modelled as UNstructured covariance matrix, simular to the case 
of the within covariance matrix assigned to marker effects.
The variance covariance for the residuals is modelled also as UNstructured by default.

```R

library(BGLR)
data(wheat)
y<-wheat.Y
X<-wheat.X
X<-scale(X)/sqrt(ncol(X))
K<-wheat.A

ETA<-list(list(X=X,model="BRR"),list(K=K,model="RKHS"))
fm<-Multitrait(y=y,ETA=ETA,nIter=1000,burnIn=500)

#Residual covariance matrix
fm$resCov

#Genetic covariance matrix (markers)
fm$ETA[[1]]$Cov

#Marker effects
fm$ETA[[1]]$beta

#Genetic covariance matrix (relationship matrix derived from pedigree)
fm$ETA[[2]]$Cov

#random effects
fm$ETA[[2]]$u


```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
