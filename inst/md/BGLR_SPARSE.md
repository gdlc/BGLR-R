## Using Sparse Matrices in BGLR

The following examples demonstrate how to use sparse matrices in BGLR. 

To do this, 
 1) In the linear predictor `ETA` we allow for ETA[[i]]$X` to be a sparse matrix (`dgCMatrix`) and
 2) Whave options for the model variable for sparse matrice:



| Dense    | Sparse |
| -------- | ------- |
| BRR  | BRR_sparse  |
...

## Advantages of using sparse matrices

Using sparse matrices can help reducing memmory requirement and can speed up computations. Because we use `dgCMatrix`, there will be sizable memmory advantages if the sparsity of the incidence matrix is more than 50%. There will be also computational advantages which will increase with the sparsity of the matrix (see examples below).

We note, however, that the creation of the sparse matrix, if not adequatly thought, may have memmory bottlnecks (peaks in mmmeory usage) that may reduce the mmemory advantage of using sparse matrices. We are working developing `apps` to efficiently create some types of sparse matrices.

## Example 1: Using a Sparse Matrix to store a Cholesky decomposition


## Example 2: Genotype by environment models

```r
rm(list=ls())
library(Matrix)
library(BGLR)
data(wheat)
 
nEnv=ncol(wheat.Y)
nGeno=nrow(wheat.Y)
y=as.vector(wheat.Y)
 
X0=scale(wheat.X,center=TRUE,scale=FALSE)
 
 
ETA=list()
# Main Effects
X=X0
for(i in 2:nEnv){
   X=rbind(X,X0)
}
ETA[[1]]=list(X=X,model='BRR')
 
# Interactions
for(i in 1:nEnv){
  X=matrix(nrow=nrow(X),ncol=ncol(X),0)
  ini=(i-1)*599+1
  end=ini+599-1
  X[ini:end,]=X0
  ETA[[i+1]]=list(X=X,model='BRR')
}
 
## Model without using sparse matrices
setwd('../output')
PWD=getwd()
 
setwd(tempdir())
  system.time(fm<-BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,verbose=FALSE))
 
## Now Sparse
ETA_SP=ETA
 
for(i in 2:(nEnv+1)){
        ETA_SP[[i]]$X=as(ETA[[i]]$X,"dgCMatrix")
        ETA_SP[[i]]$model="BRR_sparse"
}
 
object.size(ETA_SP)[1]/object.size(ETA)[1]
 
system.time(fm_SP<-BGLR(y=y,ETA=ETA_SP,nIter=12000,burnIn=2000,verbose=FALSE))

plot(fm$yHat,fm_SP$yHat);abline(a=0,b=1)

```


