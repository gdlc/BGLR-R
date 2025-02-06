## Using Sparse Matrices in BGLR

The following examples demonstrate how to use sparse matrices in BGLR. 

To do this, 
 1) In the linear predictor `ETA` we allow for `ETA[[i]]$X` to be a sparse matrix (`dgCMatrix`) and
 2) Whave options for the model variable for sparse matrixes:



| Dense    | Sparse |
| -------- | ------- |
| BRR  | BRR_sparse  |
...

## Advantages of using sparse matrices

Using sparse matrices can help reducing memory requirement and can speed up computations. Because we use `dgCMatrix`, 
there will be sizable memory advantages if the sparsity of the incidence matrix is more than 50%. There will be 
also computational advantages which will increase with the sparsity of the matrix (see examples below).

We note, however, that the creation of the sparse matrix, if not adequatly thought, 
may have memmory bottlnecks (peaks in memory usage) that may reduce the mmemory 
advantage of using sparse matrices. We are working developing `apps` to 
efficiently create some types of sparse matrices.

## Example 1: Using a Sparse Matrix to store a Cholesky decomposition

In this example we fit a GBLUP model using the Cholesky decompositon of a GRM, 
first using dense, and then sparse matrices. Approximately 50% of the entries 
of a Cholesky decomposition are equal to zero. Thus, in this case the memory advantages 
are minimal (20%) and the computational advantages are modest 
(~28% reduction in computational time when the sparse matrix is used; 
however, note that the advantage may be higher for larger n).
  
```r
rm(list=ls())
library(BGLR)
library(Matrix)

data(wheat)
X<-scale(wheat.X)/sqrt(ncol(wheat.X))
y<-wheat.Y[,1]
G<-tcrossprod(X)
diag(G)<-diag(G)+1/1e4
L<-t(chol(G))
Ls<-as(L,"dgCMatrix")

object.size(Ls)/object.size(L)

# First using a dense matrix
ETA<-list(list(X=L,model="BRR"))

set.seed(123)
setwd(tempdir())
system.time(fm<-BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,verbose=FALSE))

# Now sparse
ETA<-list(list(X=Ls,model="BRR_sparse"))
set.seed(123)
system.time(fm_SP<-BGLR(y=y,ETA=ETA,nIter=12000,burnIn=2000,verbose=FALSE))

plot(fm$yHat,fm_SP$yHat)
abline(a=0,b=1,col=2)
```

## Example 2: Genotype by environment models (with 4 environments)

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

plot(fm$yHat,fm_SP$yHat)
abline(a=0,b=1)

```

## Example 3: Genotype by environment models using the Genomes 2 Fields data set (with 136 year-locations) 

```r

# Libraries
 library(MatrixModels)
 library(BGLR)

# Data
 DATA=read.csv('/mnt/research/quantgen/projects/G2F/data/PHENO.csv')

y=DATA$yield

Z.YL=model.matrix(~year_loc-1,data=DATA)
Z.YL_SP=model.Matrix(~year_loc-1,data=DATA,sparse=TRUE)

round(100*(1-object.size(Z.YL_SP)/object.size(Z.YL)),1) # ~93% less memmory


# Here, we use the `group` variable to estimate year-location specific error variances
# Note: using groups can increase computational time significantly
system.time(
   fm<-BGLR(y=y,ETA=list(list(X=Z.YL,model='BRR')),
		verbose=FALSE,group=DATA$year_loc)
)

system.time(
	fmSP<-BGLR(y=y,ETA=list(list(X=Z.YL_SP,model='BRR_sparse')),
			verbose=FALSE,group=DATA$year_loc)
)
plot(fm$varE,fmSP$varE);abline(a=0,b=1,col=2)
plot(fm$ETA[[1]]$b,fmSP$ETA[[1]]$b);abline(a=0,b=1,col=2)

```

### Example 4: GxE using marker-by-environment interactions using markers

The following examples illustrate how to implement marker-by-environments interaction models using BGLR, 
for further details about these models see [Lopez-Cruz et al. (2015)](https://doi.org/10.1534/g3.114.016097).

``` R
rm(list=ls())

library(BGLR)
library(Matrix)

data(wheat)
Y<-wheat.Y # grain yield evaluated in 4 different environments
round(cor(Y),2)

X<-scale(wheat.X)/sqrt(ncol(wheat.X))
 
y2<-Y[,2]
y3<-Y[,3]
y<-c(y2,y3)

X0<-matrix(nrow=nrow(X),ncol=ncol(X),0) # a matrix full of zeros

X_main<-rbind(X,X)
X_1<-rbind(X,X0)
X_2<-rbind(X0,X)

set.seed(123)
fm<-BGLR(y=y,
         ETA=list(main=list(X=X_main,model='BRR'),
                  int1=list(X=X_1,model='BRR'),
                  int2=list(X=X_2,model='BRR')),
	     nIter=6000,burnIn=1000,saveAt='GxE_',groups=rep(1:2,each=nrow(X))
	   )

varU_main<-scan('GxE_ETA_main_varB.dat')[-c(1:200)]
varU_int1<-scan('GxE_ETA_int1_varB.dat')[-c(1:200)]
varU_int2<-scan('GxE_ETA_int2_varB.dat')[-c(1:200)]
varE<-read.table('GxE_varE.dat',header=FALSE)[-c(1:200),]
varU1<-varU_main+varU_int1
varU2<-varU_main+varU_int2
h2_1<-varU1/(varU1+varE[,1])
h2_2<-varU2/(varU2+varE[,2])
COR<-varU_main/sqrt(varU1*varU2)
mean(h2_1)
mean(h2_2)
mean(COR)

unlink("*.dat")

#Now sparse
X_1s<-as(X_1,"dgCMatrix")
X_2s<-as(X_2,"dgCMatrix")

set.seed(123)
fms<-BGLR(y=y,
          ETA=list(main=list(X=X_main,model='BRR'),
                   int1=list(X=X_1s,model='BRR_sparse'),
                   int2=list(X=X_2s,model='BRR_sparse')),
	      nIter=6000,burnIn=1000,groups=rep(1:2,each=nrow(X))
	   )
	   
unlink("*.dat")

par(mfrow=c(2,2))

plot(fms$yHat,fm$yHat)

plot(fms$ETA[[1]]$b,fm$ETA[[1]]$b)
cor(fms$ETA[[1]]$b,fm$ETA[[1]]$b)

plot(fms$ETA[[2]]$b,fm$ETA[[2]]$b)
cor(fms$ETA[[2]]$b,fm$ETA[[2]]$b)

plot(fms$ETA[[3]]$b,fm$ETA[[3]]$b)
cor(fms$ETA[[3]]$b,fm$ETA[[3]]$b)

```

