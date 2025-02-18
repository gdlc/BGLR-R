
## Various Ways of fitting a 'GBLUP' model using BGLR

In the following example we show how to fit a GBLUP model (i.e., a Gaussian process) using different parameterization.

<div id="menu" />
  
   * [Using oringial inputs (e.g., SNPs)](#BRR)
   * [Using a G-matrix (or kernel)](#RKHS)
   * [Using eigenvalues and eigenvectors](#RKHS2)
   * [Using scaled-principal components](#PC)
   * [Using a Cholesky decomposition](#CHOL)
   * [Using a QR decomposition](#QR)
   * [Using a Cholesky decomposition and sparse matrix](#CholSparse)
   

<div id="BRR" />

---------------------------------------------
**(i) Providing the markers, using `model='BRR'`**

In this case BGLR asigns iid normal priors to the marker effects.

``` R
 library(BGLR)
 data(wheat)

 nIter=12000
 burnIn=2000

 X=scale(wheat.X)/sqrt(ncol(wheat.X))
 y=wheat.Y[,1]

 fm1=BGLR( y=y,ETA=list(mrk=list(X=X,model='BRR')),
	   nIter=nIter,burnIn=burnIn,saveAt='brr_'
 	 )

 varE=scan('brr_varE.dat')
 varU=scan('brr_ETA_mrk_varB.dat')
 h2_1=varU/(varU+varE)
```
[Menu](#menu)


---------------------------------------------
<div id="RKHS" />

**(2) Providing the G-matrix**

BGLR Fits these Gaussian models using the eigenvalue decomposition og G. The eigenvalue decomposition is computed internally using 
`eigen()`.

```R
 G=tcrossprod(X)
 fm2=BGLR( y=y,ETA=list(G=list(K=G,model='RKHS')),
	   nIter=nIter,burnIn=burnIn,saveAt='eig_'
	 )
 varE=scan( 'eig_varE.dat')
 varU=scan('eig_ETA_G_varU.dat')
 h2_2=varU/(varU+varE)
```
[Menu](#menu)


<div id="RKHS2" />


---------------------------------------------
**(3) Providing eigenvalues and eigenvectors**

This strategy can be used to avoid computing the eigen-decomposition internally. This can be useful if a model will be fitted several times (e.g., cross-validation).

```R
 EVD=eigen(G)
 
 fm3=BGLR( y=y,ETA=list(G=list(V=EVD$vectors,d=EVD$values,model='RKHS')),
	    nIter=nIter,burnIn=burnIn,saveAt='eigb_')
 varE=scan( 'eigb_varE.dat')
 varU=scan('eigb_ETA_G_varU.dat')
 h2_3=varU/(varU+varE)
```
[Menu](#menu)


<div id="PC" />


---------------------------------------------
**(4) Providing scaled-eigenvectors and using `model='BRR'`**

```R
 PC=EVD$vectors
 for(i in 1:ncol(PC)){  PC[,i]=PC[,i]*sqrt(EVD$values[i]) }
 PC=PC[,EVD$values>1e-5]

 fm4=BGLR( y=y,ETA=list(pc=list(X=PC,model='BRR')),nIter=nIter,
	   burnIn=burnIn,saveAt='pc_')
			
 varE=scan( 'pc_varE.dat')
 varU=scan('pc_ETA_pc_varB.dat')
 h2_4=varU/(varU+varE)
```
[Menu](#menu)



<div id="CHOL" />


---------------------------------------------
**(5) Using the Cholesky decompositon and `model='BRR'`**
  This approach won't work if G is not positive definite; in our case the matrix is positive semi-definite, we can make it positive definite by adding a small constant to the diagonal.
  
```R
 diag(G)=diag(G)+1/1e4
 L=t(chol(G)) 
 
 fm5=BGLR( y=y,ETA=list(list(X=L,model='BRR',lower_tri=TRUE)),nIter=nIter, 
	   burnIn=burnIn,saveAt='chol_')
			
 varE=scan( 'chol_varE.dat')
 varU=scan('chol_ETA_1_varB.dat')
 h2_5=varU/(varU+varE)
```
[Menu](#menu)


<div id="QR" />


---------------------------------------------
**(6)Using QR-factorization**

```r
QR=qr(t(X))

Rt=t(qr.R(QR))

fm6=BGLR( y=y,ETA=list(list(X=Rt,model='BRR')),nIter=nIter,
	   burnIn=burnIn,saveAt='qr_')
			
 varE=scan( 'qr_varE.dat')
 varU=scan('qr_ETA_1_varB.dat')
 h2_6=varU/(varU+varE)

#Speeding up computations

fm7=BGLR( y=y,ETA=list(list(X=Rt,model='BRR',lower_tri=TRUE)),nIter=nIter,
           burnIn=burnIn,saveAt='qr2_')

 varE=scan( 'qr2_varE.dat')
 varU=scan('qr2_ETA_1_varB.dat')
 h2_7=varU/(varU+varE)


```
[Menu](#menu)

<div id="CholSparse" />


---------------------------------------------
**(7) Using the Cholesky decompositon and sparse matrix `model='BRR_sparse'`**
  This approach won't work if G is not positive definite; in our case the matrix is positive semi-definite, we can make it positive definite by adding a small constant to the diagonal. The resulting Cholesky factor can be represented as as sparse matrix using the 
  library Matrix.

  
```R
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

object.size(L)
object.size(Ls)

#Non zero values
Ls@x

#Cumulative number of non zero elements
Ls@p

#row index of each element
Ls@i

####################################
#WARNING: EXPERIMENTAL VERSION...
####################################

ETAs<-list(list(X=Ls,model="BRR_sparse"))

set.seed(123)
fms<-BGLR(y=y,ETA=ETAs,nIter=10000,burnIn=5000)
unlink("*.dat")

plot(y,fms$yHat)

ETA<-list(list(X=as.matrix(Ls),model="BRR"))

set.seed(123)
fm<-BGLR(y=y,ETA=ETA,nIter=10000,burnIn=5000)
unlink("*.dat")

plot(y,fm$yHat)

plot(fms$yHat,fm$yHat)
cor(fms$yHat,fm$yHat)

plot(fms$ETA[[1]]$b,fm$ETA[[1]]$b)
cor(fms$ETA[[1]]$b,fm$ETA[[1]]$b)
```
[Menu](#menu)


[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
