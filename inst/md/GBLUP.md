
### Various Ways of fitting a 'GBLUP' model using BGLR

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
 h2_2b=varU/(varU+varE)
```
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


**(5) Using the Cholesky decompositon and `model='BRR'`**
  This approach won't work if G is not positive definite. To avodi this problem we add a small constant in the diagonal.
  
```R
 diag(G)=diag(G)+1/1e4
 L=t(chol(G)) 
 
 fm5=BGLR( y=y,ETA=list(list(X=L,model='BRR')),nIter=nIter,
	   burnIn=burnIn,saveAt='chol_')
			
 varE=scan( 'chol_varE.dat')
 varU=scan('chol_ETA_1_varB.dat')
 h2_5=varU/(varU+varE)
```
[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
