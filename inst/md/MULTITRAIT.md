## Multi-trait prediction using SVD and univariate regressions on eigenvectors


### Deriving Eigenvectors of the phenotype matrix
```r
 library(BGLR)
 data(wheat)
 
 # A phenotype matrix with 4 traits
  Y=wheat.Y 

 # Molecular markers
  X=scale(wheat.X) 
 
 # Singular value decompositions Y=UDV'
  SVD=svd(Y)
  U=SVD$u
  D=diag(SVD$d)
  V=SVD$v
```

### Obtaining regression coefficients

```r
 B=matrix(nrow=ncol(X),ncol=ncol(Y))

 ETA=list(list(X=X,model='BayesB'))
 for(i in 1:ncol(Y)){
	fm=BGLR(y=U[,i],ETA=ETA,verbose=F) #use more iterations!
	B[,i]=fm$ETA[[1]]$b
 }

 # Rotating coefficients to put them in marker space
  BETA=B%*%D%*%t(SVD$v)

 # Prediction
  YHat=X%*%BETA
  
  # correlation in training
  for(i in 1:ncol(Y)){ print(cor(Y[,i],YHat[,i])) }
```


### Training-testing comparing single and multi-trait


*Training/Testing partition*

```r
  set.seed(195021)
  N=nrow(X)
  tst=sample(1:N,size=300)
  Y.TRN=Y[-tst,]
  Y.TST=Y[tst,]
  
  X.TST=X[tst,]
  
  ETA=list(list(X=X[-tst,],model='BayesB'))
```

*Prediction using single-trait models*

```r
  YHAT.ST=matrix(nrow=nrow(Y.TST),ncol=4)
	
  for(i in 1:4){
    fm=BGLR(y=Y.TRN[,i],ETA=ETA,verbose=F) #use more iterations!
    YHAT.ST[,i]=X.TST%*%fm$ETA[[1]]$b
  }	
```


*Prediction using eigenvectors*

```r
  SVD=svd(Y.TRN)
  U=SVD$u
  D=diag(SVD$d)
  V=SVD$v

  B=matrix(nrow=ncol(X), ncol=ncol(U))
	
  for(i in 1:4){
     fm=BGLR(y=U[,i],ETA=ETA,verbose=F) #use more iterations!
     B[,i]=fm$ETA[[1]]$b
  }	
  BETA=B%*%D%*%t(V)
  YHAT.EV=X.TST%*%BETA
```

*Comparison*

```r
  for(i in 1:4){
     message( round(cor(Y.TST[,i],YHAT.ST[,i]),3),"   ",round(cor(Y.TST[,i],YHAT.EV[,i]),3))
  }
```

[Return to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
