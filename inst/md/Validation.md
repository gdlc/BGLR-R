
### Three ways of producing predictions in validation data using BGLR

**(0) Generating a testing set**

```R
  library(BGLR)
  data(wheat)
  X=scale(wheat.X)
  y=wheat.Y[,2]
  n=length(y)
  nTST=150
  tst=sample(1:n,size=nTST)
```


**(1) Setting the entries of the testing set in the response (y) to NA**

```R
  yNA=y
  yNA[tst]=NA
  fm=BGLR(y=yNA,ETA=list(list(X=X,model='BRR')),nIter=6000,burnIn=1000)
  yHat_1=fm$yHat[tst]
```


**(2) Splitting both X and y into training and testing sets**

```R
  yTRN=y[-tst]   ; yTST=y[tst]
  X.TRN=X[-tst,] ; X.TST=X[tst,]
  
  fm=BGLR(y=yTRN,ETA=list(list(X=X.TRN,model='BRR')),nIter=6000,burnIn=1000)
  yHat_2=fm$mu+as.vector(X.TST%*%fm$ETA[[1]]$b)

```

**(3) Using G-matrix (only valid for GBLUP)**

```R
  G=tcrossprod(X)/ncol(X)
  
  G11=G[-tst,-tst] # genomic relationships in the training data
  G21=G[tst,-tst]
  
  fm=BGLR(y=yTRN,ETA=list(list(K=G11,model='RKHS')),nIter=6000,burnIn=1000)
  yHat_3=fm$mu+as.vector(G21%*%solve(G11)%*%fm$ETA[[1]]$u)
  cor(cbind(yHat_1,yHat_2,yHat_3))

```

**(4) Cross-validation using a loop**

```R
 folds=sample(1:5,size=n,replace=T)
 yHatCV=rep(NA,n)
 
 timeIn=proc.time()[3]
  for(i in 1:max(folds)){
  	tst=which(folds==i)
  	yNA=y
    yNA[tst]=NA
    fm=BGLR(y=yNA,ETA=list(list(X=X,model='BRR')),nIter=6000,burnIn=1000)
    yHatCV[tst]=fm$yHat[tst]
  }
  
 proc.time()[3]-timeIn
```


**(5) Cross-validation in parallel at multiple cores**

```R
  # First create a function that fits one fold
  fitFold=function(y,W,folds,fold){
  	tst=which(folds==fold)
  	yNA=y
    yNA[tst]=NA
    fm=BGLR(y=yNA,ETA=list(list(X=W,model='BRR')),nIter=6000,burnIn=1000,verbose=F)
    tmp=fm$yHat[tst]
    names(tmp)=tst
    return(tmp)
  }
  
  library(parallel)
   
   system.time( tmp<-mclapply(FUN=fitFold,y=y,W=X,folds=folds,X=1:5,mc.cores=3))
   yHatCV2=unlist(tmp)
   yHatCV2=yHatCV2[order(as.integer(names(yHatCV2)))] # need to order the vector
   cor(yHatCV,yHatCV2)
   
```



[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
