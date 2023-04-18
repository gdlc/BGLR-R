## 1. Model

## 2. Example

**Data and data**

```r
 library(BGLR)
 data(wheat)
 X=wheat.X
 tmp=svd(X,nu=5)
 group=kmeans(tmp$u,centers=2,nstart=100)$cluster
 groups=unique(group)
 nGroups=length(groups)
```

**Centering within group**

```r
  for(i in 1:nGroups){
  tmp=group==groups[i]
  means=colMeans(X[tmp,])
  X[tmp,]=sweep(x=X[tmp,],STATS=means,FUN='-',MARGIN=2L)
 }
```

**G-matrix**

```r
 G=tcrossprod(X)
 G=G/mean(diag(G))
```

## 3. Sourcing the app

```r
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/apps/RKHS.Groups.R')
```

## 4. Fitting the model

```r
 fm=RKHS.Groups(y=wheat.Y[,1],K=G,group=group)
```

## 5. Extracting results

```r
 INTERCEPTS=c(fm$mu,fm$ETA$int$b)
 
 ## predictions
 PRED=data.frame(ID=rownames(X),group=group, y=y,yHat=yHat)
 head(PRED)
 
 ## Variance components
 fm$varE
 #fm$COV #to be done

```
