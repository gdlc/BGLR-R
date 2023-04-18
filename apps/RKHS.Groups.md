```r

## Wheat data set
 library(BGLR)
 data(wheat)
 X=wheat.X
 tmp=svd(X,nu=5)
 group=kmeans(tmp$u,centers=2,nstart=100)$cluster
 groups=unique(group)
 nGroups=length(groups)

 ## Centering within group
 for(i in 1:nGroups){
  tmp=group==groups[i]
  means=colMeans(X[tmp,])
  X[tmp,]=sweep(x=X[tmp,],STATS=means,FUN='-',MARGIN=2L)
 }
 
 G=tcrossprod(X)
 G=G/mean(diag(G))

## Sourcing the app
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/apps/RKHS.Groups.R')

 fm=RKHS.Groups(y=wheat.Y[,1],K=G,group=group)

```
