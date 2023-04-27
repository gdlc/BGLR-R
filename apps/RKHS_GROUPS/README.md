## Background

This app fits GBLUP models with interactions between SNPs and groups (envitonments, sex, race, etc).

The following features are worth mentioning:

  - Separate intercepts are modeled for each group.
  - The model includes main effects (shared across groups) and group-specific interactions.  
  - The genomic variance for each group is the sum of the main effect variance plus the group-specific interaction variance.
  - Different error variances are modeled for each group.
  - Errors are assumed to be un-correlated wihtin and between groups.


## Example 1: SNP-by-ancestry group

In this example we identigy two clusters (of wheat inbred lines) and model SNP-by-group interactions
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

**Sourcing the app**

```r
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/apps/RKHS_GROUPS/RKHS_GROUPS.R')
```

**Fitting the model**

```r
 fm=RKHS.Groups(y=wheat.Y[,1],K=G,group=group)
```

**Extracting results**

```r
 INTERCEPTS=c(fm$mu,fm$ETA$int$b)
 
 ## predictions
 PRED=data.frame(ID=rownames(X),group=group, y=y,yHat=yHat)
 head(PRED)
 
 ## Variance components
 fm$varE
 #fm$COV #to be done

```
