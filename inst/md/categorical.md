
## Fitting models to binary ordinal traits (threshold model)

The examples below illustrate how to fit models for binary and ordinal traits using BGLR.

**Binary**

```R

  library(BGLR)
  data(wheat)
  y=wheat.Y[,1]
  X=scale(wheat.X)
  yBin<-ifelse(y>quantile(y,prob=.66),1,0)
  fm=BGLR(yBin,response_type='ordinal',ETA=list(list(X=X,model='BRR')),nIter=6000,burnIn=1000)

  fm$varE # fixed to 1 for identification 
  fm$mu   # same as above
  fm$threshold  
  head(fm$probs)
  boxplot(fm$probs[,2]~yBin)

```


**Ordinal**
```R
  yOrdinal=ifelse(y<quantile(y,prob=.33),1,ifelse(y<quantile(y,prob=.66),2,3))
  table(yOrdinal)
  fm=BGLR(yOrdinal,response_type='ordinal',ETA=list(list(X=X,model='BRR')),nIter=6000,burnIn=1000)

  fm$varE # fixed to 1 for identification
  fm$mu   # same as above
  fm$threshold
  head(fm$probs)
  
```
[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
