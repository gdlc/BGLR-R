### Parallel Computations with BGLR

**(1) Multi-core computing using the parallel package**

One can easily generate prallel chains at multiple cores using the ``parallel`` package. This can be conviniently done generating a wrapper to 
``BGLR``.

```R
 BGLR.wrap=function(task,seeds,...){
    seed=seeds[task]
    set.seed(seed)
    fm=BGLR(saveAt=paste0(task,'_'),...)
    return(list(fm=fm,task=task,seed=seed))
 }

```

Now we can run multiple chains in parallel.

```R
  library(parallel)
  library(BGLR)
  seeds=c(100,110,120)
  data(wheat)
  X=scale(wheat.X)
  y=wheat.Y[,1]
  ETA=list(list(X=X,model='BRR'))
  fmList=mclapply(FUN=BGLR.wrap,seeds=seeds,X=1:3,mc.cores=2,nIter=6000,burnIn=1000,verbose=F,y=y,ETA=ETA)
  fmList[[1]]$fm$varE
  
```

Using a similar approach we can conduct a [cross-validation in parallel](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Validation.md).
