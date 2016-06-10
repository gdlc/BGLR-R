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

**(2) Parallel computation in clusters**

In a High-performance computing cluster (HPCC) one can easily run multiple jobs in parallel at different nodes of the cluster. 
In high-dimensional models (e.g., hundereds of thousans of predictors) the burn-in period can be long. To overcome this problem we developed an experimental version (`BGLR2`) which allows users to:
  - save the internal environment to a file
  - run BGLR using a saved environment.

These tools can be used to efficiently collect large numbers of samples at different nodes of an HPCC.  The following example
illsutrates some of the features of `BGLR2`, relative to BGLR these are some of the additional arguments:
 - `saveEnv` (TRUE/FALSE) if TRUE a binary file containing a snapshot of the environment right at the end of the sampler is generated.
 - `BGLR_ENV` (character) a path and the name of a file containing a snapshot of a BGLR environment. If provided this environment is used to run a sampler. Since this environment contains all the elements of model specification (y, ETA, etc. ) these arguments do not need to be provided. However, a few arguments `saveAt`, `nIter`, `burnIn`,  `thin`, `rmExistingFiles`, are over-written by the call (that is, BGLR uses the values provided in the call and not the ones saved in the environment.
 - `newChain` (TRUE/FALSE) if FALSE the chain is continued (with the seed saved in `BGLR_ENV`) and samples are appended to already existing files. Otherwise, a new chain is generated. In this case the strarting values are the last ones collected in the run that generated `BGLR_ENV`, however new seed and new output files are generated.
 
**(1) Saving a snapshot of the environment at the end of the sampler.**

```R
rm(list=ls())
 dir.create('~/testBGLR2')
 setwd('~/testBGLR2')
 
 library(devtools)
 install_git('https://github.com/gdlc/BGLR-R')
 library(BGLR)
 data(wheat)
 X=wheat.X[,1:100]


 set.seed(1203)
 fm1a=BGLR2(y=wheat.Y[,1],ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),
        saveEnv=TRUE,saveAt='firstRun_',nIter=12000,burnIn=2000,thin=1)
 list.files()
```

**Let's now recover BGLR from sleep and run additional iterations**

```R
 fm1b=BGLR2(BGLR_ENV='firstRun_BGLR_ENV.RData',nIter=10000,thin=1,burnIn=0,newChain=FALSE)
 list.files()
 varE1=scan('firstRun_vare.dat')
 # Note that the number of samples in file are 22000=12000+1000
```

We can now check wheather the run in two-steps done above is equivalent to a single chain.

```R
 set.seed(1203)
 fm2=BGLR(y=wheat.Y[,1],ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),saveAt='secondRun_',nIter=22000,burnIn=2000,thin=1)
 c(fm1b$varE,fm2$varE)
 
 plot(fm1b$yHat,fm2$yHat)
 plot(scan('firstRun_mu.dat'),scan('secondRun_mu.dat'))
```

Note that you can also start a new chain (with different seed) using the saved environment. This will allow, for instance
running parallel chains all with starting values provided by a saved environment that may have been run for burn-in.


```R
 fm1c=BGLR2(BGLR_ENV='firstRun_BGLR_ENV.RData',nIter=10000,thin=1,burnIn=0,newChain=TRUE,saveAt='thirdRun_')
 list.files()
 varE1=scan('firstRun_varE.dat')[-c(1:12000)]
 varE3=scan('thirdRun_varE.dat')
 plot(varE1,varE3)
 
```



[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
