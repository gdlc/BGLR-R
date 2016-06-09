### Censored Regression

BGLR supports right, left and interval censoring. For censored outcomes the response is represented with a triplet ai,yi,bi, where ai and bi are the lower and upper bounds for the phenotype (yi). The following table describes the configuration of the tripled for un-censored, right, left and interval censored data

Case               | ai  |  yi  | bi
-------------------|-----|------|----
 Not-censored      | NA  |  yi  | NA
 Right-Censored    | ai  |  NA  | Inf
 Left-Censored     |-Inf |  NA  | bi 
 Interval-Censored | ai  |  NA  | bi


**Right Censored data BGLR and survreg**

```R
 library(BGLR)
 n=500
 x1=rnorm(n);x2=rnorm(n)
 y=100+x1*.5-x2*.8+rnorm(n)
 yNotCensored=y
 
 
 # Introducing censoring (for simplicity done at fixed points -1 and 1)
 isCensored=sample(c(T,F),size=n,replace=T)
 
 y[isCensored]=NA
 a=rep(NA,n)
 a[isCensored]=yNotCensored[isCensored]+runif(sum(isCensored))
 b=rep(NA,n)
 b[isCensored]=Inf
 head(cbind(a,y,b))
 
## Survreg (Maximum Likelihood)
 library(survival)
 time=y ; time[isCensored]=a[isCensored]
 Y=Surv(time=time,event=!isCensored)
 fm0=survreg(Y~x1+x2,dist='gaussian')
 summary(fm0)
 
## BGLR
 Z=scale(cbind(x1,x2),center=T,scale=F) # centering helps mixing
 fm=BGLR(y=y,a=a,b=b,ETA=list(list(X=Z,model='FIXED')),nIter=5100,burnIn=100)
 
 # Compare the following with summary(fm0)
 cbind(c(fm$mu,fm$ETA[[1]]$b), c(fm$SD.mu,fm$ETA[[1]]$SD.b))
 sqrt(fm$varE) # compare with 'scale' in summary(fm0)

```
