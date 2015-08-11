### Fitting models with heterogeneous error variances

```R
library(devtools)
install_git('https://github.com/gdlc/BGLR-R/')
library(BGLR)

varGroups=sample(1:3,size=1000,replace=T)
VAR=c(.1,.5,1)
error=rnorm(sd=sqrt(VAR)[varGroups],n=nrow(X))

# Example 1
#       In the example below groups defines groups of res. var.
#       The model is just an intercept model but the same
#       can be done with any of the regression models implemented
#       in BGLR, except RKHS (see 2nd example no how to fit GBLUP with heterogeneous error variances).

  fm=BGLR(y=error, groups=varGroups)
  cbind(VAR,fm$varE)
 
#Example 2 (RKHS)
  data(wheat)
  G=tcrossprod(scale(wheat.X)); G=G/mean(diag(G))
  y=wheat.Y[,1]
  varGroups=sample(1:2,n=nrow(G),replace=T)
  y=y*varGroups # error variance is twice as large in group 2
  
  # RKHS is not implemented for heteroskedastic error variances 
   fm=BGLR(y=y,ETA=list(list(K=G,model='RKHS')),groups=varGroups)

  # However, you can always use PCs!
   EVD=eigen(G)
   PC=EVD$vectors[,EVD$values>1e-4]
   for(i in 1:ncol(PC)){ PC[,i]=PC[,i]*sqrt(EVD$values[i])}

   fm=BGLR(y=y,ETA=list(list(X=PC,model='BRR')),groups=varGroups)

```R
