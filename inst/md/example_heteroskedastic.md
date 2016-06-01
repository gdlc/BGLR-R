### Fitting models with heterogeneous error variances

The following example illsutrates how to fit models with heterogeneous variances in BGLR. The same option can be used with any of the models fitted by BGLR, except RKHS.

```R
  varGroups=sample(1:3,size=1000,replace=T)
  VAR=c(.1,.5,1)
  error=rnorm(sd=sqrt(VAR)[varGroups],n=length(varGroups))
  fm=BGLR(y=error, groups=varGroups)
  cbind(VAR,fm$varE)
```

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
