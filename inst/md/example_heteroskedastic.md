### Fitting models with heterogeneous error variances

The following example illustrates how to fit models with heterogeneous residual
variances in BGLR. The `groups` option supports `FIXED`, `BRR`, `BRR_sparse`,
`BayesB`, `BayesC`, and `RKHS` terms. For RKHS, grouped residual variances are
implemented for uncensored Gaussian responses, including ordinary missing
responses (`NA`).

```R
  varGroups=sample(1:3,size=1000,replace=T)
  VAR=c(.1,.5,1)
  error=rnorm(sd=sqrt(VAR)[varGroups],n=length(varGroups))
  fm=BGLR(y=error, groups=varGroups)
  cbind(VAR,fm$varE)
```

For a single-response GBLUP model with genomic main effects and genotype by
environment interactions, see the [MET RKHS example](RKHS_groups.md) and its
[conditional-distribution derivation](RKHS_groups_derivation.md).

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
