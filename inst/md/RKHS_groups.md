## Single-response MET GBLUP with heterogeneous residual variances

`BGLR()` accepts `groups` with one or more `model="RKHS"` terms for
uncensored Gaussian responses, including ordinary missing responses (`NA`).
Each observation belongs to one residual-variance group. For multi-environment
trials (MET), use the environment as the group. The response remains a vector.

For observation $i$ of genotype $k$ in environment $h$, the example fits

$$
y_i=\mu+\alpha_h+g_k+\delta_{kh}+\varepsilon_i,
\qquad \varepsilon_i\sim N(0,\sigma_{e,h}^2).
$$

The shared genomic effects have covariance $\sigma_G^2G$ and the interaction
effects have covariance $\sigma_{GE}^2G$ within an environment, independently
between environments. Residuals are independent between observations.
This genetic structure is related to the shared and environment-specific
effects in [López-Cruz et al. (2015)](https://doi.org/10.1534/g3.114.016097).
The model here also estimates a residual variance for each environment;
that article's across-environment analyses used a common residual variance.
See the [conditional distributions and implementation derivation](RKHS_groups_derivation.md)
for the heterogeneous residual extension of BGLR's Bayesian kernel model
([Pérez and de los Campos, 2014](https://doi.org/10.1534/genetics.114.164442)).

### Constructing the observation-level kernels

Let `d` contain `y`, `genotype`, and `environment`, with one row per record,
and let `G` be a genomic relationship matrix with genotype IDs as its row and
column names. Align IDs explicitly rather than relying on factor codes.

```r
library(BGLR)
stopifnot(identical(rownames(G),colnames(G)),!anyDuplicated(rownames(G)))
id=match(as.character(d$genotype),rownames(G))
stopifnot(!anyNA(id),!anyNA(d$environment))
environment=factor(d$environment)

K_G=G[id,id,drop=FALSE]
K_GE=K_G*outer(as.character(environment),as.character(environment),"==")
X_E=model.matrix(~environment)[,-1,drop=FALSE]
```

Equivalently, $K_G=Z_GGZ_G'$ and $K_{GE}=K_G\circ(Z_EZ_E')$, where the
incidence matrices map observations to genotypes and environments. Repeated
genotypes and unequal record counts are allowed; rows of both kernels must
follow the order of `d$y`. The environment contrasts have one fewer column
than the number of environments because BGLR includes an intercept.

### Fitting and extracting results

```r
ETA=list(environment=list(X=X_E,model="FIXED"),
         G=list(K=K_G,model="RKHS"),
         GE=list(K=K_GE,model="RKHS"))

set.seed(2403)
fm=BGLR(y=d$y,ETA=ETA,groups=environment,
        nIter=12000,burnIn=2000,thin=5,saveAt="MET_")

fm$varE                 # posterior residual variance means, named by group
fm$SD.varE              # posterior standard deviations
fm$ETA$G$varU           # shared genomic variance
fm$ETA$GE$varU          # G x E variance
fm$ETA$G$u+fm$ETA$GE$u  # genetic values for each observation
fm$yHat                # fitted values and predictions at NA responses
```

`groups` specifies residual variances only; it does not add environment means
or G×E effects automatically. The example uses a common interaction variance.
Separate environment-specific interaction kernels can be supplied as additional
RKHS terms when different genetic interaction variances are required. An
unrestricted genetic covariance matrix across environments is a different model.

For prediction, set held-out responses to `NA` and retain their kernel rows
and group labels. The variance of each missing-response draw is selected from
its own group. Every environment should have observed responses to estimate
its mean and residual variance. Validation by held-out cells, new genotypes,
and new environments addresses different prediction targets.

Residual variances use the shared scalar prior parameters `S0` and `df0` supplied
to `BGLR()`. Each kernel has its own `S0` and `df0` in `ETA`. Existing default
prior calculations are unchanged. For comparison with a BRR factorization,
specify matching variance priors explicitly because automatic RKHS and BRR
prior scales differ.

With positive observation weights, residual variance is
`fm$varE[group] / weights^2`. The existing output convention is retained:
RKHS `u` and `SD.u` are on the response scale, while `yHat` and `SD.yHat`
must be divided by `weights` to obtain that scale. If supplying `V` and `d`
instead of `K`, compute the eigendecomposition of `W K W`, with
`W=diag(weights)`, before passing them to BGLR.

### Reproducible wheat example

```r
demo("RKHS_groups",package="BGLR",ask=FALSE)
```

The demonstration uses 120 wheat genotypes in four mega-environments,
withholding 10% of the genotype-environment cells. The bundled wheat responses
are standardized genotype means, so residual variances are on that scale.
Its explicit variance priors illustrate reproducibility; they should be
reviewed for a particular dataset. Assess trace stability, effective sample
sizes, and agreement across independent chains before interpreting results.

The grouped update uses the package's registered C coefficient sampler. Under
heterogeneous residual variances, kernel eigen-effects are generally correlated
conditionally, and sequential Gibbs updates can mix more slowly than the
homogeneous update. No additional package dependency is required. Grouped
RKHS with ordinal or censored responses is not supported by this extension.

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
