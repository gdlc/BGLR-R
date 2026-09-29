## RKHS regression with heterogeneous residual variances

The `groups` argument permits a separate residual variance for each group in a
single-response Gaussian model containing one or more `RKHS` terms. This
extension supports ordinary missing responses (`NA`). Grouped RKHS models with
ordinal or censored responses are outside its scope.

The Bayesian regression and kernel formulation follows
[Pérez and de los Campos (2014), Genetics 198:483–495](https://doi.org/10.1534/genetics.114.164442)
([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4196607/)).
[López-Cruz et al. (2015), G3 5:569–582](https://doi.org/10.1534/g3.114.016097)
([full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC4390573/)) provides a related
formulation combining shared and environment-specific genetic effects. Its
across-environment and interaction analyses assumed homogeneous residual
variances, as stated in that article's Software section. The heterogeneous
residual update below is derived for the current implementation.

### Model and observation weights

For observations in groups $g(i)\in\{1,\ldots,G\}$, consider

$$
y=\mathbf{1}\mu+X\beta+\sum_{j=1}^{J}u_j+\varepsilon,
\qquad u_j\mid\tau_j\sim N(0,\tau_jK_j),
\qquad \varepsilon_i\mid\sigma^2_{g(i)}
\sim N\left(0,\frac{\sigma^2_{g(i)}}{w_i^2}\right).
$$

Residuals and the distinct random-effect terms are conditionally independent.
Each $K_j$ is a symmetric positive semidefinite matrix indexed in the same order
as `y`. This derivation assumes finite, strictly positive `weights`; $w_i$
denotes the actual argument value. Thus the relative residual precision is
$w_i^2$.
This convention follows the current code's multiplication of the response and
predictors by `weights`; the weights in the 2014 article's introductory
variance expression use the reciprocal convention.

Set $W=\operatorname{diag}(w_i)$ and $R=\operatorname{diag}(\sigma^2_{g(i)})$.
On the working scale,

$$
y^*=Wy=W\mathbf{1}\mu+WX\beta+\sum_j u_j^*+\varepsilon^*,
\qquad u_j^*=Wu_j,\qquad \varepsilon^*\sim N(0,R).
$$

The `groups` argument changes the residual covariance only. Group-specific
means, shared genetic effects, and genetic interactions require their own
predictors or kernels in `ETA`.

### Kernel representation

For one RKHS term, write the retained eigendecomposition as

$$
WK_jW\simeq V_jD_jV_j^\mathsf{T},\qquad
V_j^\mathsf{T}V_j=I_{r_j},\qquad
D_j=\operatorname{diag}(d_{j1},\ldots,d_{jr_j}),
$$

where `tolD` retains eigenvalues $d_{jk}>\texttt{tolD}$. With positive retained
eigenvalues, the equivalent coefficient representation is

$$
u_j^*=V_ja_j,\qquad a_j\mid\tau_j\sim N(0,\tau_jD_j).
$$

Here $a_j$, $r_j$, and $\tau_j$ correspond to `uStar`, `levelsU`, and `varU`.
If positive eigenvalues are discarded, inference uses the resulting reduced
rank covariance. Its original-scale kernel is
$W^{-1}V_jD_jV_j^\mathsf{T}W^{-1}$, which equals $K_j$ when no nonzero
eigenvalues are discarded.

When `K` is supplied, BGLR forms and decomposes $WK_jW$. When `V` and `d` are
supplied, they must already describe that weighted kernel, with orthonormal
columns in `V`. They are not weighted again. Multiplying the eigenvectors of
$K_j$ by $W$ does not generally produce an eigendecomposition of $WK_jW$.

### Conditional coefficient updates

Holding the other terms and variance parameters fixed, define
$\rho_j=y^*-W\mathbf{1}\mu-WX\beta-\sum_{\ell\ne j}V_\ell a_\ell$.
The joint full conditional is

$$
a_j\mid\text{rest}\sim N(Q_j^{-1}h_j,Q_j^{-1}),\qquad
Q_j=V_j^\mathsf{T}R^{-1}V_j+\tau_j^{-1}D_j^{-1},\qquad
h_j=V_j^\mathsf{T}R^{-1}\rho_j.
$$

Under homogeneous residual variance,
$V_j^\mathsf{T}R^{-1}V_j=I/\sigma^2$, which permits the existing independent
normal updates. With heterogeneous residual variances, this matrix is
generally dense. The grouped implementation therefore makes a sequential
Gibbs update of the coordinates of $a_j$.

Suppress $j$ temporarily and let $e=y^*-W\mathbf{1}\mu-WX\beta-\sum_jV_ja_j$
be the current full residual. For coordinate $k$, precompute
$s_{gk}=\sum_{i:g(i)=g}v_{ik}^2$, stored as `x2`. Its conditional precision and
linear term are

$$
c_k=\sum_g\frac{s_{gk}}{\sigma_g^2}+\frac{1}{\tau d_k},
\qquad
t_k=\sum_g\frac{
\sum_{i:g(i)=g}v_{ik}e_i+a_k^{\mathrm{old}}s_{gk}
}{\sigma_g^2}.
$$

The update is

$$
a_k^{\mathrm{new}}=\frac{t_k}{c_k}+\frac{z_k}{\sqrt{c_k}},
\qquad z_k\sim N(0,1),\qquad
e\leftarrow e+v_k(a_k^{\mathrm{old}}-a_k^{\mathrm{new}}).
$$

This is the calculation performed by `sample_beta_groups` with coefficient
prior variances `varU * d`. The residual enters the routine with the old RKHS
term still subtracted; adding back the entire term before this call would
count the old coefficients twice. The routine's coordinate adjustments
preserve the full residual for subsequent terms and variance updates.
Reconstructing `u = V %*% uStar` afterward requires no additional residual
adjustment. Setting `minAbsBeta = 0` preserves the unmodified normal draw.

A complete coordinate sweep leaves the joint conditional distribution
invariant; it is not an independent draw from that multivariate distribution.
Strong conditional correlations can therefore reduce effective sample sizes.

### Variance components and missing responses

Use BGLR's parameterization
$p(v)\propto v^{-(\nu/2+1)}\exp\{-S/(2v)\}$, equivalently
$v=S/\chi^2_\nu$ for a proper prior. For an RKHS term,

$$
\tau_j\mid\text{rest}
\sim\frac{S_j+\sum_{k=1}^{r_j}a_{jk}^2/d_{jk}}
{\chi^2_{\nu_j+r_j}}.
$$

Thus `S0` and `df0` in that `ETA` term combine with the retained rank $r_j$,
not the number of observations. No residual-variance factor enters this
conditional because the kernel prior is $N(0,\tau_jK_j)$.

For group $g$, with $n_g$ records in the augmented data,

$$
\sigma_g^2\mid\text{rest}
\sim\frac{S_e+\sum_{i:g(i)=g}e_i^2}{\chi^2_{\nu_e+n_g}}.
$$

The current grouped implementation uses the same residual `S0` and `df0` for
each group. RKHS terms add no further sum of squares or degrees of freedom to
this residual-variance conditional.

For an ordinary missing response, the working-scale draw is

$$
y_i^*\mid\text{rest}\sim N(\eta_i^*,\sigma_{g(i)}^2),
\qquad e_i\leftarrow y_i^*-\eta_i^*.
$$

Consequently the standard deviations for the missing-response vector must be
selected using `sqrt(varE)[groups[whichNa]]`. Residual-variance updates include
all augmented records, including the current missing-response draws. This is
a valid data augmentation scheme: integration over a missing response
contributes one to the observed-data likelihood. Missing responses must be
excluded from observed-data likelihood summaries.

### Output scales and equivalent BRR models

For positive weights, original-scale kernel summaries are

$$
E(u_i\mid y)=\frac{E(u_i^*\mid y)}{w_i},\qquad
\operatorname{SD}(u_i\mid y)=\frac{\operatorname{SD}(u_i^*\mid y)}{w_i}.
$$

The returned RKHS `u` and `SD.u` use this original scale. The returned `V` and
`uStar` retain the working-scale representation. Under the existing BGLR
output convention, `yHat` and `SD.yHat` also remain on the working scale;
divide these two vectors by `weights` to obtain original-scale predictor
summaries. At posterior mean parameters, the Gaussian likelihood uses
original-scale residuals and standard deviations
`sqrt(varE[group]) / weights`, with the posterior mean of every group
variance.

Set $b_j=D_j^{-1/2}a_j$. Then $b_j\mid\tau_j\sim N(0,\tau_jI)$ and

$$
u_j=X_j^{\mathrm{BRR}}b_j,\qquad
X_j^{\mathrm{BRR}}=W^{-1}V_jD_j^{1/2}.
$$

Passing this matrix as the original-scale `X` of a `BRR` term gives the same
Gaussian model after BGLR applies its weights. For unit weights this reduces
to $X_j^{\mathrm{BRR}}=V_jD_j^{1/2}$. Model equivalence requires the same
retained kernel, residual model, and explicit variance priors: match RKHS
`varU` to BRR `varB`, and match both terms' `df0` and `S0`.
Their automatic prior scales differ: RKHS uses the mean retained eigenvalue,
whereas BRR uses the sum of centered working-scale predictor mean squares.
Equal `R2` values alone therefore do not ensure equal priors. The BRR
sampler's positive minimum coefficient magnitude can also preclude numerical
identity when a draw is extremely close to zero.

### Numerical validation

Useful checks distinguish the conditional calculations from posterior
simulation accuracy:

1. Compare one grouped coordinate sweep and its final residual with an
   independent calculation from the dense precision $Q_j$, using the same
   normal variates and heterogeneous group variances.
2. Compare RKHS and its BRR factor representation using explicitly matched
   priors, several kernels, non-unit weights, and missing responses. Check
   variance components and original-scale fitted values as well as
   predictions at missing records.
3. Compare `K` input with its supplied weighted `V, d` decomposition, including
   a rank-one kernel. With a single residual group, compare against the
   existing homogeneous RKHS calculation.
4. Check missing-response draws against their group-specific conditional
   variances and recompute the observed-data likelihood at posterior mean
   parameters. Simulation comparisons should allow for posterior uncertainty
   and autocorrelation-adjusted Monte Carlo error; correlation between two
   prediction vectors alone does not establish equivalence.

[Back to examples](https://github.com/gdlc/BGLR-R/blob/master/README.md)
