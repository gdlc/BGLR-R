
## BGLR: An R Package for (Bayesian) High-Dimensional Regression

The BGLR Package ([Perez & de los Campos, 2014](http://www.ncbi.nlm.nih.gov/pubmed/25009151)) implements a variety of shrinkage and variable selection regression procedures. In this repository we maintain the latest
version beta version. The latest stable release can be downloaded from [CRAN](https://cran.r-project.org/web/packages/BGLR/index.html).

#### Installation

**From CRAN (stable release)**.

```R
  install.packages(pkg='BGLR',repos='https://cran.r-project.org/')
```

**From GitHub (development version, added features)**.


```R
   install.packages(pkg='devtools',repos='https://cran.r-project.org/')  #1# install devtools
   library(devtools)                                                     #2# load the library
   install_git('https://github.com/gdlc/BGLR-R')                         #3# install BGLR from GitHub
```

*Note*: when trying to install from github on a mac you may get the following error message

```
ld: library not found for -lgfortran
```

This can be fixed it by following the following [advise](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).


## 1. Univariate responses
----------------------------------------------------------------

**Genomic Regressions**
----------------------------------------------------------------

  - [1. Parametric Bayesian Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BayesianAlphabet.md)
  
  - [2. Estimating or Fixing Hyper-parameters?](https://github.com/gdlc/BGLR-R/blob/master/inst/md/hyperparameters.md)
  
  - [3. GBLUP: various implementations with BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GBLUP.md)
  
  - [4. Prediction in testing sets: three methods](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Validation.md)
  
  - [5. Semi-parametric regression (RKHS) using BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS.md)
  
  - [6. Fitting Models with Multiple sets of Effects ("Mixed-Effects Model")](https://github.com/gdlc/BGLR-R/blob/master/inst/md/mixedModel.md)
  
  - [7. Saving Samples of Effects in Binary Files](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_saveEffects.md)
  
  - [8. Heritability Estimation: two methods](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heritability.md)

  - [10. Heterogeneous Error Variance Models](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_heteroskedastic.md)
  
  - [10. GxE Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md)
  
  - [11. Modeling Genetic Heterogeneity Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heterogeneity_interactions.md)
  
  - [12. Estimating the Proportion of Variance Explained by Principal Components](https://github.com/gdlc/BGLR-R/blob/master/inst/md/PC.md)
  
  - [13. Categorical (binary and ordinal) Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/categorical.md)
  
  - [14. Censored Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/censored.md)


  - [15. Bayesian regressions with markers sets: an example of how BGLR can be used to fit models with set-specific priors](https://github.com/gdlc/BGLR-R/blob/master/inst/md/setMethods.md)
  
  - [16. BRR-sets: Guassian prior with set-specific variances](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BRR_sets.md)
  
  - [17. Multi-trait prediction using eigenvectors](https://github.com/gdlc/BGLR-R/edit/master/inst/md/MULTITRAIT.md)


**Other Omics**
----------------------------------------------------------------

  - Metabolomics (milk-spectra): [Ferragina et al., J.D.Sci, 2015](http://www.sciencedirect.com/science/article/pii/S0022030215006645)
  
  - Multi-omic (gene expression, methylation & CNV): [Vazquez et al. (10th, WCGALP)](https://www.researchgate.net/profile/Sadeep_Shrestha/publication/266077932_Integration_of_Multi-Layer_Omic_Data_for_Prediction_of_Disease_Risk_in_Humans/links/542430f30cf26120b7a72201.pdf) & [Vazquez et al. (Genetics, 2016)](http://www.genetics.org/content/early/2016/04/27/genetics.115.185181.abstract)

**Markers or Pedigree and Environmental Covariates**
----------------------------------------------------------------
  -Wheat (SNPs and env. covariates): [Jarquin et al. (TAG, 2014)](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC3931944/)
  
  -Cotton (Pedigree and env. covariates): [Perez-Rodriguez et al.(Crop. Sci, 2015)](https://dl.sciencesocieties.org/publications/cs/abstracts/55/3/1143)
  
**Image Data**
----------------------------------------------------------------
  -Maize (Image data): [Aguate et al. (IBC, 2016)](https://dl.sciencesocieties.org/publications/cs/articles/57/5/2517)

## 2. Multi-trait models
----------------------------------------------------------------

The Multitrait function included in the BGLR package fits Bayesian multitrait models with 
arbitrary number of random effects using a Gibbs sampler a functionality similar to this 
is implemented in the [MTM](http://quantgen.github.io/MTM/vignette.html) package. In
this implementation is possible to include regression on markers directly assigning Spike-slab or 
Gaussian priors to the regression coefficients and fixed effects can be different for all the traits.
We also have improved the sampling routines to speed up computations. Next we include some examples.

- [1. Random effect model with unstructured covariance matrices](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_UN_UN.md)
- [2. Diagonal residual covariance matrix](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_UN_DIAG.md)
- [3. Factor analysis](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_FA_DIAG.md)
- [4. Recursive models]
- [5. Mixed effects model]
- [6. Ridge Regression]
- [7. Spike Slab]
- [8. Ridge Regression + Additive relationship matrix]

