
## BGLR: An R Package for (Bayesian) High-Dimensional Regression

[![CRAN status](https://www.r-pkg.org/badges/version/BGLR?color=green)](https://CRAN.R-project.org/package=BGLR)
[![CRAN checks](https://badges.cranchecks.info/worst/BGLR.svg)](https://cran.r-project.org/web/checks/check_results_BGLR.html)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/BGLR)](http://www.r-pkg.org/pkg/BGLR)
[![Downloads](http://cranlogs.r-pkg.org/badges/BGLR?color=blue)](http://www.r-pkg.org/pkg/BGLR)


The BGLR Package ([Perez & de los Campos, 2014](https://doi.org/10.1534/genetics.114.164442)) implements a variety of shrinkage and variable selection regression procedures. In this repository we maintain the latest
version beta version. The latest stable release can be downloaded from [CRAN](https://cran.r-project.org/web/packages/BGLR/index.html).


#### Citation

Please cite [Perez & de los Campos, 2014](https://doi.org/10.1534/genetics.114.164442) and [Perez & de los Campos, 2022](https://doi.org/10.1093/genetics/iyac112) for BGLR and Multitrait, respectively.


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

This can be fixed it by installing gfortan, for mac os you can use [this](https://github.com/fxcoudert/gfortran-for-macOS/releases)

Useful references:

 - [link1](https://yiqingxu.org/public/BigSurError.pdf).
 - [link2](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).


## 1. Single-Trait Models
----------------------------------------------------------------

**Examples BGLR-function**
----------------------------------------------------------------

  - [1. Parametric Bayesian Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BayesianAlphabet.md)
  
  - [2. Estimating or Fixing Hyper-parameters?](https://github.com/gdlc/BGLR-R/blob/master/inst/md/hyperparameters.md)
  
  - [3. GBLUP: various implementations with BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GBLUP.md)
  
  - [4. Prediction in testing sets: three methods](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Validation.md)
  
  - [5. Semi-parametric regression (RKHS) using BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS.md)
  
  - [6. Fitting Models with Multiple sets of Effects ("Mixed-Effects Model")](https://github.com/gdlc/BGLR-R/blob/master/inst/md/mixedModel.md)
  
  - [7. Saving Samples of Effects in Binary Files](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_saveEffects.md)
  
  - [8. Heritability Estimation: two methods](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heritability.md)

  - [9. Heterogeneous Error Variance Models](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_heteroskedastic.md)
  
  - [10. GxE Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_usingInteractions.md)
  
  - [11. Modeling Genetic Heterogeneity Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heterogeneity_interactions.md)
  
  - [12. Estimating the Proportion of Variance Explained by Principal Components](https://github.com/gdlc/BGLR-R/blob/master/inst/md/PC.md)
  
  - [13. Categorical (binary and ordinal) Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/categorical.md)
  
  - [14. Censored Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/censored.md)


  - [15. Bayesian regressions with markers sets: an example of how BGLR can be used to fit models with set-specific priors](https://github.com/gdlc/BGLR-R/blob/master/inst/md/setMethods.md)
  
  - [16. BRR-sets: Guassian prior with set-specific variances](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BRR_sets.md)
  
  - [17. Multi-trait prediction using eigenvectors](https://github.com/gdlc/BGLR-R/blob/master/inst/md/MULTITRAIT.md)

  - [18. Two-steps Finlay-Wilkinson Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/FW_BGLR.md)

**Other Omics**
----------------------------------------------------------------

  - Metabolomics (milk-spectra): [Ferragina et al., J.D.Sci, 2015](https://doi.org/10.3168/jds.2014-9143)
  
  - Multi-omic (gene expression, methylation & CNV): [Vazquez et al. (10th, WCGALP)](https://www.researchgate.net/profile/Sadeep_Shrestha/publication/266077932_Integration_of_Multi-Layer_Omic_Data_for_Prediction_of_Disease_Risk_in_Humans/links/542430f30cf26120b7a72201.pdf) & [Vazquez et al. (Genetics, 2016)](https://doi.org/10.1534%2Fgenetics.115.185181)

**Markers or Pedigree and Environmental Covariates**
----------------------------------------------------------------
  -Wheat (SNPs and env. covariates): [Jarquin et al. (TAG, 2014)](https://doi.org/10.1007%2Fs00122-013-2243-1)
  
  -Cotton (Pedigree and env. covariates): [Perez-Rodriguez et al.(Crop. Sci, 2015)](https://doi.org/10.2135/cropsci2014.08.0577)
  
**Image Data**
----------------------------------------------------------------
  -Maize (Image data): [Aguate et al. (Crop. Sci, 2017)](https://doi.org/10.2135/cropsci2017.01.0007)

## 2. Multi-trait models
----------------------------------------------------------------

The Multitrait function included in the BGLR package fits Bayesian multitrait models with 
arbitrary number of random effects using a Gibbs sampler. A functionality similar to this 
is implemented in the [MTM](http://quantgen.github.io/MTM/vignette.html) package. In
this implementation is possible to include regression on markers directly assigning Spike-slab or 
Gaussian priors for the regression coefficients and fixed effects can be different for all the traits.
We also have improved the sampling routines to speed up computations. Next we include some examples.

- [1. Random effect model with unstructured covariance matrices](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_UN_UN.md)
- [2. Diagonal residual covariance matrix](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_UN_DIAG.md)
- [3. Factor analysis](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_FA_DIAG.md)
- [4. Recursive models](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_REC_DIAG.md)
- [5. Mixed effects model](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS_FIXED_UN_UN.md)
- [6. Ridge Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BRR_UN_UN.md)
- [7. Ridge Regression + Additive relationship matrix](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BRR_UN_RKHS_UN.md)
- [8. Spike Slab](https://github.com/gdlc/BGLR-R/blob/master/inst/md/SS_UN_UN.md)
- [9. Fixed effects](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Fixed_effects_multitrait.md)

**Supplementary scripts for the draft of the paper "Multi-trait Bayesian Shrinkage and Variable Selection Models with the BGLR R-package"**

- [1. Manuscript scripts](https://github.com/gdlc/BGLR-R/blob/master/inst/md/MAUSCRIPT_SCRIPTS.Rmd)
- [2. SupplementaryFileS2_Scripts](https://github.com/gdlc/BGLR-R/blob/master/inst/md/SupplementaryFileS2_SCRIPTS.Rmd)

