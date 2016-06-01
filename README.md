## BGLR: An R Package for (Bayesian) High-Dimensional Regression

The BGLR Package ([Perez-Rodriguez & de los Campos, 2014](http://www.ncbi.nlm.nih.gov/pubmed/25009151])) implements a variety of shrinkage and variable selection regression procedures. In this repository we maintain the latest
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
   install_git('https://github.com/gdlc/BGLR/')                          #3# install BGLR from GitHub
```

*Note*: when trying to install from github on a mac you may get the following error message

```
ld: library not found for -lgfortran
```

This can be fixed it by following the following [advise](http://thecoatlessprofessor.com/programming/rcpp-rcpparmadillo-and-os-x-mavericks-lgfortran-and-lquadmath-error/).

 

###Examples

### BGLR-Examples

  - [1. Installing BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Install.md)
  
  - [2. Parametric Bayesian Regression](https://github.com/gdlc/BGLR-R/blob/master/inst/md/BayesianAlphabet.md)
  
  - [3. Estimating or Fixing Hyper-parameters?](https://github.com/gdlc/BGLR-R/blob/master/inst/md/hyperparameters.md)
  
  - [4. GBLUP: various implementations with BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GBLUP.md)
  
  - [5. Prediction in testing sets: three methods](https://github.com/gdlc/BGLR-R/blob/master/inst/md/Validation.md)
  
  - [6. Semi-parametric regression (RKHS) using BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/RKHS.md)
  
  - [5. Fitting Models with Multiple sets of Effects ("Mixed-Effects Model")](https://github.com/gdlc/BGLR-R/blob/master/inst/md/mixedModel.md)
  
  - [9. Saving Samples of Effects in Binary Files](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_saveEffects.md)
  
  - [10. Heritability Estimation: two methods]()
  
  - [11. GxE Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/GxE_interactions.md)
  
  - [12. Heterogeneous Error Variance Models](https://github.com/gdlc/BGLR-R/blob/master/inst/md/example_heteroskedastic.md)
  
  - [13. Modeling Genetic Heterogeneity Using Interactions](https://github.com/gdlc/BGLR-R/blob/master/inst/md/heterogeneity_interactions.md)
  
  - [14. Estimating the Proportion of Variance Explained by Principal Components](https://github.com/gdlc/BGLR-R/blob/master/inst/md/PC.md)
  
  - [15. Categorical (binary and ordinal) Regression]()
  
  - [16. Censored Regression]()
  
 
