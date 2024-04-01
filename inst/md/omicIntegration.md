## Integrating multiple Omics for Prediction of BC Survival

The following scripts illustrate how to fit some of the models presented in [Vazquez et al., Genetics (2016)](https://doi.org/10.1534/genetics.115.185181).

#### (1) Installing required libraries

The code below illustrates how to install BGLR and BGData from GitHub. BGLR can also be installed from CRAN using `install.packages()`.

```R
 install.packages(pkg='devtools',repos='https://cran.r-project.org/')    #1# install devtools
 library(devtools)                                                       #2# load the library
 install_git('https://github.com/gdlc/BGLR/')                            #3# install BGLR from GitHub
 install_git('https://github.com/quantgen/BGData/')                      #4# install BGLR from GitHub
```   

#### (2) Loading data and computing similarity matrices.
 
 The function `getG()`, from the [BGData](https://github.com/quantgen/BGData) R-package computes a similarity matrix of the form G=XX'. The function offers several alternatives relative to centering and scaling. The function also allows multi-core computing. For further details follow the link provided above.
 
 **Data**. The code assumes that the user has saved in the file `OMIC_DATA.RDasta` the objects that contain the phenotypic, covariates and omic information
  * XF: 
 
```R
 library(BGData)
 library(BGLR)
 load('OMIC_DATA.RData')
 
 Gge<-getG(Xge,scaleCol=T,scaleG=T) # Similarity matrix for gene expression.
 Gmt<-getG(Xmt,scaleCol=T,scaleG=T) # Similarity matrix for methylation. 
 XF<- scale(XF, scale=FALSE, center=TRUE) # centering and scaling the incidence matrix for fixed effects.
```

#### (3)  Fitting a survival model for Fixed effects using BGLR

The following code illustrates how to use BGLR to fit a fixed effects model. The matrix XF is an incidence matrix for effects. There is no column for intercept in XF because BGLR adds the intercept authomatically. The response variable `y` is assumed to be coded with two lables (e.g., 0/1), the argument `response_type` is used to indicate to BGLR that the response is ordinal (the binary case is a special case with only two levels). Predictors are given to BGLR in the form a two-level list. The argument `save_at` can be used to provide a path and a pre-fix to be added to the files saved by BGLR. For further details see [Pérez and de los Campos, Genetics, 2014](https://doi.org/10.1534/genetics.114.164442) The code also shows how to retrieve estimates of effects and of psuccess probabilities.

```R
 # Inputs
  ETA<-list( COV=list(X=XF, model='FIXED') )
  nIter=12000; burnIn=2000
 
 # Fitting the model
  fm=BGLR(y=y, ETA=LP,nIter=nIter,burnIn=burnIn, saveAt='cov_', response_type='ordinal')
 
 # Retrieving estimates
  fm$ETA$COV$b      # posterior means of fixed effects
  fm$ETA$COV$SD.b   # posteriro SD of fixed effects
  head(fm$probs)    # estimated probabilities for the 0/1 outcomes.
```
