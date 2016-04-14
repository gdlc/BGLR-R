## Integrating multiple Omics for Prediction of BC Survival

The following scripts illustrate how to fit some of the models presented in [Vazquez et al., Genetics (2016)]().

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
 
```R
 library(BGData)
 library(BGLR)
 load('OMIC_DATA.RData')
 
 Gge<-getG(Xge,scaleCol=T,scaleG=T) # Similarity matrix for gene expression.
 Gmt<-getG(Xmt,scaleCol=T,scaleG=T) # Similarity matrix for methylation. 
 XF<- scale(XF, scale=FALSE, center=TRUE) # centering and scaling the incidence matrix for fixed effects.
```

#### (3)  Fitting a survival model for Fixed effects using BGLR

```R

```
