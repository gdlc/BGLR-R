## Two-steps Finlay-Wilkinson Regression using BGLR

**Authors**: Gustavo de los Campos <gustavoc@msu.edu> ,  Paulino Perez-Rodriguez <perpdgo@gmail.com>  &  BOGARD Matthieu <M.BOGARD@arvalis.fr> 

The following data and code is reproduced from the supplementary data from [de los Campos et al., Nat. Comm., 2020](https://www.nature.com/articles/s41467-020-18480-y).

**Citation**: de los Campos G., Pérez-Rodríguez P., Bogard M., Gouache D., and Crossa J. *A data-driven simulation platform to predict cultivars’ performances under uncertain weather conditions* Nature Communications  11, Article number: 4876 (2020).

**Before you run the code below, dowload the following**: 

  - [FW.BGLR](https://github.com/gdlc/BGLR-R/blob/master/inst/md/FW_BGLR.R)- function
  - [Sample Data](https://github.com/gdlc/BGLR-R/tree/master/inst/md/sample_FW.RData)

Here we present a function to perform the Finlay-Wilkinson (1963)
analysis in two steps using the BGLR (Pérez and de los Campos) package
in R (R Core Team, 2019). We created the R function 
FW.BGLR to perform the analysis. The function takes the 
following arguments:

	- pheno: a data.frame with 3 columns, VAR (variety), 
	  ENV (environment), y (response variable).
	- X: matrix of markers coded for additive effects (e.g, 0, 1, 2).
	- G: matrix with genomic relationships between individuals.

The matrix **X** contains marker information for the varieties given in the 
data.frame. This matrix is used to compute the additive relationships between 
individuals (Lopez-Cruz et al., 2015) if matrix **G** is not given. 

The function returns a list object with the following elements:m

	- yHat: vector with predicted values for the response variable.
	- VAR: a data frame with the columns ID (varieties), 
	  int (estimated intercept), intSD (estimated standard deviation for intercept),
	  slope (estimated slope) and  slopeSD (estimated standard deviation for slope).
	- ENV: vector with environmental effects.
	
The following R code shows how to load sample data (raw means).

```r
 # Load data
  load("sample_FW.RData")

 #list objects, at this point you should have at least 2 objects: pheno and G.
  ls()

```
	
The code below shows how to load the function and perform the analysis assuming 
that objects pheno and G are already loaded in the R environment.


```r
 # Loads BGLR library
  library(BGLR)

 # Loads function for FW analysis
  source("FW_BGLR.R")

 # Fits FW-regression in two steps
  fm<-FW.BGLR(pheno=pheno,G=G,verbose=FALSE)

 # Predictions
  head(fm$yHat)

 # Intercept and slopes
  head(fm$VAR)

 # Environmental effects
  head(fm$ENV)

 # Plot Environment effect vs Genotype performance
  plot.FW(fm,pheno)

```

### References

Finlay K.W. and G.N. Wilkinson. 1963. “The Analysis of Adaptation in a Plant-Breeding Programme.” Australian Journal of Agricultural Research 14 (6). CSIRO PUBLISHING:742. https://doi.org/10.1071/AR9630742.

Lopez-Cruz M., J. Crossa, D. Bonnett, S. Dreisigacker, J. Poland, J.L. Jannink, R.P. Singh, E. Autrique, and G. de los Campos. 2015. “Increased Prediction Accuracy in Wheat Breeding Trials Using a Marker x Environment Interaction Genomic Selection Model.” G3 (Bethesda, Md.) 5 (4). G3: Genes, Genomes, Genetics:569–82. https://doi.org/10.1534/g3.114.016097.

Pérez, P. and G. de los Campos. 2014. Genome-Wide Regression and Prediction with the BGLR Statistical Package. Genetics, 198: 483-495.

R Core Team. 2019. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. URL https://www.R-project.org/.
