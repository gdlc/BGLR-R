## Fitting GBLUP model

This app fits GBLUP model. Consider the model:

$$
\boldsymbol y= \boldsymbol 1 \mu + \boldsymbol X \boldsymbol \beta + \boldsymbol Z \boldsymbol u + \boldsymbol e, 
$$

where

- $\boldsymbol y$ is the vector of phenotypes of dimmension $n \times 1$.
- $\boldsymbol 1$ is a vector of ones of dimmensions $n \times 1$.
- $\mu$ is an intercept.
- $\boldsymbol X$ is a matrix of fixed effects of dimmensions $n \times p$.
- $\boldsymbol \beta$ is a vector of fixed effects of dimmensions $p \times 1$.
- $\boldsymbol Z$ is an incidence matrix of dimensions $n \times r$.
- $\boldsymbol u$ is a vector of random effects of dimensions $r \times 1$, $\boldsymbol u \sim MN(\boldsymbol 0, \boldsymbol \sigma^2_u \boldsymbol K)$, 
  with $\boldsymbol K$ a relationship matrix of dimensions $r \times r$, $\sigma^2_u$ the variance parameter associated to $\boldsymbol u$.
-$\boldsymbol e$, the vector of residuals of dimmensions $n \times 1$, 
 $\boldsymbol e \sim MN(\boldsymbol 0, \boldsymbol \sigma^2_e \boldsymbol I)$, $\sigma^2_e$, the variance parameter associated to $\boldsymbol e$.
- We assume that $\boldsymbol e$ and $\boldsymbol u$ are independent.


The model can be rewritten as:

$$
\boldsymbol y= \boldsymbol 1 \mu + \boldsymbol X \boldsymbol \beta + \boldsymbol Z^{\ast} \boldsymbol u^{\ast} + \boldsymbol e, 
$$

where $\boldsymbol Z^\ast = \boldsymbol Z \boldsymbol \Gamma \boldsymbol \Lambda_1^{\frac{1}{2}}$, with $\boldsymbol \Gamma$, $\boldsymbol \Lambda$ 
the eigen-vectors and eigen values obtained from the eigen-value decomposition 
from $\boldsymbol K$, $\boldsymbol u^\ast \sim MN(\boldsymbol 0, \sigma^2_u \boldsymbol I)$.

The last model can be fitted in BGLR easily using ``Bayesian Ridge Regression'' and after the model is 
fitted the posterior mean of random effects is obtained as follows:

$$
\hat{\boldsymbol u} = \boldsymbol \Gamma \boldsymbol \Lambda^{\frac{1}{2}} \hat{\boldsymbol u}^\ast.
$$

Computations can be handled automatically by the app <em>RKHS.GBLUP</em> and user just needs to provide 
a matrix for fixed effects (or a formula to build the corresponding matrix), GIDs (Genotype Identifiers), 
matrix $\boldsymbol K$ and response variable.

## Example

For illustrative purposes we consider the mice dataset included in the package and, we predict the body mass index and 
we include as fixed effects the GENDER and the Litter. The relatioship between the individuals is taken into account including a random 
effect whose variance covariance matrix is proportional to the genomic relationship matrix that is built based on markers.

**Data preparation**

```r

#Load the library
library(BGLR)

#Load the data
data(mice)

#Phenotypes
pheno<-mice.pheno
pheno$GENDER<-as.factor(pheno$GENDER)
pheno$Litter<-as.factor(pheno$Litter)

#Markers
X<-mice.X

#Genomic relationship matrix
Z<-scale(X,center=TRUE,scale=TRUE)
G<-tcrossprod(Z)/ncol(Z)

```

**Sourcing the app**

```r
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/apps/GBLUP/RKHS_GBLUP.R')
```

**Fitting the model**

```r
 set.seed(456)
 fm<-RKHS.GBLUP(y=pheno$Obesity.BMI,
               FEffects=~GENDER+Litter,
               GID=pheno$SUBJECT.NAME,
               K=G,
               data=pheno,
               nIter=10000,
               burnIn=5000,
               thin=10,
               verbose=TRUE)
 ```

**Extracting results**
```r
#Posterior means for fixed effects
fm$ETA$FEffects$b


#Posterior means of random effects
fm$ETA$REffects$u

#Posterior mean of variance parameter associated to residuals
fm$varE

#Posterior mean of variance parameter associated to u
fm$ETA$REffects$varU

#Predictions of response variable
fm$yHat
```
