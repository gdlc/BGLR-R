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
$$, 

where

- $\boldsymbol Z^\ast = \boldsymbol Z \boldsymbol \Gamma \boldsymbol \Lambda_1^{\frac{1}{2}}$, with $\boldsymbol \Gamma$, $\boldsymbol \Lambda$ the eigen-vectors and eigen values obtained from the eigen-value decomposition from $\boldsymbol K$, $\boldsymbol u^\ast \sim MN(\boldsymbol 0, \sigma^2_u \boldsymbol I)$.

The last model can be fitted in BGLR easily using ``Bayesian Ridge Regression'' and after the model is fitted the posterior mean of random effects is obtained as follows:

- $\hat{\boldsymbol g}_1 = \boldsymbol \Gamma_1 \boldsymbol \Lambda_1^{\frac{1}{2}} \hat{\boldsymbol g}_1^\ast$.
- $\hat{\boldsymbol g}_2 = \boldsymbol \Gamma_2 \boldsymbol \Lambda_2^{\frac{1}{2}} \hat{\boldsymbol g}_2^\ast$.
- $\hat{\boldsymbol h}_1 = \boldsymbol \Gamma_h \boldsymbol \Lambda_h^{\frac{1}{2}} \hat{\boldsymbol h}^\ast$.

Computations can be handled automatically by the app <em>BRR.Hybrid_prediction</em> and user just needs to provide 
locations, ids for parents, ids for hybrids, genomic relationship matrices and response variable.

## Example

For illustrative purposes we consider the maize dataset described in [Covarrubias-Pazaran (2016)](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0156744), and included in the <em>sommer package</em>. The dataset contains phenotypic data for plant 
height and grain yield for 100 out of 400 possible hybrids originated from 40 inbreed lines belonging to two heterotic groups, 20 lines each,
1600 rows exist for the 400 possible hybrids evaluated in 4 locations but only 100 crosses have phenotypic information. The purpose is to
predict the other 300 crosses.

**Data preparation**

```r

#Load the library
library(BGLR)

#Load the data, you need to download the file "cornHybrid.RData" included in this app
#This loads a list, from where we extract the information
load('cornHybrid.RData')

#Extract the hybid information, a data.frame with columns: 
#1)Location, 2)GCA1, ids for parent1, 
#3)GCA2, ids fro parent 2,
#4)SCA id for hybrids, 5)Yield and 6)PlantHeight
pheno<-cornHybrid$hybrid
head(pheno)

pheno$GCA1<-as.character(pheno$GCA1)
pheno$GCA2<-as.character(pheno$GCA2)
pheno$SCA<-as.character(pheno$SCA)

#Extract relationship matrix for both parents
G<-cornHybrid$K

#Genomic relationship matrix for parent 1
GCA1<-unique(pheno$GCA1)
selected<-rownames(G)%in%GCA1
G1<-G[selected,selected]
dim(G1)
rownames(G1)

#Genomic relationship matrix for parent 2
GCA2<-unique(pheno$GCA2)
selected<-rownames(G)%in%GCA2
G2<-G[selected,selected]
dim(G2)
rownames(G2)

#Generate H
#kronecker, make.dimmanes is necessary to identify the hybrids
#with the label Parent 1:Parent 2, using the same convention in pheno data.frame
H<-kronecker(G1,G2,make.dimnames=TRUE)

#At this point we need to have 4 objects:
#1)pheno
#2)G1
#3)G2
#4)H

```

**Sourcing the app**

```r
 source('https://raw.githubusercontent.com/gdlc/BGLR-R/master/apps/Hybrid_prediction/Hybrid_prediction.R')
```

**Fitting the model**

```r
 set.seed(456)
 fm<-BRR.Hybrid_prediction(y=pheno$Yield,
                            location=pheno$Location,
                            id1=pheno$GCA1,
                            id2=pheno$GCA2,
                            idH=pheno$SCA,
                            G1=G1,
                            G2=G2,
                            H=H,
                            nIter=10000,
                            burnIn=5000,
                            thin=10,
                            verbose=TRUE)
```

**Extracting results**
````r
#Variance component for parent 1
fm$ETA[[2]]$varB

#Variance component for parent 2
fm$ETA[[3]]$varB

#Variance component for hybrids
fm$ETA[[4]]$varB

#Variance component for error
fm$varE

#predictions
predictions<-data.frame(Loc=pheno$Location, yObs=pheno$Yield,yPred=fm$yHat,hybrid=pheno$SCA)
head(predictions)

#Posterior means for random effects
#Parent1
fm$ETA[[2]]$u

#Parent2
fm$ETA[[3]]$u

#Hybrid
fm$ETA[[4]]$u

````
