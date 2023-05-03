## Hybrid prediction using GBLUP

This app fits GBLUP models to predict hybrid performance. The problem is to predict hybrid performance 
using genotipic information from parents. The data that we have available is as follows:

- Genotypic information from males and females (genomic relationship matrixes).
- Phenotypic information for some hybrids.

The desired output is the prediction of hybrids.

For more details see [Acosta-Pech et al.(2017)](https://pubmed.ncbi.nlm.nih.gov/28401254/) and references therein.

Briefly the statistical model for predicting hybrid performance is as follows:

$$
\boldsymbol y = \boldsymbol Z_E \boldsymbol \beta_E + \boldsymbol Z_1 \boldsymbol g_1 + \boldsymbol Z_2 \boldsymbol g_2 + \boldsymbol Z_h \boldsymbol h + \boldsymbol e, 
$$

where 

- $\boldsymbol Z_E$ is an incidence matrix for environments, $\boldsymbol \beta_E$ is the efect of the environments.
- $\boldsymbol Z_1$ is an incidence matrix for males.
- $\boldsymbol g_1$ is a vector or random effects for the general combining ability for males, $g_1 \sim MN(\boldsymbol 0, \sigma^2_1 \boldsymbol G_1)$, with $\boldsymbol G_1$ the genomic relationship matrix for males.
- $\boldsymbol Z_2$ is an incidence matrix for females.
- $\boldsymbol g_2$ is a vector or random effects for the general combining ability for females, $g_2 \sim MN(\boldsymbol 0, \sigma^2_2 \boldsymbol G_2)$, with $\boldsymbol G_2$ the genomic relationship matrix for females.
- $\boldsymbol Z_h$ is the incidence matrix for hybrids.
- $\boldsymbol h_f$ is a vector of random effects for specific combining ability of hybrids, $\boldsymbol h \sim MN(\boldsymbol \sigma^2_h \boldsymbol H)$, where $\boldsymbol H=\boldsymbol G_1 \otimes \boldsymbol G_2$, the symbol $\otimes$ denotes the Kronecker product of two matrices.
- $\boldsymbol e \sim MN(\boldsymbol 0, \sigma^2_e \boldsymbol I)$, with $\boldsymbol I$, the identity matrix.

The model can be rewritten as:

$$
\boldsymbol y = \boldsymbol Z_E \boldsymbol \beta_E + \boldsymbol Z_1^\ast \boldsymbol g_1^\ast + \boldsymbol Z_2^\ast \boldsymbol g_2^\ast + \boldsymbol Z_h^\ast \boldsymbol h^\ast + \boldsymbol e, 
$$

where

- $\boldsymbol Z_1^\ast = \boldsymbol Z_1 \boldsymbol \Gamma_1 \boldsymbol \Lambda_1^{\frac{1}{2}}$, with $\boldsymbol \Gamma_1$, $\boldsymbol \Lambda_1$ the eigen-vectors and eigen values obtained from the eigen-value decomposition from $\boldsymbol G_1$, $g_1^\ast \sim MN(\boldsymbol 0, \sigma^2_1 \boldsymbol I)$.
- $\boldsymbol Z_2^\ast = \boldsymbol Z_2 \boldsymbol \Gamma_2 \boldsymbol \Lambda_2^{\frac{1}{2}}$, with $\boldsymbol \Gamma_2$, $\boldsymbol \Lambda_2$ the eigen-vectors and eigen values obtained from the eigen-value decomposition from $\boldsymbol G_2$, $g_2^\ast \sim MN(\boldsymbol 0, \sigma^2_2 \boldsymbol I)$.

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
