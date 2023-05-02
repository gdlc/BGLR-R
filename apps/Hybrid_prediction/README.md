## Hybrid prediction using GBLUP

This app fits GBLUP models to predict hybrid performance. The problem is to predict hybrid performance 
using genotipic information from parents. The data that we have available is as follows:

- Genotypic information from males and females (genomic relationship matrixes).
- Phenotypic information for some hybrids.

The desired output is the prediction of hybrids.

For more details see [Acosta-Pech et al., 2017](https://pubmed.ncbi.nlm.nih.gov/28401254/) and references therein.

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
