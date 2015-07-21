## BGLR: An R Package for (Bayesian) High-Dimensional Regression

The BGLR Package implements a variety of shrinkage and variable selection regression procedures. In this repository we maintain the latest
version beta version. The latest stable release can be downloaded from [CRAN](https://cran.r-project.org/web/packages/BGLR/index.html).

Technical details about the software and the methods implemented, as well as several examples can be found in the following article:

[Genome-wide regression and prediction with the BGLR statistical package (Genetics, 2014)](http://www.ncbi.nlm.nih.gov/pubmed/25009151)


Recent changes:
   - In V.1.0.4 We have implemented models for heterogeneous error variances (example).
   - More recently, in this repository, we:
            - added the possiblity of saving samples of effects in binary files (example),
            - implemented a new method (BRR-groups), this model is similar to BRR but uses one variance parameter per set of predictors (example),
            - added a function that computes, from samples of effects, the contribution of sets of SNPs to variance (example).
            
            
            
