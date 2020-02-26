\name{simulated3t}
\alias{simulated3t}
\title{simulated data for 3 traits}
\description{

Simulated dataset for three traits. Markers, QTL and phenotypes are
simulated for three traits. Here we extend the simulation scheme 
described by Cheng et al. (2018) for the case of three traits.  So we
simulated 100 evenly spaced loci on each of 4 chromosomes of length
10 cM. We selected 10 loci on each chromosome as QTL. We sampled
states from Bernoulli distribution with p=0.5. After that
we simulated 500 generations to obtain linkage disequilibrium using 
500 males and 500 females that were mated at random. 
Random mating was continued for 5 more generations and 
population size was increased to 4000 males and 4000 females. 
The QTL on chromosome 1 has effect on trait 1, wherehas 
chomosomes 1 and 2 had effects on traits 2 and 3 respectively.
The QTL on chromosome 4 had effects on the 3 traits. The effects 
of QTL on chromosome 4 were simulated from a multivariate
normal distribution with null mean and variance covariance matrix:

\preformatted{
1.00	 0.75  0.50
0.75	 1.00  0.75  
0.50	 0.75  1.00
 } 

The genetic values were scaled to have variance 1.0.
The phenotypes  for these traits were obtained by 
adding residuals to genetic values, residuals were 
simulated from a multivariate normal distribution 
with null mean and covariance matrix: 
 
\preformatted{
6.0  6.0  1.0
6.0  8.0  2.0  
1.0  2.0  1.0
 } 
 
 In total, 8000 individuals were simulated and the genetic covariance 
 matrix is:

\preformatted{
1.00  0.34  0.07
0.34  1.00  0.21 
0.07  0.21  1.00 
}

}

\usage{
  data(simulated3t)
}

\format{

 Matrix simulated3t.X contains the marker information.
 Matrix simulated3t.pheno contains the phenotypical information.
 
}

\references{
  
  Cheng, H.,  K. Kadir, J., Zeng, D. Garrick and R. Fernando. 2018. Genomic Prediction from Multiple-Trait 
  Bayesian Regression Methods Using Mixture Priors. \emph{Genetics}, \bold{209(1)}: 89-103.
  
}

\keyword{datasets}
