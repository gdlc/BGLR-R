\name{wheat.X}
\alias{wheat.X}
\title{Molecular markers}
\description{

 Is a matrix (599 x 1279) with DArT genotypes; data are
from pure lines and genotypes were coded as 0/1 denoting
the absence/presence of the DArT. Markers with a
minor allele frequency lower than 0.05 were removed, and
missing genotypes were imputed with samples from the
marginal distribution of marker genotypes, that is, \eqn{x_{ij}=Bernoulli(\hat p_j)}, where  \eqn{\hat p_j}  
is the estimated allele frequency computed from the non-missing genotypes. The number of DArT 
MMs after edition was 1279.

}

\source{
  International Maize and Wheat Improvement Center (CIMMYT), Mexico.
}

\keyword{datasets}
