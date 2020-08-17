\name{getGCovar}
\alias{getGCovar}
\title{getGCovar}
\description{
  Genetic co-variance matrix using MCMC samples.
}
\usage{
  getGCovar(X,B)
}
\arguments{
  \item{X}{(numeric) matrix of covariates.}
  \item{B}{(numeric) matrix that contains samples for regression coefficients, 3D array, 
  	            with dim=c(nRow,p,traits), where nRow number of MCMC samples saved, p 
  	            is the number of predictors and traits is the number of traits.}
}
\value{
  Genetic co-variance matrix.
}

\references{
	
Lehermeier, C., G. de los Campos, V. Wimmer and C.-C. Schon. Genomic Variance Estimates: With or without Disequilibrium Covariances?. \emph{J Anim Breed Genet.}, \bold{134} (3): 232-241.

}

\author{
Gustavo de los Campos.
}

