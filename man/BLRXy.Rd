\name{BLRXy}
\alias{BLRXy}
\title{Bayesian Linear Regression}
\description{
  The BLRXy(`Bayesian Linear Regression') function
  fits various types of parametric Bayesian regressions to continuos outcomes.
}

\usage{

BLRXy(y,
      X,centerX=TRUE,centerY=TRUE,scaleX=FALSE,imputeX=TRUE,
      nIter=1500,burnIn=500,
      thin=5,R2=0.5,
      S0=NULL,df0=5,
      priors=NULL,
      idPriors=NULL,
      verbose=TRUE)
}

\arguments{
   \item{y}{A numeric vector of length n.}
   \item{X}{An incidence matrix of dimension n times p.}
   \item{centerX}{Logical, if TRUE the columns of incidence matrix X are centered.}
   \item{centerY}{Logical, if TRUE the vector y is centered.}
   \item{scaleX}{Logical, if TRUE the columns of X are scaled.}
   \item{imputeX}{Logical, if TRUE missing values in X are imputed.}
   \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
   \item{R2}{(numeric, \code{0<R2<1}) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if
             the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes 
             are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode.}
   \item{S0, df0}{(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. 
                  In the parameterization of the scaled-inverse chi square in BLR the expected values is \code{S0/(df0-2)}. The default value for the 
                  df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual 
                  variance  equals \code{var(y)*R2}.}
   \item{priors}{(list) This is a two-level list used to specify the regression function (or linear predictor). 
              Regression on covariates and other types of random effects are specified in this two-level list.

	      For linear regressions the following options are implemented: BRR (Gaussian prior), BayesA (scaled-t prior), BayesC (two component mixture prior with a point of
              mass at zero and a Gaussian slab). 
             }
   \item{idPriors}{(numeric) an integer vector that allow us to specify the priors for the columns of matrix X.}
   \item{verbose}{(logical) if TRUE the iteration history is printed, default TRUE.}
}

\details{

}

\value{

}

\references{

de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

de los Campos, G., D. Gianola, G. J. M., Rosa, K. A., Weigel, and J. Crossa. 2010.  Semi-parametric genomic-enabled prediction of genetic values using 
reproducing kernel Hilbert spaces methods. \emph{Genetics Research}, \bold{92}:295-308.

}

\author{
Gustavo de los Campos, Paulino Perez Rodriguez.
}
\examples{

\dontrun{

}

}
\keyword{models}