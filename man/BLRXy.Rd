\name{BLRXy}
\alias{BLRXy}
\title{Bayesian Linear Regression}
\description{
  The BLRXy(`Bayesian Linear Regression') function
  fits various types of parametric Bayesian regressions to continuos outcomes.
  This is a wrapper for function BLRCross.
}

\usage{

BLRXy(y, intercept=TRUE, ETA, 
             nIter = 1500, burnIn = 500, thin = 5, 
             S0 = NULL, df0 = 5, R2 = 0.5, 
             verbose = TRUE, saveAt="",
             rmExistingFiles = TRUE)
}

\arguments{
	\item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs allowed).}
	\item{intercept}{(logical) indicates if an intercept is included.}
   \item{ETA}{ (list) This is a two-level list used to specify the regression function (or linear predictor). 
              Regression on covariates and other types of random effects are specified in this two-level list. For instance:

\preformatted{
ETA=list(list(X=W, model="FIXED"), 
              list(X=Z,model="BRR")),
} 

	      specifies that the linear predictor should include: an intercept (included by default), 
	      a linear regression on W with regression coefficients treated as fixed effects (i.e., flat prior), 
	      plus regression on Z, with regression coefficients modeled as in the 
          Bayesian Ridge Regression.

	      The following options are implemented for linear regressions: FIXED (Flat prior),
	      BayesA, BayesB, BRR (Gaussian prior), BayesC and SSVS. 
             }
   \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
   \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
   \item{S0, df0}{ (numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance, only used with Gaussian outcomes. 
                  In the parameterization of the scaled-inverse chi square in BGLR the expected values is \code{S0/(df0-2)}. The default value for the 
                  df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual 
                  variance  equals \code{var(y)*R2} (see below). For further details see the vignettes in the package.}
   \item{R2}{(numeric, \code{0<R2<1}) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if 
             the hyper-parameters are not specified; if that is the case, internaly, hyper-paramters are set so that the prior modes are consistent with the 
             variance partition specified by R2 and the prior distribution is relatively flat at the mode. For 
             further details see the vignettes in the package.}
   \item{verbose}{(logical) if TRUE the iteration history is printed, default TRUE.}
   \item{rmExistingFiles}{(logical) if TRUE removes existing output files from previous runs, default TRUE.}
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

library(BGLR)

p=1000
n=1500

data(mice)
X=scale(mice.X[1:n,1:p],center=TRUE)
QTL=seq(from=50,to=p-50,by=80)

b=rep(0,p)
b[QTL]=1
signal=as.vector(X\%*\%b)
 
error=rnorm(sd=sd(signal),n=n)
y=error+signal
y=2+y

#Example 1

#BayesA
ETA=list(list(X=X,model="BayesA"))

fm1=BLRXy(y=y,ETA=ETA)
plot(fm1$yHat,y)

}

}
\keyword{models}
