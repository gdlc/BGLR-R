\name{BGLR}
\alias{BGLR}
\title{Bayesian Generalized Linear Regression}
\description{
  The BGLR (`Bayesian Generalized Linear Regression') function
  fits various types of parametric and semi-parametric Bayesian regressions to continuos 
  (censored or not), binary and ordinal outcomes.
}
\usage{
  BGLR(y, response_type = "gaussian", a=NULL, b=NULL,ETA = NULL, nIter = 1500,
       burnIn = 500, thin = 5, saveAt = "", S0 = NULL, 
       df0 =5, R2 = 0.5, weights = NULL,
       verbose = TRUE, rmExistingFiles = TRUE, groups=NULL)

}
\arguments{
   \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs allowed).}
   \item{response_type}{(string) admits values \code{"gaussian"} or \code{"ordinal"}. The Gaussian outcome may be censored or not (see below). 
                        If \code{response_type="gaussian"}, \code{y} should be coercible to numeric. If \code{response_type="ordinal"}, \code{y}
                         should be coercible to character, and the order of the outcomes is determined based on the
                        alphanumeric order \code{(0<1<2..<a<b..)}. For ordinal traits the probit link is used.}
   \item{a,b}{(numeric, \eqn{n}{n}) only requiered for censored outcomes, \code{a} and \code{b} are vectors specifying lower and 
              upper bounds for censored observations, respectively. The default value, for non-censored and ordinal 
              outcomes, is \code{NULL} (see details).}
   \item{ETA}{ (list) This is a two-level list used to specify the regression function (or linear predictor). By default the linear predictor 
              (the conditional expectation function in case of Gaussian outcomes) includes only an intercept.
              Regression on covariates and other types of random effects are specified in this two-level list. For instance:

\preformatted{
ETA=list(list(X=W, model="FIXED"), 
              list(X=Z,model="BL"), 
              list(K=G,model="RKHS")),
} 

	      specifies that the linear predictor should include: an intercept (included by default) plus a linear regression on W with regression 
              coefficients treated as fixed effects (i.e., flat prior), plus regression on Z, with regression coefficients modeled as in the 
              Bayesian Lasso of Park and Casella (2008) plus and a random effect with co-variance structure G.


	      For linear regressions the following options are implemented: FIXED (Flat prior), BRR (Gaussian prior), BayesA (scaled-t prior), BL (Double-Exponential prior), 
              BayesB (two component mixture prior with a point of mass at zero and a sclaed-t slab), BayesC (two component mixture prior with a point of
              mass at zero and a Gaussian slab). In linear regressions X can be the incidence matrix for effects or a formula (e.g. \code{~factor(sex) + age})
              in which case the incidence matrix is created internally using the \code{model.matrix} function of R. For Gaussian processes (RKHS) a co-variance matrix
              (K) must be provided. Further details about the models implemented in BGLR see the vignettes in the package. 
             }
   \item{weights}{(numeric, \eqn{n}{n}) a vector of weights, may be \code{NULL}. If weights is not \code{NULL}, the 
                  residual variance of each data-point is set to be proportional to the inverse of the squared-weight. Only
                  used with Gaussian outcomes.}
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
   \item{groups}{(factor) a vector of the same length of y that associates observations with groups, each group will have an associated variance component for the error term.}
}

\details{

BGLR implements a Gibbs sampler for a Bayesian regresion model. The linear predictor (or regression function)
includes an intercept (introduced by default) plus a number of user-specified regression components (X) and 
random effects (u), that is:

  \ifelse{html}{\out{&eta;=1&mu; + X<sub>1</sub>&beta;<sub>1</sub>+...+X<sub>p</sub>&beta;<sub>p</sub>+u<sub>1</sub>+...+u<sub>q</sub>}}{\deqn{\eta=1\mu +X_1 \beta_1+...+X_p \beta_p +u_1+...+u_q}}
  
The components of the linear predictor are specified in the argument ETA (see above). 
The user can specify as many linear terms as desired, and for each component the user can choose the 
prior density to be assigned. The distribution of the response is modeled as a function of the 
linear predictor.

For Gaussian outcomes, the linear predictor is the conditional expectation, and censoring is allowed. For censored data points the actual response value (\eqn{y_i}{y_i})
is missing, and the entries of the vectors a and b (see above) give the lower an upper vound for \eqn{y_i}{y_i}. The
following table shows the configuration of the triplet (y, a, b) for un-censored, right-censored, 
left-censored and interval censored.

\tabular{lccc}{
                     \tab a             \tab y          \tab b  \cr
   Un-censored       \tab NULL          \tab \eqn{y_i}  \tab NULL \cr
   Right censored    \tab \eqn{a_i}     \tab NA         \tab \eqn{\infty} \cr
   Left censored     \tab \eqn{-\infty} \tab NA         \tab \eqn{b_i} \cr
   Interval censored \tab \eqn{a_i}     \tab NA         \tab \eqn{b_i}
 } 

Internally, censoring is dealt with as a missing data problem.

\emph{Ordinal outcomes} are modelled using the probit link, implemented via data augmentation. 
In this case the linear predictor becomes the mean of the underlying liability variable which is 
normal with mean equal to the linear predictor and variance equal to one. In case of only two classes 
(binary outcome) the threshold is set equal to zero, for more than two classess thresholds are estimated 
from the data. Further details about this approach can be found in Albert and Chib (1993).

}

\value{

A list with estimated posterior means, estimated posterior standard deviations, and the parameters used to fit the model.
See the vignettes in the package for further details.

}

\references{

Albert J,. S. Chib. 1993. Bayesian Analysis of Binary and Polychotomus Response Data. \emph{JASA}, \bold{88}: 669-679.

de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

de los Campos, G., D. Gianola, G. J. M., Rosa, K. A., Weigel, and J. Crossa. 2010.  Semi-parametric genomic-enabled prediction of genetic values using 
reproducing kernel Hilbert spaces methods. \emph{Genetics Research}, \bold{92}:295-308.

Park T. and G. Casella. 2008. The Bayesian LASSO. \emph{Journal of the American Statistical Association} \bold{103}: 681-686.

Spiegelhalter, D.J., N.G. Best, B.P. Carlin and A. van der Linde. 2002. Bayesian measures of model complexity and 
fit (with discussion). \emph{Journal of the Royal Statistical Society}, Series B (Statistical Methodology) \bold{64} (4): 583-639.
}

\author{
Gustavo de los Campos, Paulino Perez Rodriguez,
}
\examples{

\dontrun{
#Demos
library(BGLR)

#BayesA
demo(BA)

#BayesB
demo(BB)

#Bayesian LASSO
demo(BL)

#Bayesian Ridge Regression
demo(BRR)

#BayesCpi
demo(BayesCpi)

#RKHS
demo(RKHS)

#Binary traits
demo(Bernoulli)

#Ordinal traits
demo(ordinal)

#Censored traits
demo(censored)

}

}
\keyword{models}
