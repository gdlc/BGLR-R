\name{BLR}
\alias{BLR}
\title{Bayesian Linear Regression}
\description{
  The BLR (`Bayesian Linear Regression') function
  was designed to fit parametric regression models using different
  types of shrinkage methods. An earlier version of this program was presented in de los Campos \emph{et al.} (2009).
}
\usage{
  BLR(y, XF, XR, XL, GF, prior, nIter, burnIn, thin,thin2,saveAt,
      minAbsBeta,weights)
}
\arguments{
  \item{y}{(numeric, \eqn{n}{n}) the data-vector (NAs allowed).}
  \item{XF}{(numeric, \eqn{n \times pF}{n x pF}) incidence matrix for \eqn{\boldsymbol \beta_F}{bF}, may be NULL.}
  \item{XR}{(numeric, \eqn{n \times pR}{n x pR}) incidence matrix for \eqn{\boldsymbol \beta_R}{bR}, may be NULL.}
  \item{XL}{(numeric, \eqn{n \times pL}{n x pL}) incidence matrix for \eqn{\boldsymbol \beta_L}{bL}, may be NULL.}
  \item{GF}{(list) providing an \eqn{\$}{$}ID (integer, \eqn{n}{n}) linking observations to groups
            (e.g., lines or sires) and a (co)variance structure (\eqn{\$}{$}A, numeric, \eqn{pU \times pU}{pU x pU}) between effects of the grouping factor
            (e.g., line or sire effects). Note: ID must be an integer taking values from 1 to \eqn{pU}{pU}; ID[i]=\eqn{q}{q} indicates that
            the ith observation in \eqn{\boldsymbol y}{y} belongs to cluster \eqn{q}{q} whose (co)variance function is in the qth row (column) of \eqn{\boldsymbol A}{A}.
            GF may be NULL.}
  \item{weights}{(numeric, \eqn{n}{n}) a vector of weights, may be NULL.}
  \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
  \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
  \item{prior}{(list) containing the following elements,
                \itemize{
                  \item prior$varE, prior$varBR, prior$varU: (list) each providing degree of freedom ($df) and scale ($S). These are the parameters of the scaled inverse-\eqn{\chi^2}
                  distributions assigned to variance components, see Eq. (2) below. In the parameterization used by BLR() the prior expectation of variance parameters is \eqn{S/(df-2)}{S/(df-2)}.
                  \item prior$lambda: (list) providing $value (initial value for \eqn{\lambda}{lambda}); $type (`random' or `fixed') this argument specifies
                        whether \eqn{\lambda}{lambda} should be kept fixed at the value provided by $value or updated with samples from the posterior
                        distribution; and, either $shape and $rate (this when a Gamma prior is desired on \eqn{\lambda^2}{lambda^2}) or $shape1, $shape2 and
                        $max, in this case \eqn{p(\lambda |\max, \alpha_1, \alpha_2) \propto Beta \left(\frac{\lambda}{\max} | \alpha_1, \alpha_2 \right)}{p(lambda|max,alpha1,alpha2) = Beta(lamda/max |alpha1,alpha2)}. For detailed description of these priors see de los Campos \emph{et al.} (2009).
                }
              }
  \item{thin2}{This value controls wether the running means are saved to disk or not. If thin2 is greater than nIter the running 
               means are not saved (default, thin2=\eqn{1 \times 10^{10}}{1e10}).}
  \item{minAbsBeta}{The minimum absolute value of the components of \eqn{\boldsymbol \beta_L}{bL} to avoid numeric problems when sampling from \eqn{\boldsymbol \tau^2}, default \eqn{1 \times 10^{-9}}{1e-9} }
}
\details{

The program runs a Gibbs sampler for the Bayesian regression model described below.

\bold{Likelihood}. The equation for the data is:

\deqn{ 
 \begin{array}{lr}
 \boldsymbol y= \boldsymbol 1 \mu + \boldsymbol X_F \boldsymbol \beta_F + \boldsymbol X_R \boldsymbol \beta_R + \boldsymbol X_L \boldsymbol \beta_L + \boldsymbol{Zu}  + \boldsymbol \varepsilon & (1) 
 \end{array}
     }{y=1*mu + XF*bF + XR*bR + XL*bL + Z*u + e ...(1)}

where \eqn{\boldsymbol y}{y}, the response is a \eqn{n \times 1}{n x 1} vector (NAs allowed); \eqn{\mu}{mu} is
an intercept; \eqn{\boldsymbol X_F, \boldsymbol X_R, \boldsymbol X_L}{XF, XR, XL} and \eqn{\boldsymbol Z}{Z} are incidence matrices 
used to accommodate different 
types of effects (see below), and; \eqn{\boldsymbol \varepsilon}{e} is a vector of model residuals assumed to be 
distributed as \eqn{\boldsymbol \varepsilon \sim N(\boldsymbol 0,Diag(\sigma_{\boldsymbol \varepsilon}^2/w_i^2))}{e ~ MVN(0,Diag(varE/wi^2))}, 
here \eqn{\sigma_{\boldsymbol \varepsilon}^2}{varE} is an (unknown) 
variance parameter and \eqn{w_i} are (known) weights that allow for heterogeneous-residual variances.

Any of the elements in the right-hand side of the linear predictor, except \eqn{\mu}{mu} and \eqn{\boldsymbol \varepsilon}{e} , can be omitted; 
by default the program runs an intercept model.

\bold{Prior}. The residual variance is assigned a scaled inverse-\eqn{\chi^2} prior with degree of freedom and scale parameter
provided by the user, that is, \eqn{\sigma_{\boldsymbol \varepsilon}^2 \sim \chi^{-2} (\sigma_{\boldsymbol \varepsilon}^2 | df_{\boldsymbol \varepsilon}, S_{\boldsymbol \varepsilon})}{varE ~ Inv_Scaled_chisq(varE | dfe,Se)}. The regression coefficients \eqn{\left\{\mu, \boldsymbol \beta_F, \boldsymbol \beta_R, \boldsymbol \beta_L, \boldsymbol u \right\}}{{mu,bF,bR,bL,u}} are assigned priors
that yield different type of shrinkage. The intercept and the vector of regression coefficients \eqn{\boldsymbol \beta_F}{bF} are assigned flat priors
(i.e., estimates are not shrunk). The vector of regression coefficients \eqn{\boldsymbol \beta_R}{bR} is assigned a
Gaussian prior with variance common to all effects, that is,
\eqn{\beta_{R,j} \mathop \sim \limits^{iid} N(0, \sigma_{\boldsymbol \beta_R}^2)}{bRj ~ NIID(0,varBr)}. This prior is 
the Bayesian counterpart of Ridge Regression. The variance parameter \eqn{\sigma_{\boldsymbol \beta_R}^2}{varBr}, 
is treated as unknown and it is assigned a scaled inverse-\eqn{\chi^2} prior, that is,
\eqn{\sigma_{\boldsymbol \beta_R}^2 \sim \chi^{-2} (\sigma_{\boldsymbol \beta_R}^2 | df_{\boldsymbol \beta_R}, S_{\boldsymbol \beta_R})}{varBr ~ Inv_Scaled_chisq(varBr | dfBr,SBr)} with degrees of freedom 
\eqn{df_{\boldsymbol \beta_R}}{dfBr}, and scale \eqn{S_{\boldsymbol \beta_R}}{SBr} provided by the user.

The vector of regression coefficients \eqn{\boldsymbol \beta_L}{bL} is treated as in
the Bayesian LASSO of Park and Casella (2008). Specifically,

\deqn{p(\boldsymbol \beta_L, \boldsymbol \tau^2, \lambda | \sigma_{\boldsymbol \varepsilon}^2) = \left\{\prod_k N(\beta_{L,k} | 0, \sigma_{\boldsymbol \varepsilon}^2 \tau_k^2) Exp\left(\tau_k^2 | \lambda^2\right) \right\} p(\lambda),}{p(bL,tau^2,lambda) | varE) = prod(N(bL_k|0,varE*tau_k^2)*Exp(tau_k^2 | lamda^2), k) p(lambda),}

where, \eqn{Exp(\cdot|\cdot)}{Exp(.|.)} is an exponential prior and \eqn{p(\lambda)}{p(lambda)} can either be: (a)
a mass-point at some value (i.e., fixed \eqn{\lambda}{lambda}); (b) \eqn{p(\lambda^2) \sim Gamma(r,\delta)}{p(lambda^2)~Gamma(r,delta)} this 
is the prior suggested by Park and Casella (2008); or, (c) \eqn{p(\lambda | \max, \alpha_1, \alpha_2) \propto Beta\left(\frac{\lambda}{\max} | \alpha_1,\alpha_2 \right)}{p(lambda|max,alpha1,alpha2) = Beta(lamda/max |alpha1,alpha2)}, see de los Campos \emph{et al}. (2009) for details. It can be shown that the marginal prior of regression coefficients \eqn{\beta_{L,k}, \int N(\beta_{L,k} | 0, \sigma_{\boldsymbol \varepsilon}^2 \tau_k^2) Exp\left(\tau_k^2 | \lambda^2\right) \partial \tau_k^2}{bL_k, Integrate(N(bL_k|0,varE * tau_k^2) * Exp(tau_k^2 | lambda^2) d tau_k^2)}, is Double-Exponential. This prior has thicker tails and higher peak of mass at zero than the Gaussian prior used for \eqn{\boldsymbol \beta_R}{bR}, inducing a different type of shrinkage.

The vector \eqn{\boldsymbol u }{u} is used to model the so called `infinitesimal effects', and is assigned a prior \eqn{\boldsymbol u \sim N(\boldsymbol 0, \boldsymbol A\sigma_{\boldsymbol u}^2)}{u~MVN(0,varU)}, 
where, \eqn{\boldsymbol A}{A} is a positive-definite matrix (usually a relationship matrix computed from a pedigree) and \eqn{\sigma_{\boldsymbol u}^2}{varU} is an unknow variance, whose prior is
\eqn{\sigma_{\boldsymbol u}^2 \sim \chi^{-2} (\sigma_{\boldsymbol u}^2 | df_{\boldsymbol u}, S_{\boldsymbol u})}{varU ~ Inv_Scaled_chisq(varU | dfu,Su)}.

Collecting the above mentioned assumptions, the posterior distribution of model unknowns, 
\eqn{\boldsymbol \theta= \left\{\mu, \boldsymbol \beta_F, \boldsymbol \beta_R, \sigma_{\boldsymbol \beta_R}^2, \boldsymbol \beta_L, \boldsymbol \tau^2, \lambda, \boldsymbol u, \sigma_{\boldsymbol u}^2, \sigma_{\boldsymbol \varepsilon}^2, \right\}}{theta={mu, bF, bR, varBr, bL, tau^2, lambda, u, varU, varE}}, is, 

\deqn{\begin{array}{lclr}
        p(\boldsymbol \theta | \boldsymbol y) & \propto & N\left( \boldsymbol y | \boldsymbol 1 \mu + \boldsymbol X_F \boldsymbol \beta_F + \boldsymbol X_R \boldsymbol \beta_R + \boldsymbol X_L \boldsymbol \beta_L + \boldsymbol{Zu}; Diag\left\{ \frac{\sigma_\varepsilon^2}{w_i^2}\right\}\right) & \\
                                              &         & \times  \left\{ \prod\limits_j N\left(\beta_{R,j} | 0, \sigma_{\boldsymbol \beta_R}^2 \right) \right\} \chi^{-2} \left(\sigma_{\boldsymbol \beta_R}^2  | df_{\boldsymbol \beta_R}, S_{\boldsymbol \beta_R}\right)  & \\
                                              &         & \times \left\{ \prod\limits_k N \left( \beta_{L,k} |0,\sigma_{\boldsymbol \varepsilon}^2 \tau_k^2 \right) Exp \left(\tau_k^2 | \lambda^2\right)\right\} p(\lambda) & (2)\\
                                              &         & \times N(\boldsymbol u | \boldsymbol 0,\boldsymbol A\sigma_{\boldsymbol u}^2) \chi^{-2} (\sigma_{\boldsymbol u}^2 | df_{\boldsymbol u}, S_{\boldsymbol u}) \chi^{-2} (\sigma_{\boldsymbol \varepsilon}^2 | df_{\boldsymbol \varepsilon}, S_{\boldsymbol \varepsilon}) & 
      \end{array}
     }{p(theta|y)=N(y|1*mu + XF*bF + XR*bR + XL*bL + Z*u) 
                  * prod(N(bR_j|0,varBr),j) * Inv_Scaled_chisq(varBr | dfBr,SBr)
                  * prod(N(bL_k|0,varE * tau_k^2)* Exp(tau_k^2 | lambda^2),k) * p(lambda)
                  * MVN(u|0,varU) * Inv_Scaled_chisq(varU | dfu,Su) * Inv_Scaled_chisq(varE | dfe,Se) ...(2)
      }
}
\value{
A list with posterior means, posterior standard deviations, and the parameters used to fit the model:
\item{$yHat}{the posterior mean of \eqn{\boldsymbol 1 \mu + \boldsymbol X_F \boldsymbol \beta_F + \boldsymbol X_R \boldsymbol \beta_R + \boldsymbol X_L \boldsymbol \beta_L + \boldsymbol{Zu}  + \boldsymbol\varepsilon}{1*mu + XF*bF + XR*bR + XL*bL + Z*u}.}
\item{$SD.yHat}{the corresponding posterior standard deviation.}
\item{$mu}{the posterior mean of the intercept.}
\item{$varE}{the posterior mean of \eqn{\sigma_{\boldsymbol \varepsilon}^2}{varE}.}
\item{$bR}{the posterior mean of \eqn{\boldsymbol \beta_R}{bR}.}
\item{$SD.bR}{the corresponding posterior standard deviation.}
\item{$varBr}{the posterior mean of \eqn{\sigma_{\boldsymbol \beta_R}^2}{varBr}.}
\item{$bL}{the posterior mean of \eqn{\boldsymbol \beta_L}{bL}.}
\item{$SD.bL}{the corresponding posterior standard deviation.}
\item{$tau2}{the posterior mean of \eqn{\boldsymbol \tau^2}{tau^2}.}
\item{$lambda}{the posterior mean of \eqn{\lambda}{lambda}.}
\item{$u}{the posterior mean of \eqn{\boldsymbol u}{u}.}
\item{$SD.u}{the corresponding posterior standard deviation.}
\item{$varU}{the posterior mean of \eqn{\sigma_{\boldsymbol u}^2}{varU}.}
\item{$fit}{a list with evaluations of effective number of parameters and DIC (Spiegelhalter \emph{et al.}, 2002).}
\item{$whichNa}{a vector indicating which entries in \eqn{\boldsymbol y} were missing.}
\item{$prior}{a list containig the priors used during the analysis.}
\item{$weights}{vector of weights.}
\item{$fit}{list containing the following elements,
            \itemize{
             \item $logLikAtPostMean: log-likelihood evaluated at posterior mean.
             \item $postMeanLogLik: the posterior mean of the Log-Likelihood.
             \item $pD: estimated effective number of parameters, Spiegelhalter \emph{et al.} (2002).
             \item $DIC: the deviance information criterion, Spiegelhalter \emph{et al.} (2002).                        
            }}
\item{$nIter}{the number of iterations made in the Gibbs sampler.}
\item{$burnIn}{the nuber of iteratios used as burn-in.}
\item{$thin}{the thin used.}
\item{$y}{original data-vector.}

The posterior means returned by BLR are calculated after burnIn is
passed and at a thin as specified by the user.

\bold{Save}. The routine will save samples of \eqn{\mu}{mu}, variance components and \eqn{\lambda}{lambda} and running means
(rm*.dat). Running means are computed using the thinning specified by
the user (see argument thin above); however these running means are
saved at a thinning specified by argument thin2 (by default, thin2=\eqn{1 \times 10^{10}}{1e10} so that running means are computed 
as the sampler runs but not saved to the disc).

}
\references{
de los Campos G., H. Naya, D. Gianola, J. Crossa, A. Legarra, E. Manfredi, K. Weigel and J. Cotes. 2009.
Predicting Quantitative Traits with Regression Models for Dense Molecular Markers and Pedigree. \emph{Genetics} \bold{182}: 375-385.

Park T. and G. Casella. 2008. The Bayesian LASSO. \emph{Journal of the American Statistical Association} \bold{103}: 681-686.

Spiegelhalter, D.J., N.G. Best, B.P. Carlin and A. van der Linde. 2002. Bayesian measures of model complexity and 
fit (with discussion). \emph{Journal of the Royal Statistical Society}, Series B (Statistical Methodology) \bold{64} (4): 583-639.
}
\author{
Gustavo de los Campos, Paulino Perez Rodriguez,
}
\examples{
\dontrun{
########################################################################
##Example 1:
########################################################################

rm(list=ls())
setwd(tempdir())
library(BGLR)
data(wheat)     #Loads the wheat dataset

y=wheat.Y[,1]
### Creates a testing set with 100 observations
whichNa<-sample(1:length(y),size=100,replace=FALSE)
yNa<-y
yNa[whichNa]<-NA

### Runs the Gibbs sampler
fm<-BLR(y=yNa,XL=wheat.X,GF=list(ID=1:nrow(wheat.A),A=wheat.A),
                           prior=list(varE=list(df=3,S=0.25),
                           varU=list(df=3,S=0.63),
                           lambda=list(shape=0.52,rate=1e-4,
                           type='random',value=30)),
                           nIter=5500,burnIn=500,thin=1)

MSE.tst<-mean((fm$yHat[whichNa]-y[whichNa])^2)
MSE.tst
MSE.trn<-mean((fm$yHat[-whichNa]-y[-whichNa])^2)
MSE.trn
COR.tst<-cor(fm$yHat[whichNa],y[whichNa])
COR.tst
COR.trn<-cor(fm$yHat[-whichNa],y[-whichNa])
COR.trn

plot(fm$yHat~y,xlab="Phenotype",
     ylab="Pred. Gen. Value" ,cex=.8)
points(x=y[whichNa],y=fm$yHat[whichNa],col=2,cex=.8,pch=19)

x11()
plot(scan('varE.dat'),type="o",
        ylab=expression(paste(sigma[epsilon]^2)))

########################################################################
#Example 2: Ten fold, Cross validation, environment 1,
########################################################################

rm(list=ls())
setwd(tempdir())
library(BGLR)
data(wheat)     #Loads the wheat dataset
nIter<-1500     #For real data sets more samples are needed
burnIn<-500     
thin<-10
folds<-10
y<-wheat.Y[,1]
A<-wheat.A

priorBL<-list(
               varE=list(df=3,S=2.5),
               varU=list(df=3,S=0.63),
               lambda = list(shape=0.52,rate=1e-5,value=20,type='random')
             )
             
set.seed(123)  #Set seed for the random number generator
sets<-rep(1:10,60)[-1]
sets<-sets[order(runif(nrow(A)))]
COR.CV<-rep(NA,times=(folds+1))
names(COR.CV)<-c(paste('fold=',1:folds,sep=''),'Pooled')
w<-rep(1/nrow(A),folds) ## weights for pooled correlations and MSE
yHatCV<-numeric()

for(fold in 1:folds)
{
   yNa<-y
   whichNa<-which(sets==fold)
   yNa[whichNa]<-NA
   prefix<-paste('PM_BL','_fold_',fold,'_',sep='')
   fm<-BLR(y=yNa,XL=wheat.X,GF=list(ID=(1:nrow(wheat.A)),A=wheat.A),prior=priorBL,
               nIter=nIter,burnIn=burnIn,thin=thin)
   yHatCV[whichNa]<-fm$yHat[fm$whichNa]
   w[fold]<-w[fold]*length(fm$whichNa)
   COR.CV[fold]<-cor(fm$yHat[fm$whichNa],y[whichNa])
}

COR.CV[11]<-mean(COR.CV[1:10])
COR.CV

########################################################################
}
}
\keyword{models}
