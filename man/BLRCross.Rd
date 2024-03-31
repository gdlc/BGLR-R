\name{BLRCross}
\alias{BLRCross}
\title{Bayesian Linear Regression}
\description{
  The BLR (`Bayesian Linear Regression') function
  fits various types of parametric Bayesian regressions to continuos outcomes.
}

\usage{

BLRCross(y=NULL,my=NULL,vy=NULL,n=NULL,
         XX,Xy,nIter=1500,burnIn=500,
         thin=5,R2=0.5,
         S0=NULL,df0=5,
         priors=NULL,
         idPriors=NULL,
         verbose=TRUE,
         saveAt = "", 
         rmExistingFiles=TRUE)
}

\arguments{
   \item{y}{A numeric vector of length n, NAs not allowed, if NULL you must provide my, vy and n.}
   \item{my}{numeric, sample mean of y. If NULL you must provide y.}
   \item{vy}{numeric, sample variance of y. If NULL you must provide y.}
   \item{n}{integer, sample size. If NULL you must provide y.}
   \item{XX}{A matrix, XX=crossprod(X), with X an incidence matrix of dimension n times p.}
   \item{Xy}{A numeric vector of length p, Xy=crossprod(X,y).}
   \item{nIter,burnIn, thin}{(integer) the number of iterations, burn-in and thinning.}
   \item{R2}{(numeric, \code{0<R2<1}) The proportion of variance that one expects, a priori, to be explained by the regression. Only used if
             the hyper-parameters are not specified; if that is the case, internally, hyper-parameters are set so that the prior modes 
             are consistent with the variance partition specified by R2 and the prior distribution is relatively flat at the mode.}
   \item{S0, df0}{(numeric) The scale parameter for the scaled inverse-chi squared prior assigned to the residual variance.
                  In the parameterization of the scaled-inverse chi square in BLRCross the expected values is \code{S0/(df0-2)}. The default value for the 
                  df parameter is 5. If the scale is not specified a value is calculated so that the prior mode of the residual 
                  variance  equals \code{var(y)*R2}.}
   \item{priors}{(list) This is a two-level list used to specify the regression function (or linear predictor). 
              Regression on covariates and other types of random effects are specified in this two-level list.

	      For linear regressions the following options are implemented: FIXED (flat prior),
	      BayesA, BayesB, BRR (Gaussian prior), BayesC, SSVS and RKHS. 
             }
   \item{idPriors}{(numeric) an integer vector that allow us to specify the priors for the columns of matrix X.}
   \item{verbose}{(logical) if TRUE the iteration history is printed, default TRUE.}
   \item{saveAt}{(string) this may include a path and a pre-fix that will be added to the name of the files that are saved as the program runs.}
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
y=y-mean(y)

XX=crossprod(X)
Xy=as.vector(crossprod(X,y))

#Example 1

## BayesA ################################################################################
priors=list(list(model="BayesA"))
idPriors=rep(1,p)
fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)

## BayesB ################################################################################
priors=list(list(model="BayesB"))
idPriors=idPriors=rep(1,p)
fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)
plot(fm1$ETA[[1]]$b)
points(QTL,b[QTL],col="red")

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)
plot(fm1s$ETA[[1]]$b)
points(QTL,b[QTL],col="red")

## BayesC ################################################################################
priors=list(list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=1/100,counts=1E6))
idPriors=rep(1,p)

fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)
plot(fm1$ETA[[1]]$b)
points(QTL,b[QTL],col="red")

plot(fm1$ETA[[1]]$d)
points(QTL,b[QTL],col="red")

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)
plot(fm1s$ETA[[1]]$b)
points(QTL,b[QTL],col="red")


## SSVS (Absolutely Continuous Spike Slab) ###############################################
priors=list(list(model="SSVS",R2=NULL,df0=NULL,S0=NULL,probIn=NULL,counts=NULL))
idPriors=rep(1,p)

fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)
plot(fm1$ETA[[1]]$b)

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)
plot(fm1s$ETA[[1]]$b)

priors=list(list(model="SSVS",R2=NULL,df0=NULL,S0=NULL,probIn=NULL,
                      counts=NULL,cprobIn=0.5,ccounts=2))
idPriors=rep(1,p)

fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)
plot(fm1$ETA[[1]]$b)

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)

plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)
plot(fm1s$ETA[[1]]$b)

## Ridge Regression ######################################################################
priors=list(list(model="BRR",R2=NULL,df0=NULL,S0=NULL))
idPriors=rep(1,p)

fm1=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1$ETA[[1]]$b),y)

#summary statistics
fm1s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
plot(as.vector(X\%*\%fm1s$ETA[[1]]$b),y)

#Example 2
priors=list(list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=NULL,counts=NULL),
            list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=NULL,counts=NULL))

idPriors=c(rep(1,p/2),rep(2,p/2))

fm2=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
bHat=c(fm2$ETA[[1]]$b,fm2$ETA[[2]]$b)
plot(as.vector(X\%*\%bHat),y)
plot(bHat)

#summary statistics
fm2s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
bHat=c(fm2s$ETA[[1]]$b,fm2s$ETA[[2]]$b)
plot(as.vector(X\%*\%bHat),y)
plot(bHat)

#Example 3
priors=list(list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=NULL,counts=NULL),
            list(model="BRR",R2=NULL,df0=NULL,S0=NULL))
idPriors=c(rep(1,p/2),rep(2,p/2))

fm3=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
bHat=c(fm3$ETA[[1]]$b,fm3$ETA[[2]]$b)
plot(as.vector(X\%*\%bHat),y)
plot(bHat)

#summary statistics
fm3s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,priors=priors,idPriors=idPriors)
bHat=c(fm3s$ETA[[1]]$b,fm3s$ETA[[2]]$b)
plot(as.vector(X\%*\%bHat),y)
plot(bHat)

#Example 4 Fixed effects

y=y+3

X=cbind(1,X)

XX=crossprod(X)
Xy=as.vector(crossprod(X,y))

priors=list(list(model="FIXED"),
            list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=1/100,counts=1e6),
            list(model="BRR",R2=NULL,df0=NULL,S0=NULL))
idPriors=c(1,rep(2,p/2),rep(3,p/2))

fm0=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0$ETA[[1]]$b,fm0$ETA[[2]]$b,fm0$ETA[[3]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

#summary statistics
fm0s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0s$ETA[[1]]$b,fm0s$ETA[[2]]$b,fm0s$ETA[[3]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

priors=list(list(model="FIXED"),
            list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=1/100,counts=1e6),
            list(model="BRR",R2=NULL,df0=NULL,S0=NULL),
            list(model="BayesA"))
idPriors=c(1,rep(2,333),rep(3,333),rep(4,334))

fm0=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0$ETA[[1]]$b,fm0$ETA[[2]]$b,fm0$ETA[[3]]$b,fm0$ETA[[4]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

#summary statistics
fm0s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)

bHat=c(fm0s$ETA[[1]]$b,fm0s$ETA[[2]]$b,fm0s$ETA[[3]]$b,fm0$ETA[[4]]$b)
plot(y,X\%*\%bHat)
plot(bHat)
              

priors=list(list(model="FIXED"),
            list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=1/100,counts=1e6),
            list(model="BRR",R2=NULL,df0=NULL,S0=NULL),
            list(model="BayesA"),
            list(model="BayesB"))
idPriors=c(1,rep(2,250),rep(3,250),rep(4,250),rep(5,250))

fm0=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0$ETA[[1]]$b,fm0$ETA[[2]]$b,fm0$ETA[[3]]$b,fm0$ETA[[4]]$b,fm0$ETA[[5]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

#summary statistics
fm0s=BLRCross(my=mean(y),vy=var(y),n=length(y),
              XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0s$ETA[[1]]$b,fm0s$ETA[[2]]$b,fm0s$ETA[[3]]$b,fm0s$ETA[[4]]$b,fm0s$ETA[[5]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

priors=list(list(model="FIXED"),
            list(model="BayesC",R2=NULL,df0=NULL,S0=NULL,probIn=1/100,counts=1e6),
            list(model="BRR",R2=NULL,df0=NULL,S0=NULL),
            list(model="BayesA"),
            list(model="BayesB"),
            list(model="SSVS"))
idPriors=c(1,rep(2,200),rep(3,200),rep(4,200),rep(5,200),rep(6,200))

fm0=BLRCross(y=y,XX=XX,Xy=Xy,nIter=10000,burnIn=5000,priors=priors,idPriors=idPriors)
bHat=c(fm0$ETA[[1]]$b,fm0$ETA[[2]]$b,fm0$ETA[[3]]$b,fm0$ETA[[4]]$b,fm0$ETA[[5]]$b,fm0$ETA[[6]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

#summary statistics
fm0s=BLRCross(my=mean(y),vy=var(y),n=length(y),XX=XX,Xy=Xy,nIter=10000,burnIn=5000,
              priors=priors,idPriors=idPriors)
bHat=c(fm0s$ETA[[1]]$b,fm0s$ETA[[2]]$b,fm0s$ETA[[3]]$b,fm0s$ETA[[4]]$b,
       fm0s$ETA[[5]]$b,fm0s$ETA[[6]]$b)
plot(y,X\%*\%bHat)
plot(bHat)

}

}
\keyword{models}
