\name{getVariances}
\alias{getVariances}
\title{getVariances}
\description{
  Computes the sample-variance (var()) for sets of markers as well as the total variance.
}
\usage{

	getVariances(X, B, sets, verbose=TRUE)

}
\arguments{

	\item{X}{(numeric, \eqn{n \times p}{n x p}) incidence matrix for \eqn{\boldsymbol \beta}{b}.}
	\item{B}{(numeric), object returned by the function  readBinMat().}
	\item{sets}{(numeric).}
	\item{verbose}{(logical), if TRUE it shows progress information in the console.}
}

\value{
	A matrix with variances for markers as well as the total.	
}

\author{
Gustavo de los Campos.
}
\examples{

\dontrun{
#Demos

library(BGLR)
data(wheat)
y=wheat.Y[,1] ; X=scale(wheat.X)
dir.create('test_saveEffects')
setwd('test_saveEffects')
fm=BGLR(y=y,ETA=list(list(X=X,model='BayesB',saveEffects=TRUE)),nIter=12000,thin=2,burnIn=2000)
B=readBinMat('ETA_1_b.bin')
plot(B[,1],type='o',col=4)
VAR=getVariances(B=B,X=X,sets=sample(1:20,size=1279,replace=T))
head(VAR)
plot(VAR[,"total"])

}

}
\keyword{models}
