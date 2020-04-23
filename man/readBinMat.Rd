\name{readBinMat}
\alias{readBinMat}
\title{readBinMat}
\description{
  Function to read effects saved by BGLR when ETA[[j]]$saveEffects=TRUE.
}
\usage{

	readBinMat(filename,byrow=TRUE)

}
\arguments{

	\item{filename}{(string), the name of the file to be read.}
	\item{byrow}{(logical), if TRUE the matrix is created by filling its corresponding elements by rows.}
}

\value{
	A matrix with samples of regression coefficients.	
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

}

}
\keyword{models}
