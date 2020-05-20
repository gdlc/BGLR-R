\name{readBinMatMultitrait}
\alias{readBinMatMultitrait}
\title{readBinMatMultitrait}
\description{
  Function to read effects saved by Multitrait when ETA[[j]]$saveEffects=TRUE.
}
\usage{

	readBinMatMultitrait(filename,storageMode="double")

}
\arguments{

	\item{filename}{(string), the name of the file to be read.}
	\item{storageMode}{(character), the storage mode used to save effects via ETA[[j]]$storageMode: 'double' (default) or 'single'.}
}

\value{
	A 3D array, with dim=c(nRow,p,traits), where nRow number of MCMC samples saved,
	p is the number of predictors and traits is the number of traits. 
}

\author{
Gustavo de los Campos, Paulino Perez-Rodriguez.
}
\examples{

\dontrun{
library(BGLR)

data(wheat)
y<-wheat.Y 
X<-wheat.X

fm<-Multitrait(y=y,ETA=list(list(X=X,model='BRR',saveEffects=TRUE)),
                      nIter=1000,thin=10,burnIn=500)
                      
 B<-readBinMatMultitrait('ETA_1_beta.bin')

}

}

\keyword{models}
