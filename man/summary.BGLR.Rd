\name{summary.BGLR}
\alias{summary.BGLR}
\title{summary for BGLR fitted models}
\description{
    Gives a summary for a fitted model using BGLR function.
}
\usage{
\method{summary}{BGLR}(object, ...)
}
\arguments{
  \item{object}{An object of class \code{BGLR}.}
  \item{...}{Further arguments passed to or from other methods.}
}
\author{
  Gustavo de los Campos, Paulino Perez Rodriguez,
}
\seealso{
  \code{BGLR}.
}
\examples{

\dontrun{

setwd(tempdir())
library(BGLR)
data(wheat)
out=BLR(y=wheat.Y[,1],XL=wheat.X)
summary(out)

}

}
\keyword{summary}
