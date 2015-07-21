\name{residuals.BGLR}
\alias{residuals.BGLR}
\title{Extracts models residuals}
\description{
    extracts model residuals from objects returned by BGLR function.
}
\usage{
\method{residuals}{BGLR}(object, ...)
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
residuals(out)

}

}
\keyword{residuals}
