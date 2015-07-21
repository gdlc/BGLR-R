\name{plot.BGLR}
\alias{plot.BGLR}
\title{Plots for BGLR Analysis}
\description{
    Plots observed vs predicted values for objects of class BGLR.
}
\usage{
\method{plot}{BGLR}(x, ...)
}
\arguments{
  \item{x}{An object of class \code{BGLR}.}
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
plot(out)

}

}
\keyword{plot}
