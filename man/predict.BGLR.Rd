\name{predict.BGLR}
\alias{predict.BGLR}
\title{Model Predictions}
\description{
    extracts predictions from the results
    of BGLR function.
}
\usage{
\method{predict}{BGLR}(object, newdata, ...)

}
\arguments{
  \item{object}{An object of class \code{BGLR}.}
  \item{newdata}{Currently not supported, for new data you should assing missing value indicator (NAs) to the corresponding entries in the 
                 response vector (y).}
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
predict(out)

}

}
