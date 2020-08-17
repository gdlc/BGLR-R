% Function taken from MCMCpack version 1.4.5
% 
\name{xpnd}
\alias{xpnd}
\title{Expand a Vector into a Symmetric Matrix}
\usage{
xpnd(x, nrow = NULL)
}
\arguments{
\item{x}{A list of elements to expand into symmetric matrix.}

\item{nrow}{The number of rows (and columns) in the returned matrix.  Look
into the details.}
}
\value{
An \eqn{(nrows \times nrows)} symmetric matrix.
}
\description{
This function takes a vector of appropriate length (typically created using
\code{vech}) and creates a symmetric matrix.
}
\details{
This function is particularly useful when dealing with variance covariance
matrices. Note that R stores matrices in column major order, and that the
items in \code{x} will be recycled to fill the matrix if need be.

The number of rows can be specified or automatically computed from the
number of elements in a given object via \eqn{(-1 + \sqrt{(1 + 8 *
length(x))}) / 2}.
}
\examples{

  xpnd(c(1,2,3,4,4,5,6,7,8,9),4)
  xpnd(c(1,2,3,4,4,5,6,7,8,9))

}
\seealso{
\code{\link{vech}}
}
\keyword{manip}
