%Function included in the MCMCpack, version 1.4-5

\name{vech}
\alias{vech}
\title{Extract Lower Triangular Elements from a Symmetric Matrix}
\usage{
vech(x)
}
\arguments{
\item{x}{A symmetric matrix.}
}
\value{
A list of the lower triangular elements.
}
\description{
This function takes a symmetric matrix and extracts a list of all lower
triangular elements.
}
\details{
This function checks to make sure the matrix is square, but it does not
check for symmetry (it just pulls the lower triangular elements).  The
elements are stored in column major order.  The original matrix can be
restored using the \code{xpnd} command.
}
\examples{

   symmat <- matrix(c(1,2,3,4,2,4,5,6,3,5,7,8,4,6,8,9),4,4)
   vech(symmat)

}
\seealso{
\code{\link{xpnd}}
}
\keyword{manip}
