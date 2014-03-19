\name{medianFilter}
\alias{medianFilter}
\title{Apply a median filter to a matrix}
\description{
  For each element in a matix, replace it with the median of the
  values around it.
}
\usage{
medianFilter(x, mrad, nrad)
}
\arguments{
  \item{x}{numeric matrix to median filter}
  \item{mrad}{
    number of rows on either side of the value to use for median
    calculation
  }
  \item{nrad}{
    number of rows on either side of the value to use for median
    calculation
  }
}
\value{
  A matrix whose values have been median filtered
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\examples{
mat <- matrix(1:25, nrow=5)
mat
medianFilter(mat, 1, 1)
}
\keyword{manip}
