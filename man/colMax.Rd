\name{colMax}
\alias{colMax}
\alias{rowMax}
\title{Find row and column maximum values}
\description{
  Find row and column maximum values for numeric arrays.
}
\usage{
colMax(x, na.rm = FALSE, dims = 1)
rowMax(x, na.rm = FALSE, dims = 1)
}
\arguments{
  \item{x}{
    an array of two or more dimensions, containing numeric values
  }
  \item{na.rm}{
    logical. Should missing values (including 'NaN') be omitted
    from the calculations? (not currently implemented)
  }
  \item{dims}{
    Which dimensions are regarded as "rows" or "columns" to maximize.
    For \code{rowMax}, the maximum is over dimensions \code{dims+1,
    ...}; for \code{colMax} it is over dimensions \code{1:dims}.
  }
}
\details{
  These functions are designed to act like the \code{colSums} series
  of functions except that they only currently handle real arrays
  and will not  remove \code{NA} values.
}
\value{
  A numeric array of suitable size, or a vector if the result is
  one-dimensional. The \code{dimnames} (or \code{names} for a
  vector result) are taken from the original array.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{colSums}}
}
\keyword{array}
\keyword{manip}
