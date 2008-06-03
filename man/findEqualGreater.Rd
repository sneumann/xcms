\name{findEqualGreater}
\alias{findEqualGreater}
\alias{findEqualLess}
\alias{findEqualGreaterM}
\alias{findRange}
\title{Find values in sorted vectors}
\description{
  Find values in sorted vectors.
}
\usage{
findEqualGreater(x, value)
findEqualLess(x, value)
findEqualGreaterM(x, values)
findRange(x, values, NAOK = FALSE)
}
\arguments{
  \item{x}{numeric vector sorted in increasing order}
  \item{value}{value to find in \code{x}}
  \item{values}{numeric values to find in \code{x}}
  \item{NAOK}{don't check for NA values in \code{x}}
}
\details{
  \code{findEqualGreater} finds the index of the first value in
  \code{x} that is equal or greater than \code{value}. \code{findEqualLess}
  does same except that it finds equal or less. \code{findEqualGreaterM}
  creates an index of a vector by finding specified values.
  \code{findRange} locates the start and stop indicides of a range
  of two \code{x} values.

  The only things that save time at this point are
  \code{findeEqualGreaterM} (when the length of values approaches
  the lenght of x) and \code{findRange} (when \code{NAOK} is set
  to TRUE). They run in log(N) and N time, respectively.
}
\value{
  An integer vector with the position(s) of the values(s).
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{utilities}
\keyword{internal}
