\name{descendZero}
\alias{descendZero}
\alias{descendMin}
\alias{descendValue}
\title{Find start and end points of a peak}
\description{
  Decends down the sides of a data peak and finds either the points
  greater than or equal to the zero intercept, the intercept with
  a given value, or the bottom of the first valley on each side.
}
\usage{
descendZero(y, istart = which.max(y))
descendValue(y, value, istart = which.max(y))
descendMin(y, istart = which.max(y))
}
\arguments{
  \item{y}{numeric vector with values}
  \item{istart}{starting point for descent}
  \item{value}{numeric value to descend to}
}
\value{
  An integer vector of length 2 with the starting and ending indicies
  of the peak start and end points.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{descendValue}}
}
\examples{
normdist <- dnorm(seq(-4, 4, .1)) - .1
xcms:::descendZero(normdist)
normdist[xcms:::descendZero(normdist)]
xcms:::descendValue(normdist, .15)
normdist[xcms:::descendValue(normdist, .15)]
xcms:::descendMin(normdist)
}
\keyword{iteration}
\keyword{internal}
