\name{descendZero}
\alias{descendZero}
\alias{descendMin}
\title{Find start and end points of a peak}
\description{
  Decends down the sides of a data peak and finds either the zero
  intercept or the first valley on each side.
}
\usage{
descendZero(y, istart = max(y))
descendMin(y, istart = max(y))
}
\arguments{
  \item{y}{numeric vector with values}
  \item{istart}{starting point for descent}
}
\value{
  A vector of length 2 with the starting and ending indicies of the
  peak start and end points.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\keyword{manip}
\keyword{internal}
