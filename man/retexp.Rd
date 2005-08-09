\name{retexp}
\alias{retexp}
\title{Set retention time window to a specified width}
\description{
  Expands (or contracts) the retention time window in each row of
  a matrix as defined by the \code{retmin} and \code{retmax} columns.
}
\usage{
retexp(peakrange, width = 200)
}
\arguments{
  \item{peakrange}{maxtrix with columns \code{retmin} and \code{retmax}}
  \item{width}{new width for the window}
}
\value{
  The altered matrix.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{getEIC}}
}
\keyword{manip}
