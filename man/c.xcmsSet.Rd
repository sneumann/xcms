\name{c.xcmsSet}
\alias{c.xcmsSet}
\title{Combine xcmsSet objects}
\description{
  Combines the samples and peaks from multiple \code{xcmsSet} objects
  into a single object. Group and retention time correction data
  are discarded. The \code{profinfo} list is set to be equal to the
  first object.
}
\usage{
c.xcmsSet(...)
}
\arguments{
  \item{\dots}{\code{xcmsSet} objects}
}
\value{
  A \code{xcmsSet} object.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{xcmsSet-class}}
}
\keyword{manip}
