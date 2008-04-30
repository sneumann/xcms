\name{image.xcmsRaw}
\alias{image.xcmsRaw}
\alias{image.xcmsProfile}
\title{Plot log intensity image of a xcmsRaw/xcmsProfile object}
\description{
  Create log intensity false-color image of a xcmsRaw object plotted
  with m/z and retention time axes
}
\usage{
\method{image}{xcmsRaw}(x, col = rainbow(256), ...)
\method{image}{xcmsProfile}(x, col = rainbow(256), ...)
}
\arguments{
  \item{x}{xcmsRaw or xcmsProfile object}
  \item{col}{vector of colors to use for for the image}
  \item{\dots}{arguments for \code{profRange}}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{xcmsRaw-class}}, \code{\link{xcmsProfile-class}}
}
\keyword{hplot}
