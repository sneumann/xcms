\name{levelplot-methods}
\docType{methods}
\alias{levelplot,xcmsRaw-method}
\title{Plot log intensity image of a xcmsRaw object}
\description{
  Create an image of the raw (profile) data m/z against retention
  time, with the intensity color coded.
}
\section{Methods}{
\describe{
\item{x = "xcmsRaw"}{
  \code{
    image(x, log=TRUE, col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(256), ...)
  }
}
}}
\arguments{
  \item{x}{xcmsRaw object.}
  \item{log}{Whether the intensity should be log transformed.}
  \item{col.regions}{The color ramp that should be used for encoding of
    the intensity.}
  \item{\dots}{Arguments for \code{profRange}.}
}
\author{Johannes Rainer, \email{johannes.rainer@eurac.edu}}
\seealso{
  \code{\link{xcmsRaw-class}}, \code{\link{xcmsSet-class}}
}
\keyword{hplot}

