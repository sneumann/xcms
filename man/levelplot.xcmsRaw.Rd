\name{levelplot-methods}
\docType{methods}
\alias{levelplot,xcmsRaw-method}
\alias{levelplot,xcmsSet-method}
\title{Plot log intensity image of a xcmsRaw object}
\description{
  Create an image of the raw (profile) data m/z against retention
  time, with the intensity color coded.
}
\section{Methods}{
\describe{
\item{x = "xcmsRaw"}{
  \code{
    levelplot(x, log=TRUE, col.regions=colorRampPalette(brewer.pal(9,
    "YlOrRd"))(256), ...)
  }
}
\item{x = "xcmsSet"}{
  \code{
    levelplot(x, log=TRUE, col.regions=colorRampPalette(brewer.pal(9,
    "YlOrRd"))(256), rt="raw", ...)
  }
}
}}
\arguments{
  \item{x}{xcmsRaw object.}
  \item{log}{Whether the intensity should be log transformed.}
  \item{col.regions}{The color ramp that should be used for encoding of
    the intensity.}
  \item{rt}{wheter the original (\code{rt="raw"}) or the corrected
    (\code{rt="corrected"}) retention times should be used.}
  \item{\dots}{Arguments for \code{profRange}.}
}
\author{Johannes Rainer, \email{johannes.rainer@eurac.edu}}
\seealso{
  \code{\link{xcmsRaw-class}}, \code{\link{xcmsSet-class}}
}
\keyword{hplot}

