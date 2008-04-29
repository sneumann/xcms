\name{image-methods}
\docType{methods}
\alias{image.xcmsRaw}
\alias{image, image-methods}
\title{Plot log intensity image of a xcmsRaw object}
\description{
  Create log intensity false-color image of a xcmsRaw object plotted
  with m/z and retention time axes
}
\section{Methods}{
\describe{
\item{x = "xcmsRaw"}{
  \code{
    image(x, col = rainbow(256), ...)
  }
}
}}
\arguments{
  \item{x}{xcmsRaw object}
  \item{col}{vector of colors to use for for the image}
  \item{\dots}{arguments for \code{profRange}}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{hplot}
