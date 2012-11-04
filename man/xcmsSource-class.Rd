\name{xcmsSource-class}
\docType{class}
\alias{xcmsSource-class}
\title{Virtual class for raw data sources}
\description{
  This virtual class provides an implementation-independent way to load
  mass spectrometer data from various sources for use in an
  \code{\linkS4class{xcmsRaw}} object. Subclasses can be defined to
  enable data to be loaded from user-specified sources. The virtual
  class \code{\linkS4class{xcmsFileSource}} is included out of the box
  which contains a file name as a character string.

  When implementing child classes of \code{xcmsSource}, a corresponding
  \code{\link{loadRaw-methods}} method must be provided which accepts
  the \code{xcmsSource} child class and returns a list in the format
  described in \code{\link{loadRaw-methods}}.
}
\section{Objects from the Class}{
  A virtual Class: No objects may be created from it.
}

\author{Daniel Hackney, \email{dan@haxney.org}}
\seealso{
  \code{\link{xcmsSource-methods}} for creating \code{xcmsSource}
  objects in various ways.
}
\keyword{classes}
