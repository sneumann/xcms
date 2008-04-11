\name{plotScan-methods}
\docType{methods}
%\alias{plotScan-methods}
\alias{plotScan}
\alias{plotScan,xcmsRaw-method}
\title{Plot a single mass scan}
\description{
  Plot a single mass scan using the impulse representation. Most
  useful for centroided data.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotScan(object, scan, massrange = numeric(), ident = FALSE)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{scan}{
    integer with number of scan to plot
  }
  \item{massrange}{
    numeric vector of length >= 2 whose range will be used to select
    masses to plot
  }
  \item{ident}{
    logical, use mouse to interactively identify and label individual
    masses
  }
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{hplot}
