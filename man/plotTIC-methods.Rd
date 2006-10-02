\name{plotTIC-methods}
\docType{methods}
%\alias{plotTIC-methods}
\alias{plotTIC}
\alias{plotTIC,xcmsRaw-method}
\title{Plot total ion count}
\description{
  Plot chromatogram of total ion count. Optionally allow identification
  of target peaks and viewing/identification of individual spectra.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotTIC(object, ident = FALSE, msident = FALSE)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{ident}{
    logical, use mouse to identify and label chromatographic peaks
  }
  \item{ident}{
    logical, use mouse to identify and label spectral peaks
  }
}
\value{
  If \code{ident == TRUE}, an integer vector with the indecies of
  the points that were identified. Otherwise a two-column matrix
  with the plotted points.
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{hplot}
\keyword{iplot}
