\name{plotSpec-methods}
\docType{methods}
%\alias{plotSpec-methods}
\alias{plotSpec}
\alias{plotSpec,xcmsRaw-method}
\title{Plot mass spectra from the profile matrix}
\description{
  Uses the pre-generated profile mode matrix to plot mass spectra
  over a specified retention time range.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotSpec(object, ident = FALSE, vline = numeric(0), ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{ident}{logical, use mouse to identify and label peaks}
  \item{vline}{numeric vector with locations of vertical lines}
  \item{...}{arguments passed to \code{\link{profRange}}}
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
