\name{plotChrom-methods}
\docType{methods}
%\alias{plotChrom-methods}
\alias{plotChrom}
\alias{plotChrom,xcmsRaw-method}
\title{Plot extracted ion chromatograms from the profile matrix}
\description{
  Uses the pre-generated profile mode matrix to plot averaged or
  base peak extracted ion chromatograms over a specified mass range.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{plotChrom(object, base = FALSE, ident = FALSE,
                  fitgauss = FALSE, vline = numeric(0), ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{base}{logical, plot a base-peak chromatogram}
  \item{ident}{logical, use mouse to identify and label peaks}
  \item{fitgauss}{logical, fit a gaussian to the largest peak}
  \item{vline}{numeric vector with locations of vertical lines}
  \item{...}{arguments passed to \code{\link{profRange}}}
}
\value{
  If \code{ident == TRUE}, an integer vector with the indecies of
  the points that were identified. If \code{fitgauss == TRUE}, a
  \code{nls} model with the fitted gaussian. Otherwise a two-column
  matrix with the plotted points.
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{hplot}
\keyword{iplot}
