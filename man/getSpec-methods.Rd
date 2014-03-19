\name{getSpec-methods}
\docType{methods}
%\alias{getSpec-methods}
\alias{getSpec}
\alias{getSpec,xcmsRaw-method}
\title{Get average m/z and intensity values for multiple mass scans}
\description{
  Return full-resolution averaged data from multiple mass scans.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{getSpec(object, ...)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{...}{
    arguments passed to \code{\link{profRange}} used to sepecify the
    spectral segments of interest for averaging
  }
}
\details{
  Based on the mass points from the spectra selected, a master
  unique list of masses is generated. Every spectra is interpolated
  at those masses and then averaged.
}
\value{
  A matrix with two columns:
  \item{mz}{m/z values}
  \item{intensity}{intensity values}
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profRange}},
  \code{\link{getScan}}
}
\keyword{methods}
\keyword{manip}
