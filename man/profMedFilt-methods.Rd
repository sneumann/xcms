\name{profMedFilt-methods}
\docType{methods}
%\alias{profMedFilt-methods}
\alias{profMedFilt}
\alias{profMedFilt,xcmsRaw-method}
\title{Median filtering of the profile matrix}
\description{
  Apply a median filter of given size to a profile matrix.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{profMedFilt(object, massrad = 0, scanrad = 0)}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} object}
  \item{massrad}{
    number of m/z grid points on either side to use for median
    calculation
  }
  \item{scanrad}{
    number of scan grid points on either side to use for median
    calculation
  }
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{medianFilter}}
}
\keyword{methods}
\keyword{manip}
