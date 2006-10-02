\name{profMedFilt-methods}
\docType{methods}
%\alias{profMedFilt-methods}
\alias{profMedFilt}
\alias{profMedFilt,xcmsRaw-method}
\title{Apply a median filter to the profile matrix}
\description{
  
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
