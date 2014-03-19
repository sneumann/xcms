\name{profMethod-methods}
\docType{methods}
%\alias{profMethod-methods}
%\alias{profMethod<--methods}
\alias{profMethod}
\alias{profMethod<-}
\alias{profMethod,xcmsRaw-method}
\alias{profMethod<-,xcmsRaw-method}
\title{Get and set method for generating profile data}
\description{
  These methods get and set the method for generating profile
  (matrix) data from raw mass spectral data. It can currently be
  \code{bin}, \code{binlin}, \code{binlinbase}, or \code{intlin}.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{profMethod(object)}
}
}}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profMethod}},
  \code{\link{profBin}},
  \code{\link{plotSpec}},
  \code{\link{plotChrom}},
  \code{\link{findPeaks}}
}
\keyword{methods}
\keyword{manip}
