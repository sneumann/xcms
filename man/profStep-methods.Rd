\name{profStep-methods}
\docType{methods}
%\alias{profStep-methods}
%\alias{profStep<--methods}
\alias{profStep}
\alias{profStep<-}
\alias{profStep,xcmsRaw-method}
\alias{profStep<-,xcmsRaw-method}
\title{Get and set m/z step for generating profile data}
\description{
  These methods get and set the m/z step for generating profile
  (matrix) data from raw mass spectral data. Smaller steps yield
  more precision at the cost of greater memory usage.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{profStep(object)}
}
}}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{profMethod}}
}
\keyword{methods}
\keyword{manip}
