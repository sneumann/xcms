\name{fillPeaks-methods}
\docType{methods}
%\alias{fillPeaks-methods}
\alias{fillPeaks}
\alias{fillPeaks,xcmsSet-method}
\title{Integrate areas of missing peaks}
\description{
  For each sample, identify peak groups where that sample is not
  represented. For each of those peak groups, integrate the signal
  in the region of that peak group and create a new peak.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{fillPeaks(object)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
}
\details{
  After peak grouping, there will always be peak groups that do not
  include peaks from every sample. This method produces intensity
  values for those missing samples by integrating raw data in peak
  group region. In a given group, the start and ending retention
  time points for integration are defined by the median start and
  end points of the other detected peaks. The start and end m/z
  values are similarly determined.

  Importantly, if retention time correction data is available, the
  alignment information is used to more precisely integrate the
  propper region of the raw data.
}
\value{
  A \code{xcmsSet} objects with filled in peak groups.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{getPeaks}}
}
\keyword{methods}
\keyword{file}
