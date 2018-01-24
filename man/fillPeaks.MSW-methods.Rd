\name{fillPeaks.MSW-methods}
\docType{methods}
%\alias{fillPeaks.MSW-methods}
\alias{fillPeaks.MSW}
\alias{fillPeaks.MSW,xcmsSet-method}
\title{Integrate areas of missing peaks in FTICR-MS data}
\description{
  For each sample, identify peak groups where that sample is not
  represented. For each of those peak groups, integrate the signal
  in the region of that peak group and create a new peak.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{fillPeaks.MSW(object)}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
}
\details{
  After peak grouping, there will always be peak groups that do not
  include peaks from every sample. This method produces intensity
  values for those missing samples by integrating raw data in peak
  group region. In a given group, the start and ending m/z values for
  integration are defined by the median start and
  end points of the other detected peaks.

}
\note{
  In contrast to the \code{\link{fillPeaks.chrom}} method the maximum
  intensity reported in column \code{"maxo"} is not the maximum
  intensity measured in the expected peak area (defined by columns
  \code{"mzmin"} and \code{"mzmax"}), but the largest intensity of mz
  value(s) closest to the \code{"mzmed"} of the feature.
}
\value{
  A \code{xcmsSet} objects with filled in peak groups.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{getPeaks}}
  \code{\link{fillPeaks}}
}
\keyword{methods}
\keyword{file}
