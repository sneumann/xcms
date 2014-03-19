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
  \code{fillPeaks(object, method="")}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{method}{the filling method}
}
\details{
  After peak grouping, there will always be peak groups that do not
  include peaks from every sample. This method produces intensity
  values for those missing samples by integrating raw data in peak
  group region. According to the type of raw-data there are 2 different
  methods available. for filling gcms/lcms data the method "chrom" integrates
  raw-data in the chromatographic domain, whereas "MSW" is used for
  peaklists without retention-time information like those from
  direct-infusion spectra.
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
