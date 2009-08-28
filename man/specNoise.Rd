\name{specNoise}
\alias{specNoise}
\title{Calculate noise for a sparse continuum mass spectrum}
\description{
  Given a sparse continuum mass spectrum, determine regions where
  no signal is present, substituting half of the minimum intensity
  for those regions. Calculate the noise level as the weighted mean
  of the regions with signal and the regions without signal.
  If there is only one raw peak, return zero.
}
\usage{
specNoise(spec, gap = quantile(diff(spec[, "mz"]), 0.9))
}
\arguments{
  \item{spec}{matrix with named columns \code{mz} and \code{intensity}}
  \item{gap}{
    threshold above which to data points are considerd to be separated
    by a blank region and not bridged by an interpolating line
  }
}
\details{
  The default gap value is determined from the 90th percentile of
  the pair-wise differences between adjacent mass values.
}
\value{
  A numeric noise level
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{getSpec}},
  \code{\link{specPeaks}}
}
\keyword{manip}
