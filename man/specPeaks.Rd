\name{specPeaks}
\alias{specPeaks}
\title{Identify peaks in a sparse continuum mode spectrum}
\description{
  Given a spectrum, identify and list significant peaks as determined
  by several criteria.
}
\usage{
specPeaks(spec, sn = 20, mzgap = 0.2)
}
\arguments{
  \item{spec}{matrix with named columns \code{mz} and \code{intensity}}
  \item{sn}{minimum signal to noise ratio}
  \item{mzgap}{
    minimal distance between adjacent peaks, with smaller peaks
    being excluded
  }
}
\details{
  Peaks must meet two criteria to be considered peaks: 1) Their
  s/n ratio must exceed a certain threshold. 2) They must not be
  within a given distance of any greater intensity peaks.
}
\value{
  A matrix with columns:

  \item{mz}{m/z at maximum peak intensity}
  \item{intensity}{maximum intensity of the peak}
  \item{fwhm}{full width at half max of the peak}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{getSpec}},
  \code{\link{specNoise}}
}
\keyword{manip}
