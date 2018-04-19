\name{fillPeaks.chrom-methods}
\docType{methods}
%\alias{fillPeaks.chrom-methods}
\alias{fillPeaks.chrom}
\alias{fillPeaks.chrom,xcmsSet-method}
\title{Integrate areas of missing peaks}
\description{
  For each sample, identify peak groups where that sample is not
  represented. For each of those peak groups, integrate the signal
  in the region of that peak group and create a new peak.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{fillPeaks.chrom(object, nSlaves=0,expand.mz=1,expand.rt=1,
    BPPARAM = bpparam())}
}
}}

\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{nSlaves}{(DEPRECATED): number of slaves/cores to be used for
    parallel peak filling.
    MPI is used if installed, otherwise the snow package is employed for
    multicore support. If none of the two packages is available it uses
    the parallel package for parallel processing on multiple CPUs of the
    current machine. Users are advised to use the \code{BPPARAM}
    parameter instead.}
  \item{expand.mz}{Expansion factor for the m/z range used for integration.}
  \item{expand.rt}{Expansion factor for the rentention time range used
    for integration.}
  \item{BPPARAM}{allows to define a specific parallel processing setup
  for the current task (see \code{\link{bpparam}} from the
  \code{BiocParallel} package help more information). The default uses
  the globally defined parallel setup.}
}
\details{
  After peak grouping, there will always be peak groups that do not
  include peaks from every sample. This method produces intensity
  values for those missing samples by integrating raw data in peak
  group region. In a given group, the start and ending retention
  time points for integration are defined by the median start and
  end points of the other detected peaks. The start and end m/z
  values are similarly determined. Intensities
  can be still be zero, which is a rather unusual intensity for a peak.
  This is the case if e.g. the raw data was threshholded, and the
  integration area contains no actual raw intensities, or if one sample
  is miscalibrated, such thet the raw data points are (just) outside
  the integration area.

  Importantly, if retention time correction data is available, the
  alignment information is used to more precisely integrate the
  propper region of the raw data. If the corrected retention time
  is beyond the end of the raw data, the value will be not-a-number (NaN).

}
\value{
  A \code{xcmsSet} objects with filled in peak groups (into and maxo).
}

\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{getPeaks}}
  \code{\link{fillPeaks}}
}
\keyword{methods}
\keyword{file}
