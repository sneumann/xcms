\name{getEIC-methods}
\docType{methods}
%\alias{getEIC-methods}
\alias{getEIC}
%\alias{getEIC,xcmsRaw-method}
%\alias{getEIC,xcmsSet-method}
\title{Get extracted ion chromatograms for specified m/z ranges}
\description{
  Generate multiple extracted ion chromatograms for m/z values of
  interest.  For \code{xcmsSet} objects, reread original raw data
  and apply precomputed retention time correction, if applicable.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{getEIC(object, mzrange, step = 0.1)}
}
\item{object = "xcmsSet"}{
  \code{getEIC(object, mzrange, sampidx = seq(along = sampnames(object)))}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} or \code{xcmsSet} object}
  \item{mzrange}{
    either a two column matrix with minimum or maximum m/z or a
    matrix of any dimensions containing columns \code{mzmin} and
    \code{mzmax}
  }
  \item{step}{step size to use for profile generation}
  \item{sampidx}{
    integer vector with indicies of samples for which to get EICs
  }
}
\value{
  For \code{xcmsRaw} objects, a matrix with a column for each
  \code{mzmin}, \code{mzmax} pair. Rows correspond to individual
  scans.

  For \code{xcmsSet} objects, a list of matricies as described
  above, one for each sample.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{xcmsSet-class}}
}
\keyword{methods}
\keyword{file}
