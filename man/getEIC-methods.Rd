\name{getEIC-methods}
\docType{methods}
%\alias{getEIC-methods}
\alias{getEIC}
\alias{getEIC,xcmsRaw-method}
\alias{getEIC,xcmsSet-method}
\title{Get extracted ion chromatograms for specified m/z ranges}
\description{
  Generate multiple extracted ion chromatograms for m/z values of
  interest.  For \code{xcmsSet} objects, reread original raw data
  and apply precomputed retention time correction, if applicable.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{getEIC(object, mzrange, rtrange = NULL, step = 0.1)}
}
\item{object = "xcmsSet"}{
  \code{getEIC(object, mzrange, rtrange = 200, groupidx, 
               sampleidx = sampnames(object), rt = c("corrected", "raw"))}
}
}}
\arguments{
  \item{object}{the \code{xcmsRaw} or \code{xcmsSet} object}
  \item{mzrange}{
    either a two column matrix with minimum or maximum m/z or a
    matrix of any dimensions containing columns \code{mzmin} and
    \code{mzmax}

    for \code{xcmsSet} objects, if left blank the group data will
    be used instead
  }
  \item{rtrange}{
    a two column matrix the same size as \code{mzrange} with minimum
    and maximum retention times between which to return EIC data
    points
    
    for \code{xcmsSet} objects, it may also be a single number
    specifying the time window around the peak to return EIC data
    points
  }
  \item{step}{step size to use for profile generation}
  \item{groupidx}{
    either character vector with names or integer vector with
    indicies of peak groups for which to get EICs
  }
  \item{sampleidx}{
    either character vector with names or integer vector with
    indicies of samples for which to get EICs
  }
  \item{rt}{
    \code{"corrected"} for using corrected retention times, or
    \code{"raw"} for using raw retention times
  }
}
\value{
  For \code{xcmsRaw} objects, if \code{rtrange} is NULL, a matrix
  with a column for each \code{mzmin}, \code{mzmax} pair. Rows
  correspond to individual scans. Otherwise, a list of two column
  (retention time/intensity) matricies, one for each \code{mzmin},
  \code{mzmax} pair.

  For \code{xcmsSet} objects, an \code{xcmsEIC} object.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{xcmsSet-class}},
  \code{\link{xcmsEIC-class}}
}
\keyword{methods}
\keyword{file}
