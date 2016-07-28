\name{getEIC-methods}
\docType{methods}
%\alias{getEIC-methods}
\alias{getEIC}
\alias{getEIC,xcmsRaw-method}
\alias{getEIC,xcmsSet-method}
\title{Get extracted ion chromatograms for specified m/z ranges}
\description{
  Generate multiple extracted ion chromatograms for m/z values of
  interest. For \code{xcmsSet} objects, reread original raw data
  and apply precomputed retention time correction, if applicable.

  Note that this method will \emph{always} return profile, not raw
  data (with profile data being the binned data along M/Z). See details
  for further information.
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
    Either a two column matrix with minimum or maximum m/z or a
    matrix of any dimensions containing columns \code{mzmin} and
    \code{mzmax}. If not specified, the method for \code{xcmsRaw} returns
    the base peak chromatogram (BPC, i.e. the most intense signal for
    each RT across all m/z).

    For \code{xcmsSet} objects the group data will be used if
    \code{mzrange} is not provided.
  }
  \item{rtrange}{
    A two column matrix the same size as \code{mzrange} with minimum
    and maximum retention times between which to return EIC data
    points. If not specified, the method returns the chromatogram for
    the full RT range.

    For \code{xcmsSet} objects, it may also be a single number
    specifying the time window around the peak to return EIC data
    points
  }
  \item{step}{step (bin) size to use for profile generation. Note that a
    value of \code{step = 0} is not supported.}
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
\details{
  In contrast to the \code{\link{rawEIC}} method, that extracts the
  actual raw values, this method extracts them from the object's profile
  matrix (or if the provided \code{step} argument does not match the
  \code{profStep} of the \code{object} the profile matrix is calculated
  on the fly and the values returned).
}
\value{
  For \code{xcmsSet} and \code{xcmsRaw} objects, an \code{xcmsEIC} object.
}
\seealso{
  \code{\link{xcmsRaw-class}},
  \code{\link{xcmsSet-class}},
  \code{\link{xcmsEIC-class}},
  \code{\link{rawEIC}}
}
\keyword{methods}
\keyword{file}
