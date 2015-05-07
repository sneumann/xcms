\name{getXcmsRaw-methods}
\docType{methods}
%\alias{getEIC-methods}
\alias{getXcmsRaw}
\alias{getXcmsRaw,xcmsSet-method}
\title{Load the raw data for one or more files in the xcmsSet}
\description{
  Reads the raw data applies evential retention time corrections and
  waters Lock mass correction and
  returns it as an \code{xcmsRaw} object (or list of \code{xcmsRaw}
  objects) for one or more files of the \code{xcmsSet} object.
}
\section{Methods}{
\describe{
\item{object = "xcmsSet"}{
  \code{getXcmsRaw(object, sampleidx=1,
    profmethod=profinfo(object)$method, profstep=profinfo(object)$step,
    rt=c("corrected", "raw"), ...
  )}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{sampleidx}{
    The index of the sample for which the raw data should be
    returned. Can be a single number or a numeric vector with the
    indices. Alternatively, the file name can be specified.
  }
  \item{profmethod}{
    The profile method.
  }
  \item{profstep}{
    The profile step.
  }
  \item{rt}{
    Whether corrected or raw retention times should be returned.
  }
  \item{...}{
    Additional arguments submitted to the \code{\link{xcmsRaw}} function.
  }
}
\value{
  A single \code{xcmsRaw} object or a list of \code{xcmsRaw} objects.
}
\author{Johannes Rainer, \email{johannes.rainer@eurac.edu}}
\seealso{
  \code{\link{xcmsRaw-class}},
}
\keyword{methods}
\keyword{file}
