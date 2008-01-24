\name{getPeaks-methods}
\docType{methods}
%\alias{getPeaks-methods}
\alias{getPeaks}
\alias{getPeaks,xcmsRaw-method}
\title{Get peak intensities for specified regions}
\description{
  Integrate extracted ion chromatograms in pre-defined defined
  regions. Return output similar to \code{\link{findPeaks}}.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{getPeaks(object, peakrange, step = 0.1, pipeline = new("xcmsPipelineProfile"))}
}
}}
\arguments{
  \item{object}{the \code{xcmsSet} object}
  \item{peakrange}{
    matrix or data frame with 4 columns: \code{mzmin}, \code{mzmax},
    \code{rtmin}, \code{rtmax} (they must be in that order or named)
  }
  \item{step}{step size to use for profile generation}
  \item{pipeline}{an \code{\linkS4class{xcmsPipeline}} for producing a
    profile matrix. In the canonical workflow, this would usually come from an
    \code{\linkS4class{xcmsProfile}} object. The \code{step} parameter
    above will override the corresponding parameter in the pipeline. If
    \code{pipeline} is omitted, the \code{object} is searched for a
    profile generation pipeline. If not found, a
    default pipeline is constructed.}
}
\value{
  A matrix with columns:
  \item{i}{
    rank of peak identified in merged EIC (<= \code{max}), always \code{NA}
  }
  \item{mz}{
    weighted (by intensity) mean of peak m/z across scans
  }
  \item{mzmin}{
    m/z of minimum step
  }
  \item{mzmax}{
    m/z of maximum step
  }
  \item{ret}{
    retention time of peak midpoint
  }
  \item{retmin}{
    leading edge of peak retention time
  }
  \item{retmax}{
    trailing edge of peak retention time
  }
  \item{into}{
    integrated area of original (raw) peak
  }
  \item{intf}{
    integrated area of filtered peak, always \code{NA}
  }
  \item{maxo}{
    maximum intensity of original (raw) peak
  }
  \item{maxf}{
    maximum intensity of filtered peak, always \code{NA}
  }
}
\seealso{
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
\keyword{manip}
