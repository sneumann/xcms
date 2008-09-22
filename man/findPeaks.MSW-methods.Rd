\name{findPeaks.MSW-methods}
\docType{methods}
\alias{findPeaks.MSW}
\alias{findPeaks.MSW,xcmsRaw-method}
\title{Feature detection for single-spectrum non-chromatography MS data}
\description{
  Processing Mass Spectrometry direct-injection
  spectrum by using wavelet based algorithm.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.MSW(object, snthresh=3, verbose.columns = FALSE, ...)
    }
}
}}
\details{
  This is a wrapper around the peak picker in the bioconductor package
  MassSpecWavelet calling \code{peakDetectionCWT} and \code{tuneInPeakInfo}.
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{snthresh}{signal to noise ratio cutoff}
  \item{verbose.columns}{
    additional peak meta data columns are returned 
  }
  \item{...}{Parameters passed to \code{peakDetectionCWT}}
}
\value{
  A matrix with columns:

  \item{mz}{
    weighted (by intensity) mean of peak m/z across scans
  }
  \item{mzmin}{
    same as mz
  }
  \item{mzmax}{
    same as mz
  }
  \item{rt}{
    always -1
  }
  \item{rtmin}{
    always -1
  }
  \item{rtmax}{
    always -1
  }
  \item{into}{
    integrated area of original (raw) peak
  }
  \item{maxo}{
    intensity of original (raw) peak at the centroid position
  }
  \item{sn}{Signal/Noise ratio}
}
\author{Steffen Neumann, \email{sneumann@ipb-halle.de}}
\seealso{
  \code{\link[MassSpecWavelet]{peakDetectionCWT}}
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
