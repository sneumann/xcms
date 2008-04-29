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
    findPeaks.MSW(object, snthresh=3,
    scales=seq(1,22,3), nearbyPeak=TRUE,
    peakScaleRange=5, amp.Th=0.01,
    minNoiseLevel=amp.Th/SNR.Th,
    ridgeLength=24, tuneIn=FALSE,
    sleep=0, verbose.columns = FALSE)
    }
}
}}
\details{
  This is a wrapper around the peak picker in the bioconductor package
  MassSpecWavelet calling 'cwt', 'get.localMaximum.cwt',
     'get.ridge', 'identify.majorPeaks' and tuneIn.peakInfo.
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{snthresh}{signal to noise ratio cutoff}

  \item{scales}{scales of CWT}
  \item{nearbyPeak}{Determine whether to include the nearby small peaks of
          major peaks. TRUE by default}

  \item{sleep}{
    number of seconds to pause between plotting peak finding cycles
  }
  \item{verbose.columns}{
    additional peak meta data columns are returned 
  }
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
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
