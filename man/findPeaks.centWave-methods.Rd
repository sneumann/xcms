\name{findPeaks.centWave-methods}
\docType{methods}
\alias{findPeaks.centWave}
\alias{findPeaks.centWave,xcmsRaw-method}
\title{Feature detection for centroided LC/MS data}
\description{
  Peak density and wavelet based feature detection for centroided LC/MS data
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.centWave(object, scanrange=c(1,length(object@scantime)), minEntries=4, dev=140e-6, snthresh=10, noiserange=10, minPeakWidth=7, scales=c(5,7,9,12), maxGaussErr=0.25, maxGaussOverlap = 0.5, minPtsAboveBaseLine=3, rbwidth=minPeakWidth * 2, scRangeTol=2, maxDescOutlier=floor(minPeakWidth/2), mzdiff=-0.02, rtdiff=-round(2/3 *minPeakWidth *mean(diff(object@scantime))), sleep=0, verbose.columns = FALSE)
  }
}
}}
\details{
  This algorithm is useful only for centroid LC/MS data. In the first phase of the method 
  areas of high peak density (characterised by having at least \code{minEntries} peaks
  within \code{dev} ppm in consecutive scans) in the LC/MS map are located.
  In the second phase these \emph{regions of interest} (ROI) are further analysed. 
  Gaussian curves are fitted to the chromatographic peaks,
  and areas larger than \code{2*minPeakWidth} scans are checked
  whether they contain more than one peak (these may also overlap) using
  wavelet analysis. 
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{scanrange}{scan range to check for peaks}
  \item{minEntries}{minimum number of peaks(centroids) within the range of \code{dev} ppm}
  \item{dev}{see above}
  \item{snthresh}{signal to noise ratio cutoff}
  \item{noiserange}{
    neighbourhood size (number of scans) for each ROI used for noise/baseline analysis 
  }
  \item{minPeakWidth}{
    Minimum extent (number of scans) of a peak. Areas bigger than \code{2*minPeakWidth} 
    are checked whether they contain more than one peak.
  }
  \item{scales}{scales to use for wavelet analysis, the scale range should correspond to the width range of the chromatographic peaks (number of scans) }
  \item{maxGaussErr}{maximum root-mean-square error (RMSE) of Gauss fit, RMSE is close to zero if peak is gaussian-shaped}
  \item{maxGaussOverlap}{maximal allowed overlapping (peak area)}
  \item{minPtsAboveBaseLine}{used for baseline estimation}
  \item{rbwidth}{width of running median for baseline estimation}
  \item{scRangeTol}{tolerance (in scans) if gaussian peak center is not within the initial ROI}
  \item{maxDescOutlier}{Maximum number of outliers on the descent to find peak borders}
  \item{mzdiff}{
    minimum difference in m/z for peaks with overlapping retention
    times, can be negative to allow overlap
  }
  \item{rtdiff}{
    minimum difference in retention time for peaks with overlapping m/z,
    can be negative to allow overlap
  }
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
    m/z of minimum step
  }
  \item{mzmax}{
    m/z of maximum step
  }
  \item{rt}{
    retention time of peak midpoint
  }
  \item{rtmin}{
    leading edge of peak retention time
  }
  \item{rtmax}{
    trailing edge of peak retention time
  }
  \item{into}{
    integrated area of original (raw) peak
  }
  \item{maxo}{
    maximum intensity of original (raw) peak
  }
  \item{sn}{Signal/Noise ratio}
  \item{egauss}{RMSE of Gaussian fit}
  \item{}{
      if \code{verbose.columns} is \code{TRUE} additionally :
  }
  \item{mu}{Gauss parameter $\mu$}
  \item{sigma}{Gauss parameter $\sigma$}
  \item{h}{Gauss parameter h}
  \item{f}{ROI number}
  \item{scpos}{Peak position found by wavelet analysis}
  \item{scmin}{Left peak limit found by wavelet analysis (scan number)}
  \item{scmax}{Right peak limit found by wavelet analysis (scan number)}
  \item{lmin}{leading edge of peak (scan number)}
  \item{lmax}{trailing edge of peak (scan number)}

}
\author{Ralf Tautenhahn, \email{rtautenh@ipb-halle.de}}
\seealso{
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
