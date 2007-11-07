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
    findPeaks.centWave(object, scanrange=c(1,length(object@scantime)), minEntries=4, dev=140e-6, snthresh=20, 
    minPeakWidth=7, noiserange=c(minPeakWidth*3,minPeakWidth*6), scales=c(5,7,9,12,16,20), maxGaussOverlap = 0.5,
    minPtsAboveBaseLine=4, scRangeTol=2, maxDescOutlier=floor(minPeakWidth/2), mzdiff=-0.001,
    rtdiff=-round(2/3 *minPeakWidth *mean(diff(object@scantime))),
    integrate=1, sleep=0, fitgauss = FALSE, verbose.columns = FALSE) 
  }
}
}}
\details{
  This algorithm is most suitable for high resolution centroid LC/TOF-MS data. In the first phase of the method 
  areas of high peak density (characterised by having at least \code{minEntries} peaks
  within \code{dev} ppm in consecutive scans) in the LC/MS map are located.
  In the second phase these \emph{regions of interest} (ROI) are further analysed. 
  Continuous wavelet transform (CWT) is used to locate chromatographic peaks on different scales. 
  
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{scanrange}{scan range to check for peaks}
  \item{minEntries}{minimum number of peaks(centroids) within the range of \code{dev} ppm}
  \item{dev}{see above}
  \item{snthresh}{signal to noise ratio cutoff}
  \item{noiserange}{
    neighbourhood size (number of scans) for each ROI used for noise/baseline analysis, two region sizes can be used
    for better local basline estimation
  }
  \item{minPeakWidth}{
    Minimum extent (number of scans) of a peak. Only used to estimate other parameters.
  }
  \item{scales}{scales to use for wavelet analysis, the scale range should correspond to the width range of the chromatographic peaks (number of scans) }
  \item{maxGaussOverlap}{maximal allowed overlapping (peak area)}
  \item{minPtsAboveBaseLine}{minimum number of data points above baseline}
  \item{scRangeTol}{tolerance (in scans) if gaussian peak center is not within the initial ROI}
  \item{maxDescOutlier}{Maximum number of outliers on the descent to find peak borders (used if integrate=2)}
  \item{mzdiff}{
    minimum difference in m/z for peaks with overlapping retention
    times, can be negative to allow overlap
  }
  \item{rtdiff}{
    minimum difference in retention time for peaks with overlapping m/z,
    can be negative to allow overlap
  }
  \item{integrate}{
    Integration method. If \code{=1} peak limits are found through descent on the mexican hat filtered data,
    if \code{=2} the descent is done on the real data. Method 2 is honest, while method 1 is more robust to noise.
  }
  \item{sleep}{
    number of seconds to pause between plotting peak finding cycles
  }
  \item{fitgauss}{
    logical, if TRUE a gaussian is fitted to each peak
  }
  
  \item{verbose.columns}{
    logical, if TRUE additional peak meta data columns are returned 
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
  \item{sn}{Signal/Noise ratio, calculated as (maxo / baseline value) }
  \item{egauss}{RMSE of Gaussian fit}
  \item{}{
      if \code{verbose.columns} is \code{TRUE} additionally :
  }
  \item{mu}{Gauss parameter $\mu$}
  \item{sigma}{Gauss parameter $\sigma$}
  \item{h}{Gauss parameter h}
  \item{f}{ROI number}
  \item{scale}{Scale on which the peak was localised}
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
