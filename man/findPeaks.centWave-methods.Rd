\name{findPeaks.centWave-methods}
\docType{methods}
\alias{findPeaks.centWave}
\alias{findPeaks.centWave,xcmsRaw-method}
\title{Feature detection for high resolution LC/MS data}
\description{
  Peak density and wavelet based feature detection for high resolution LC/MS data in centroid mode
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.centWave(object, ppm=25, peakwidth=c(20,50), snthresh=20,
    prefilter=c(3,100), integrate=1, mzdiff=-0.001, fitgauss=FALSE,
    scanrange= c(1, length(object@scantime)), sleep=0, verbose.columns=FALSE)
  }
}
}}
\details{
  This algorithm is most suitable for high resolution LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase of the method mass traces (characterised as regions with less than \code{ppm} m/z deviation in consecutive scans) in the LC/MS map are located.
  In the second phase these mass traces are further analysed. 
  Continuous wavelet transform (CWT) is used to locate chromatographic peaks on different scales. 
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{ppm}{maxmial tolerated m/z deviation in consecutive scans, in ppm (parts per million)}
  \item{peakwidth}{Chromatographic peak width, given as range (min,max) in seconds}
  \item{snthresh}{signal to noise ratio cutoff, definition see below.}
  \item{prefilter}{\code{prefilter=c(k,I)}. Prefilter step for the first phase. Mass traces are only retained if they contain at least \code{k} peaks with intensity >= \code{I}.}
  
  \item{integrate}{
    Integration method. If \code{=1} peak limits are found through descent on the mexican hat filtered data, if \code{=2} the descent is done on the real data. Method 2 is very accurate but prone to noise, while method 1 is more robust to noise but less exact.
   }
  \item{mzdiff}{
    minimum difference in m/z for peaks with overlapping retention
    times, can be negative to allow overlap
  }
  \item{fitgauss}{
    logical, if TRUE a gaussian is fitted to each peak
  }
  \item{scanrange}{scan range to process}
  \item{sleep}{
    number of seconds to pause between plotting peak finding cycles
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
    integrated peak intensity
  }
   \item{intb}{
    baseline corrected integrated peak intensity
  }
  \item{maxo}{
    maximum peak intensity 
  }
  \item{sn}{
    Signal/Noise ratio, defined as \code{(maxo - baseline)/sd}, where \cr
    \code{maxo} is the maximum peak intensity, \cr
    \code{baseline} the estimated baseline value and \cr
    \code{sd} the standard deviation of local chromatographic noise.
  }
  \item{egauss}{RMSE of Gaussian fit}
  \item{}{
      if \code{verbose.columns} is \code{TRUE} additionally :
  }
  \item{mu}{Gauss parameter mu}
  \item{sigma}{Gauss parameter sigma}
  \item{h}{Gauss parameter h}
  \item{f}{number of m/z ROI where the peak was localised}
  \item{dppm}{m/z deviation of mass trace in ppm}
  \item{scale}{Scale on which the peak was localised}
  \item{scpos}{Peak position found by wavelet analysis}
  \item{scmin}{Left peak limit found by wavelet analysis (scan number)}
  \item{scmax}{Right peak limit found by wavelet analysis (scan number)}
}
\author{Ralf Tautenhahn, \email{rtautenh@ipb-halle.de}}
\seealso{
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
