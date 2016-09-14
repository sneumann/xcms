\name{findPeaks.addPredictedIsotopeFeatures-methods}
\docType{methods}
\alias{findPeaks.addPredictedIsotopeFeatures}
\alias{findPeaks.addPredictedIsotopeFeatures,xcmsRaw-method}
\title{Feature detection based on predicted isotope features for high resolution LC/MS data}
\description{
  Peak density and wavelet based feature detection aiming at isotope peaks for high resolution LC/MS data in centroid mode
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.centWave(object, ppm=25, peakwidth=c(20,50), 
    prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001, fitgauss=FALSE,
    scanrange= numeric(), noise=0, sleep=0, verbose.columns=FALSE, xcmsPeaks, snthresh=6.25, maxcharge=3, maxiso=5, mzIntervalExtension=TRUE)
  }
}
}}
\details{
  This algorithm is most suitable for high resolution LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase of the method isotope ROIs (regions of interest) in the LC/MS map are predicted.
  In the second phase these mass traces are further analysed.
  Continuous wavelet transform (CWT) is used to locate chromatographic peaks on different scales.
  The resulting peak list and the given peak list (\code{xcmsPeaks}) are merged and redundant peaks are removed.
}
\arguments{
  \item{object}{\code{xcmsSet} object}
  \item{ppm}{maxmial tolerated m/z deviation in consecutive scans, in ppm (parts per million)}
  \item{peakwidth}{Chromatographic peak width, given as range (min,max) in seconds}
  \item{prefilter}{\code{prefilter=c(k,I)}. Prefilter step for the first phase. Mass traces are only retained if they contain at least \code{k} peaks with intensity >= \code{I}.}
  \item{mzCenterFun}{Function to calculate the m/z center of the feature: \code{wMean} intensity weighted mean of the feature m/z values, \code{mean} mean of the feature m/z values, \code{apex} use m/z value at peak apex,
        \code{wMeanApex3} intensity weighted mean of the m/z value at peak apex and the m/z value left and right of it, \code{meanApex3} mean of the m/z value at peak apex and the m/z value left and right of it.
    }
  \item{integrate}{
    Integration method. If \code{=1} peak limits are found through descent on the mexican hat filtered data, if \code{=2} the descent is done on the real data. Method 2 is very accurate but prone to noise, while method 1 is more robust to noise but less exact.
   }
  \item{mzdiff}{
    minimum difference in m/z for peaks with overlapping retention
    times, can be negative to allow overlap
  }
  \item{fitgauss}{
    logical, if TRUE a Gaussian is fitted to each peak
  }
  \item{scanrange}{scan range to process}
  \item{noise}{optional argument which is useful for data that was centroided without any intensity threshold,
    centroids with intensity < \code{noise} are omitted from ROI detection}
  \item{sleep}{
    number of seconds to pause between plotting peak finding cycles
  }
  \item{verbose.columns}{
    logical, if TRUE additional peak meta data columns are returned
  }
  \item{xcmsPeaks}{peak list picked using the \code{centWave} algorithm with parameter \code{verbose.columns} set to TRUE (columns \code{scmin} and \code{scmax} needed)}
  \item{snthresh}{signal to noise ratio cutoff, definition see below.}
  \item{maxcharge}{max. number of the isotope charge.}
  \item{maxiso}{max. number of the isotope peaks to predict for each detected feature.}
  \item{mzIntervalExtension}{logical, if TRUE predicted isotope ROIs (regions of interest) are extended in the m/z dimension to increase the detection of low intensity and hence noisy peaks.}
}
\value{
  A matrix with columns:

  \item{mz}{
    weighted (by intensity) mean of peak m/z across scans
  }
  \item{mzmin}{
    m/z peak minimum
  }
  \item{mzmax}{
    m/z peak maximum
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
  \item{mu}{Gaussian parameter mu}
  \item{sigma}{Gaussian parameter sigma}
  \item{h}{Gaussian parameter h}
  \item{f}{Region number of m/z ROI where the peak was localised}
  \item{dppm}{m/z deviation of mass trace across scans in ppm}
  \item{scale}{Scale on which the peak was localised}
  \item{scpos}{Peak position found by wavelet analysis}
  \item{scmin}{Left peak limit found by wavelet analysis (scan number)}
  \item{scmax}{Right peak limit found by wavelet analysis (scan number)}
}
\author{Ralf Tautenhahn}
\encoding{UTF-8}
\references{
  Ralf Tautenhahn, Christoph B\"{o}ttcher, and Steffen Neumann
  "Highly sensitive feature detection for high resolution LC/MS"
  BMC Bioinformatics 2008, 9:504\\
  Hendrik Treutler and Steffen Neumann. 
  "Prediction, detection, and validation of isotope clusters in mass spectrometry data"
  Submitted to Metabolites 2016, Special Issue "Bioinformatics and Data Analysis"
}
\seealso{
  \code{\link{findPeaks.centWave}}
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
