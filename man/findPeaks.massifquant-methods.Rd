\name{findPeaks.massifquant-methods}
\docType{methods}
\alias{findPeaks.massifquant}
\alias{findPeaks.massifquant,xcmsRaw-method}
\title{Feature detection for high resolution LC/MS data}
\description{
  Kalman filter based feature detection for high resolution LC/MS data in centroid mode
  (currently experimental).
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  For Orbitrap Data with poor to acceptable chromatography, suggested default parameters.
  \code{
        findPeaks.massifquant(object, scanrange = c(1, length(object@scantime)), 
                 minIntensity = 6400, minCentroids = 12, 
                 consecMissedLim = 2, criticalVal = 1.7321, 
                 ppm = 10,  segs = 1, scanBack = 1)
    }
  
  For TOF Data with perfect chromatography, suggested default parameters.
    \code{
    findPeaks.massifquant(object, scanrange = c(1, length(object@scantime), 
                 minIntensity = 1800, minCentroids = 6, 
                 consecMissedLim = 1, criticalVal = 0.7111, 
                 ppm = 10,  segs = 1, scanBack = 1)
    }
}
}}
\details{ 
  This algorithm is most suitable for high resolution LC/\{OrbiTrap, TOF\}-MS data in centroid mode. Simultaneous kalman filters identify features and calculate their area under the curve. 
  Originally developed on LTQ Orbitrap data with much less than perfect chromatography, the default parameters are set to that specification. Users will find it useful to do some simple exploratory data analysis to find out where to set a minimum intensity, and identify how many scans an average feature may be. May we suggest using TOPPView as a visualization tool. Historicaly, the consecutiveMissedLim parameter should be set to (2) on Orbitrap data and (1) on TOF data, but never should exceed (4). The criticalVal parameter is perhaps most dificult to dial in appropriately and visual inspection of peak identification is the best suggested tool for quick optimization. The ppm, sets, and scanBack parameters have shown less influence than the other parameters and exist to give users flexibility and better accuracy.
  }
\arguments{
  \item{object}{\code{xcmsRaw} object}
  \item{scanrange}{scan range to process \code{scanrange = c(1, lastScan)} where lastScan is an integer}
   \item{minIntensity}{All real features should exceed this height.}
  \item{minCentroids}{A lower bound for how many scans a feature spans; a feature only incorporates one centroid per scan}
  \item{consecMissedLim}{As a feature is detected, the Kalman Filter may not find a centroid in every scan; After 1 or more misses, this consecutive missed limit informs massifquant when to stop a Kalman Filter to stop looking for a feature.}
  \item{criticalVal}{criticalVal helps determine the error bounds +/- of the Kalman Filter estimate. If the data has very fine mass resolution, a smaller critical val might be better and vice versa. A centroid apart of the feature should fall within these bounds on any given scan. Much like in the construction of a confidence interval, criticalVal loosely translates to be a  multiplier of the standard error estimate reported by the Kalman Filter. It is a relaxed application of the confidence interval because it doesn't change as more data is incorporated into the estimation proces, which would change the degrees of freedom and hence the critical value t.     }
  \item{ppm}{maximum m/z deviation in consecutive scans by ppm (parts per million).}
  \item{segs}{(segs = 1 #if turned on  segs = 0 #if turned off) With very few data points, sometimes a Kalman Filter "falls off" and stops tracking a feature prematurely. Another Kalman Filter is instantiated and begins following the rest of the signal. Because tracking is done backwards to forwards, this algorithmic defect leaves a real feature divided into two segments (segs for segmentation). With this option turned on, the program identifies segmented features and combines them into one with two sample t-test. The only danger is that samples with time consecutive features that appear conjoined to form a saddle will also be combined.}
  \item{scanBack}{(segs = 1 #if turned on  segs = 0 #if turned off) The convergence of a Kalman Filter to a feature's precise m/z mapping  is very fast, but sometimes it incorporates erroneous centroids as part of a feature (especially early on). The "scanBack" option removes the occasional outlier that lies beyond the converged bounds of the Kalman Filter. The option does not directly affect identification of a feature because it is a postprocessing measure; nonetheless, can potentially improve the quantitation by removing unlikely elements of an established feature.}
}
\value{
  A matrix with columns:

  \item{mz}{
    weighted mean (by intensity) of feature m/z across scans
  }
  \item{mzmin}{
    m/z peak minimum
  }
  \item{mzmax}{
    m/z peak maximum 
  }
  \item{rt}{
    retention time of peak midpoint estimate
  }
  \item{rtmin}{
    leading edge of peak retention time
  }
  \item{rtmax}{
    trailing edge of peak retention time
  }
  \item{into}{
    integrated peak intensity without any normalization
  }
  \item{maxo}{
    maximum peak intensity 
  }
}
\author{Chris Conley}
\encoding{UTF-8}
\references{
    yet another peak finder (still needing a title). High Impact Journal. Nov. 2011. 
}
\seealso{
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
