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
\item{Orbitrap Data}{

  Please note that most of the parameter names are the same as findPeaks.centWave. Only additional parameters to massifquant are described in detail.

   \code{
    findPeaks.massifquant(object, ppm=5, peakwidth=c(16,700), snthresh=5,
    prefilter=c(15,1e5), mzCenterFun="wMean", integrate=1, mzdiff=-0.001, fitgauss=FALSE,
    scanrange= numeric(), noise=0, sleep=0, verbose.columns=FALSE,
    criticalValue = 1.7321, consecMissedLimit = 2, unions = 1, checkBack = 1, withWave = 1)
  }
}
}
}
\details{
  This algorithm is most suitable for high resolution LC/\{OrbiTrap, TOF\}-MS data in centroid mode. Simultaneous kalman filters identify features and calculate their area under the curve.
  Originally developed on LTQ Orbitrap data with much less than perfect chromatography, the default parameters are set to that specification. Users will find it useful to do some simple exploratory data analysis to find out where to set a minimum intensity, and identify how many scans an average feature may be. May we suggest using TOPPView as a visualization tool. Historicaly, the "consecMissedLimit" parameter should be set to (2) on Orbitrap data and (1) on TOF data, but never should exceed (4). The "criticalValue" parameter is perhaps most dificult to dial in appropriately and visual inspection of peak identification is the best suggested tool for quick optimization. The ppm and checkBack parameters have shown less influence than the other parameters and exist to give users flexibility and better accuracy.
  }
\arguments{
  \item{consecMissedLimit}{As a feature is detected, the Kalman Filter may not find a centroid in every scan; After 1 or more misses, this consecutive missed limit informs massifquant when to stop a Kalman Filter to stop looking for a feature.}
  \item{criticalValue}{criticalVal helps determine the error bounds +/- of the Kalman Filter estimate. If the data has very fine mass resolution, a smaller critical val might be better and vice versa. A centroid apart of the feature should fall within these bounds on any given scan. Much like in the construction of a confidence interval, criticalVal loosely translates to be a  multiplier of the standard error estimate reported by the Kalman Filter. It is a relaxed application of the confidence interval because it doesn't change as more data is incorporated into the estimation proces, which would change the degrees of freedom and hence the critical value t.     }
  \item{unions}{(unions = 1 #if turned on  unions = 0 #if turned off) With very few data points, sometimes a Kalman Filter "falls off" and stops tracking a feature prematurely. Another Kalman Filter is instantiated and begins following the rest of the signal. Because tracking is done backwards to forwards, this algorithmic defect leaves a real feature divided into two segments. With this option turned on, the program identifies segmented features and combines them (unites them) into one with two sample t-test. The only danger is that samples with time consecutive features that appear conjoined to form a saddle will also be combined.}
  \item{checkBack}{(checkBack = 1 #if turned on  checkBack = 0 #if turned off) The convergence of a Kalman Filter to a feature's precise m/z mapping  is very fast, but sometimes it incorporates erroneous centroids as part of a feature (especially early on). The "scanBack" option removes the occasional outlier that lies beyond the converged bounds of the Kalman Filter. The option does not directly affect identification of a feature because it is a postprocessing measure; nonetheless, can potentially improve the quantitation by removing unlikely elements of an established feature.}
  \item{withWave}{
    withWave = 1; # turns on centWave estimation
    withWave = 0; # otherwise  }
}
\value{
  In the case you choose to not use centWave for intensity estimation. It returns a matrix with columns as follows.

  \item{mz}{
    weighted mean (by intensity) of feature m/z across scans
  }
  \item{mzmin}{
    m/z peak minimum
  }
  \item{mzmax}{
    m/z peak maximum
  }
  \item{scmin}{
    starting scan number of feature
  }
  \item{scmax}{
     starting scan number of feature
  }
  \item{length}{
    the number of centroids in a feature
  }
  \item{area}{
    feature area that is not normalized by the scan rate.
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
