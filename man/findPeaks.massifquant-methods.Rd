\name{findPeaks.massifquant-methods}
\docType{methods}
\alias{findPeaks.massifquant}
\alias{findPeaks.massifquant,xcmsRaw-method}
\title{Feature detection for XC-MS data.}
\description{
Massifquant is a Kalman filter (KF) based feature detection for XC-MS data in centroid mode (currently in experimental stage). Optionally allows for calling the method "centWave" on features discovered by Massifquant to further refine the feature detection; to do so, supply any additional parameters specific to centWave (even more experimental). The method may be conveniently called through the xcmsSet(...) method.
}

\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
findPeaks.massifquant(object, ppm=10, peakwidth=c(20,50), snthresh=10,
     prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001,
     fitgauss=FALSE, scanrange= numeric(), noise=0,
     sleep=0, verbose.columns=FALSE, criticalValue = 1.125, consecMissedLimit = 2,
     unions = 1, checkBack = 0, withWave = 0)
}
}
}}

\details{
This algorithm's performance has been tested rigorously on high resolution LC/\{OrbiTrap, TOF\}-MS data in centroid mode. Simultaneous kalman filters identify features and calculate their area under the curve. The default parameters are set to operate on a complex LC-MS Orbitrap sample. Users will find it useful to do some simple exploratory data analysis to find out where to set a minimum intensity, and identify how many scans an average feature spans. The "consecMissedLimit" parameter has yielded good performance on Orbitrap data when set to (2) and on TOF data it was found best to be at (1). This may change as the algorithm has yet to be tested on many samples. The "criticalValue" parameter is perhaps most dificult to dial in appropriately and visual inspection of peak identification is the best suggested tool for quick optimization.  The "ppm" and "checkBack" parameters have shown less influence than the other parameters and exist to give users flexibility and better accuracy.
  }
\arguments{
The following arguments are specific to Massifquant. Any additional arguments supplied must correspond as specified by the method findPeaks.centWave.
\item{object}{An xcmsRaw object.}
\item{criticalValue}{Numeric: Suggested values: (0.1-3.0). This setting helps determine the the Kalman Filter prediciton margin of error. A real centroid belonging to a bonafide feature must fall within the KF prediction margin of error. Much like in the construction of a confidence interval, criticalVal loosely translates to be a  multiplier of the standard error of the prediction reported by the Kalman Filter. If the features in the XC-MS sample have a small mass deviance in ppm error, a smaller critical value might be better and vice versa.}
\item{consecMissedLimit}{Integer: Suggested values:(1,2,3). While a feature is in the proces of being detected by a Kalman Filter, the Kalman Filter may not find a predicted centroid in every scan. After 1 or more consecutive failed predictions, this setting informs Massifquant when to stop a Kalman Filter from following a candidate feature.}
\item{prefilter}{Numeric Vector: (Positive Integer, Positive Numeric): The first argument is only used if (withWave = 1); see centWave for details. The second argument specifies the minimum threshold for the maximum intensity of a feature that must be met.}
\item{peakwidth}{Integer Vector: (Positive Integer, Positive Integer): Only the first argument is used for Massifquant, which specifices the minimum feature length in time scans. If centWave is used, then the second argument is the maximum feature length  subject to being greater than the mininum feature length.}
\item{ppm}{The minimum estimated parts per million mass resolution a feature must possess.}
\item{unions}{Integer: set to 1 if apply t-test union on segmentation; set to 0 if no t-test to be applied on chromatographically continous features sharing same m/z range. Explanation: With very few data points, sometimes a Kalman Filter stops tracking a feature prematurely. Another Kalman Filter is instantiated and begins following the rest of the signal. Because tracking is done backwards to forwards, this algorithmic defect leaves a real feature divided into two segments or more. With this option turned on, the program identifies segmented features and combines them (merges them) into one with a two sample t-test. The potential danger of this option is that some truly distinct features may be merged.}
\item{withWave}{Integer: set to 1 if turned on; set to 0 if turned off. Allows the user to find features first with Massifquant and then filter those features with the second phase of centWave, which includes wavelet estimation.}
\item{checkBack}{Integer: set to 1 if turned on; set to 0 if turned off. The convergence of a Kalman Filter to a feature's precise m/z mapping  is very fast, but sometimes it incorporates erroneous centroids as part of a feature (especially early on). The "scanBack" option is an attempt to remove the occasional outlier that lies beyond the converged bounds of the Kalman Filter. The option does not directly affect identification of a feature because it is a postprocessing measure; it has not shown to be a extremely useful thus far and the default is set to being turned off.}
}
\value{
If the method findPeaks.massifquant(...) is used, then a matrix is returned with rows corresponding to features, and properties of the features listed with the following column names. Otherwise,
if centWave feature is used also (withWave = 1), or Massifquant is called through the xcmsSet(...) method, then their corresponding return values are used.
  \item{mz}{
    weighted m/z mean (weighted by intensity) of the feature
  }
  \item{mzmin}{
    m/z  lower boundary of the feature
    }
  \item{mzmax}{
    m/z  upper boundary of the feature
  }
  \item{rtmin}{
    starting scan time of the feature
  }
  \item{rtmax}{
     starting scan time of the feature
  }
  \item{into}{
   the raw quantitation (area under the curve) of the feature.
  }
  \item{area}{
    feature area that is not normalized by the scan rate.
  }
}
\examples{
library(faahKO)
library(xcms)
#load all the wild type and Knock out samples
cdfpath <- system.file("cdf", package = "faahKO")
## Subset to only the first 2 files.
cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)[1:2]

## Run the massifquant analysis. Setting the noise level to 10000 to speed up
## execution of the examples - in a real use case it should be set to a reasoable
## value.
xset <- xcmsSet(cdffiles, method = "massifquant",
                consecMissedLimit = 1,
                snthresh = 10,
                criticalValue = 1.73,
                ppm = 10,
                peakwidth= c(30, 60),
                prefilter= c(1,3000),
                noise = 10000,
                withWave = 0)
}
\author{Christopher Conley}
\encoding{UTF-8}
\references{
    Submitted for review. Christopher Conley, Ralf J .O Torgrip. Ryan Taylor, and John T. Prince.  "Massifquant: open-source Kalman filter based XC-MS feature detection". August 2013.
}
\seealso{
  \code{\link{centWave}} for the new user interface.
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsSet}}
  \code{\link{xcmsRaw}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
