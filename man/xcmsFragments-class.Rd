\name{xcmsFragments-class}
\docType{class}
\alias{xcmsFragments-class}
\alias{collect}


\title{Class xcmsFragments, a class for handling Tandem MS and MS$^n$ data}
\description{  
  This class is similar to \code{\link{xcmsSet}} because it stores peaks
  from a number of individual files. However,  xcmsFragments keeps
  Tandem MS and e.g. Ion Trap or Orbitrap MS$^n$ peaks, including the
  parent ion relationships.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{xcmsFragments}}
  constructor and filled with peaks using the collect method.  
}
\section{Slots}{
  \describe{
    \item{\code{peaks}:}{
      matrix with colmns peakID (MS1 parent in corresponding xcmsSet),
      MSnParentPeakID (parent peak within this xcmsFragments), msLevel
      (e.g. 2 for Tandem MS), rt (retention time in case of LC data), mz
      (fragment mass-to-charge), intensity (peak intensity as reported
      by XXX) 
    }
  }
}
\section{Methods}{
  \describe{
    \item{\link{findPeaks}}{
      \code{signature(object = "xcmsFragments")}: feature detection using
      matched filtration in the chromatographic time domain
    }
    \item{\link{getEIC}}{
      \code{signature(object = "xcmsFragments")}: get extracted ion
      chromatograms in specified m/z ranges
    }
    \item{\link{getPeaks}}{
      \code{signature(object = "xcmsRaw")}: get data for peaks in
      specified m/z and time ranges
    }
    \item{\link{getScan}}{
      \code{signature(object = "xcmsRaw")}: get m/z and intensity
      values for a single mass scan
    }
    \item{\link{getSpec}}{
      \code{signature(object = "xcmsRaw")}: get average m/z and
      intensity values for multiple mass scans
    }
    \item{\link[xcms:image.xcmsRaw]{image}}{
      \code{signature(x = "xcmsRaw")}: get data for peaks in
      specified m/z and time ranges
    }
    \item{\link{plotChrom}}{
      \code{signature(object = "xcmsRaw")}: plot a chromatogram
      from profile data
    }
    \item{\link{plotRaw}}{
      \code{signature(object = "xcmsRaw")}: plot locations of raw
      intensity data points
    }
    \item{\link{plotScan}}{
      \code{signature(object = "xcmsRaw")}: plot a mass spectrum
      of an individual scan from the raw data
    }
    \item{\link{plotSpec}}{
      \code{signature(object = "xcmsRaw")}: plot a mass spectrum
      from profile data
    }
    \item{\link{plotSurf}}{
      \code{signature(object = "xcmsRaw")}: experimental method for
      plotting 3D surface of profile data with \code{rgl}.
    }
    \item{\link{plotTIC}}{
      \code{signature(object = "xcmsRaw")}: plot total ion count
      chromatogram
    }
    \item{\link{profMedFilt}}{
      \code{signature(object = "xcmsRaw")}: median filter profile
      data in time and m/z dimensions
    }
    \item{\link{profMethod<-}}{
      \code{signature(object = "xcmsRaw")}: change the method of
      generating the \code{profile} matrix
    }
    \item{\link{profMethod}}{
      \code{signature(object = "xcmsRaw")}: get the method of
      generating the \code{profile} matrix
    }
    \item{profMz}{
      \code{signature(object = "xcmsRaw")}: get vector of m/z values
      for each row of the \code{profile} matrix
    }
    \item{\link{profRange}}{
      \code{signature(object = "xcmsRaw")}: interpret flexible ways
      of specifying subsets of the \code{profile} matrix
    }
    \item{\link{profStep<-}}{
      \code{signature(object = "xcmsRaw")}: change the m/z step
      used for generating the \code{profile} matrix
    }
    \item{\link{profStep}}{
      \code{signature(object = "xcmsRaw")}: get the m/z step used
      for generating the \code{profile} matrix
    }
    \item{revMz}{
      \code{signature(object = "xcmsRaw")}: reverse the order of the
      data points for each scan
    }
    \item{sortMz}{
      \code{signature(object = "xcmsRaw")}: sort the data points
      by increasing m/z for each scan
    }
  }
}
\references{
  A parallel effort in metabolite profiling data sharing:
  \url{http://metlin.scripps.edu/}
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{xcmsRaw}}
}
\keyword{classes}
