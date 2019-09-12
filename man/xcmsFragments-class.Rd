\name{xcmsFragments-class}
\docType{class}
\alias{xcmsFragments-class}
\alias{plotTree}
\alias{plotTree,xcmsFragments-method}
\alias{show,xcmsFragments-method}

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
      (fragment mass-to-charge), intensity (peak intensity extracted
      from the original \code{xcmsSet}), sample (the index of the rawData-file).
    }
    \item{\code{MS2spec}:}{This is a list of matrixes. Each matrix in the list is a single collected spectra from \code{collect}. The column ID's are mz, intensity, and full width half maximum(fwhm). The fwhm column is only relevant if the spectra came from profile data.}
    \item{\code{specinfo}:}{This is a matrix with reference data for the spectra in MS2spec. The column id's are preMZ, AccMZ, rtmin, rtmax, ref, CollisionEnergy. The preMZ is precursor mass from the MS1 scan. This mass is given by the XML file. With some instruments this mass is only given as nominal mass, therefore a AccMZ is given which is a weighted average mass from the MS1 scan of the collected spectra. The retention time is given by rtmin and rtmax. The ref column is a pointer to the MS2spec matrix spectra. The collisionEnergy column is the collision Energy for the spectra.}
  }
}
\section{Methods}{
  \describe{
    \item{\link{collect}}{
      \code{signature(object = "xcmsFragments")}: gets a xcmsSet-object, collects ms1-peaks from it and the msn-peaks from the corresponding xcmsRaw-files.
     }
    \item{plotTree}{
      \code{signature(object = "xcmsFragments")}: prints a (text based) pseudo-tree of the peaktable to display the dependencies of the peaks among each other.
    }
    \item{\link{show}}{
      \code{signature(object = "xcmsFragments")}: print a human-readable
      description of this object to the console.
    }
  }
}
\author{S. Neumann, J. Kutzera}
\seealso{
  \code{\link{xcmsRaw}}
}
\keyword{classes}
