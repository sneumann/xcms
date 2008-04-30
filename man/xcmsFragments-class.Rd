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
      (fragment mass-to-charge), intensity (peak intensity as reported
      by XXX), sample (The ID of the rawData-file).
    }
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
\references{
  A parallel effort in metabolite profiling data sharing:
  \url{http://metlin.scripps.edu/}
}
\author{S. Neumann, J. Kutzera}
\note{
  No notes yet.
}
\seealso{
  \code{\link{xcmsRaw}}
}
\keyword{classes}
