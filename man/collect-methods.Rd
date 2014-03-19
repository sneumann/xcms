\name{collect-methods}
\docType{methods}
\alias{collect-methods}
\alias{collect}
\alias{collect,xcmsRaw-method}
\alias{collect,xcmsFragments-method}
\title{Collect MS^n peaks into xcmsFragments}
\description{
  Collecting Peaks into \code{\linkS4class{xcmsFragments}}s from several
  MS-runs using \code{\linkS4class{xcmsSet}} and
  \code{\linkS4class{xcmsRaw}}.
}
\section{Methods}{
\describe{
\item{object = "xcmsFragments"}{
  \code{
    collect(object, ...)
   }
}
}}
\arguments{
  \item{object}{(empty) \code{\link{xcmsFragments-class}} object}
  \item{xs}{A \code{\link{xcmsSet-class}} object which contains
    picked ms1-peaks from several experiments}
  \item{compMethod}{("floor", "round", "none"): compare-method which is
    used to find the parent peak of a MSnpeak through comparing the
    MZ-values of the  MS1peaks with the MSnParentPeaks.}
  \item{snthresh, mzgap, uniq}{these are the parameters for the
    getspec-peakpicker included in xcmsRaw.}
}
\details{
  After running collect(xFragments,xSet) The peak table of the
  xcmsFragments includes the ms1Peaks from all experiments stored in a
  xcmsSet-object. Further it contains the relevant msN-peaks from the
  xcmsRaw-objects, which were created temporarily with the paths in
  xcmsSet.
}

\value{
  A matrix with columns:
  \item{peakID}{
    unique identifier of every peak
  }
  \item{MSnParentPeakID}{
    PeakID of the parent peak of a msLevel>1 - peak, it is 0 if the peak is msLevel 1.
  }
  \item{msLevel}{
    The msLevel of the peak.
  }
  \item{rt}{
    retention time of the peak midpoint
  }
  \item{mz}{
    the mz-Value of the peak
  }
  \item{intensity}{
    the intensity of the peak
  }
 \item{sample}{
    the number of the sample from the xcmsSet
  }
 \item{GroupPeakMSn}{Used for grouped xcmsSet groups}
 \item{CollisionEnergy}{The collision energy of the fragment}
}
\keyword{methods}
