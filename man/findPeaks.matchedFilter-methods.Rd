\name{findPeaks.matchedFilter-methods}
\docType{methods}
\alias{findPeaks.matchedFilter}
\alias{findPeaks.matchedFilter,xcmsRaw-method}
\title{Feature detection in the chromatographic time domain}

\description{
  Find peaks in extracted the chromatographic time domain of the profile
  matrix.
}

\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks.matchedFilter(object, fwhm = 30, sigma = fwhm/2.3548, max = 5, 
    snthresh = 10, step = 0.1, steps = 2, mzdiff = 0.8 - step*steps, 
    index = FALSE, sleep = 0)
  }
}
}}
\arguments{
  \item{object}{\code{xcmsRaw} object}
  \item{fwhm}{
    full width at half maximum of matched filtration gaussian model
    peak
  }
  \item{sigma}{
    standard deviation of matched filtration model peak
  }
  \item{max}{
    maximum number of peaks per extracted ion chromatogram
  }
  \item{snthresh}{signal to noise ratio cutoff}
  \item{step}{step size to use for profile generation}
  \item{steps}{number of steps to merge prior to filtration}
  \item{mzdiff}{
    minimum difference in m/z for peaks with overlapping retention
    times
  }
  \item{index}{
    return indicies instead of values for m/z and retention times
  }
  \item{sleep}{
    number of seconds to pause between plotting peak finding cycles
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
  \item{intf}{
    integrated area of filtered peak
  }
  \item{maxo}{
    maximum intensity of original (raw) peak
  }
  \item{maxf}{
    maximum intensity of filtered peak
  }
  \item{i}{
    rank of peak identified in merged EIC (<= \code{max})
  }
  \item{sn}{
    signal to noise ratio of the peak
  }
  
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\seealso{
  \code{\link{findPeaks-methods}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
