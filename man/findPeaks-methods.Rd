\name{findPeaks-methods}
\docType{methods}
\alias{findPeaks-methods}
\alias{findPeaks}
\alias{findPeaks,xcmsRaw-method}
\title{Feature detection for GC/MS and LC/MS Data - methods}
\description{
  A number of peak pickers exist in XCMS. \code{findPeaks}
  is the generic method.
}
\section{Methods}{
\describe{
\item{object = "xcmsRaw"}{
  \code{
    findPeaks(object, ...)
  }
}
}}
\arguments{
  \item{object}{\code{\link{xcmsRaw-class}} object}
  \item{method}{Method to use for peak detection. See details.}
  \item{...}{Optional arguments to be passed along}
}
\details{
  Different algorithms can be used by specifying them with the
  \code{method} argument. For example to use the matched filter
  approach described by Smith et al (2006) one would use:
  \code{findPeaks(object, method="matchedFilter")}.

  Further arguments given by \code{...} are
  passed through to the function implementing
  the \code{method}. 
  
  A character vector of \emph{nicknames} for the 
  algorithms available is returned by 
  \code{getOption("BioC")$xcms$findPeaks.methods}.
  If the nickname of a method is called "centWave",
  the help page for that specific method can
  be accessed with \code{?findPeaks.centWave}. 
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
  \item{maxo}{
    maximum intensity of original (raw) peak
  }
  
  and additional columns depending on the choosen method.

}
\seealso{
  \code{\link{findPeaks.matchedFilter}}
  \code{\link{findPeaks.centWave}}
  \code{\link{xcmsRaw-class}}
}
\keyword{methods}
