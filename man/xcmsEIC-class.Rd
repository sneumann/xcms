\name{xcmsEIC-class}
\docType{class}
\alias{xcmsEIC-class}
\alias{mzrange,xcmsEIC-method}
\alias{mzrange}
\alias{rtrange,xcmsEIC-method}
\alias{rtrange}
\alias{show,xcmsEIC-method}

\title{Class xcmsEIC, a class for multi-sample extracted ion chromatograms}
\description{
  This class is used to store and plot parallel extracted ion
  chromatograms from multiple sample files. It integrates with the
  \code{xcmsSet} class to display peak area integrated during peak
  identification or fill-in.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{getEIC}} method of
  the \code{xcmsSet} class. Objects can also be created by calls
  of the form \code{new("xcmsEIC", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{eic}:}{
      list containing named entries for every sample. for each entry, a 
      list of two column EIC matricies with retention time and intensity
    }
    \item{\code{mzrange}:}{
      two column matrix containing starting and ending m/z for each EIC
    }
    \item{\code{rtrange}:}{
      two column matrix containing starting and ending time for
      each EIC
    }
    \item{\code{rt}:}{
      either \code{"raw"} or \code{"corrected"} to specify retention
      times contained in the object
    }
    \item{\code{groupnames}:}{
      group names from \code{xcmsSet} object used to generate EICs
    }
  }
}
\section{Methods}{
  \describe{
    \item{\link{groupnames}}{
      \code{signature(object = "xcmsEIC")}: get \code{groupnames} slot
    }
    \item{mzrange}{
      \code{signature(object = "xcmsEIC")}: get \code{mzrange} slot
    }
    \item{\link[xcms:plot-methods]{plot}}{
      \code{signature(x = "xcmsEIC")}: plot the extracted ion
      chromatograms
    }
    \item{rtrange}{
      \code{signature(object = "xcmsEIC")}: get \code{rtrange} slot
    }
    \item{\link{sampnames}}{
      \code{signature(object = "xcmsEIC")}: get sample names
    }
  }
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  No notes yet.
}
\seealso{
  \code{\link{getEIC}}
}
\keyword{classes}
