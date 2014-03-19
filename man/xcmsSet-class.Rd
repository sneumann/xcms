\name{xcmsSet-class}
\docType{class}
\alias{xcmsSet-class}
\alias{filepaths<-,xcmsSet-method}
\alias{filepaths<-}
\alias{filepaths,xcmsSet-method}
\alias{filepaths}
\alias{groupidx<-,xcmsSet-method}
\alias{groupidx<-}
\alias{groupidx,xcmsSet-method}
\alias{groupidx}
\alias{groups<-,xcmsSet-method}
\alias{groups<-}
\alias{groups,xcmsSet-method}
\alias{groups}
\alias{peaks<-,xcmsSet-method}
\alias{peaks<-}
\alias{peaks,xcmsSet,missing-method}
\alias{peaks}
\alias{phenoData}
\alias{phenoData,xcmsSet-method}
\alias{phenoData<-}
\alias{phenoData<-,xcmsSet-method}
\alias{profinfo<-,xcmsSet-method}
\alias{profinfo<-}
\alias{profinfo,xcmsSet-method}
\alias{profinfo}
\alias{progressCallback}
\alias{progressCallback,xcmsSet-method}
\alias{progressCallback<-}
\alias{progressCallback<-,xcmsSet-method}
\alias{sampclass<-,xcmsSet-method}
\alias{sampclass<-}
\alias{sampclass,xcmsSet-method}
\alias{sampclass}
\alias{sampnames<-,xcmsSet-method}
\alias{sampnames<-}
\alias{show,xcmsSet-method}

\title{Class xcmsSet, a class for preprocessing peak data}
\description{
  This class transforms a set of peaks from multiple LC/MS or GC/MS
  samples into a matrix of preprocessed data. It groups the peaks
  and does nonlinear retention time correction without internal
  standards. It fills in missing peak values from raw data. Lastly,
  it generates extracted ion chromatograms for ions of interest.
}
\section{Objects from the Class}{
  Objects can be created with the \code{\link{xcmsSet}} constructor
  which gathers peaks from a set NetCDF files. Objects can also be
  created by calls of the form \code{new("xcmsSet", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{peaks}:}{
      matrix containing peak data
    }
    \item{\code{filled}:}{
      a vector with peak indices of peaks which have been added by a
      \code{\link{fillPeaks}} method,
    }
    \item{\code{groups}:}{
      matrix containing statistics about peak groups
    }
    \item{\code{groupidx}:}{
      list containing indices of peaks in each group
    }
    \item{\code{phenoData}:}{
      a data frame containing the experimental design factors
    }
    \item{\code{rt}:}{
      list containing two lists, \code{raw} and \code{corrected},
      each containing retention times for every scan of every sample
    }
    \item{\code{filepaths}:}{
      character vector with absolute path name of each NetCDF file
    }
    \item{\code{profinfo}:}{
      list containing two values, \code{method} - profile generation
      method, and \code{step} - profile m/z step size
    }
	\item{dataCorrection}{
	  logical vector filled if the waters Lock mass correction
		parameter is used.
	}
    \item{\code{polarity}:}{
      a string ("positive" or "negative" or NULL) describing whether only
      positive or negative scans have been used reading the raw data.
    }
    \item{\code{progressInfo}:}{
       progress informations for some xcms functions (for GUI)
    }
    \item{\code{progressCallback}:}{
       function to be called, when progressInfo changes (for GUI)
    }
  }
}
\section{Methods}{
  \describe{
    \item{\link[xcms:c.xcmsSet]{c}}{
      \code{signature("xcmsSet")}: combine objects together
    }
    \item{filepaths<-}{
      \code{signature(object = "xcmsSet")}: set \code{filepaths} slot
    }
    \item{filepaths}{
      \code{signature(object = "xcmsSet")}: get \code{filepaths} slot
    }
    \item{\link{diffreport}}{
      \code{signature(object = "xcmsSet")}: create report of
      differentially regulated ions including EICs
    }
    \item{\link{fillPeaks}}{
      \code{signature(object = "xcmsSet")}: fill in peak data for
      groups with missing peaks
    }
    \item{\link{getEIC}}{
      \code{signature(object = "xcmsSet")}: get list of EICs for
      each sample in the set
    }
    \item{groupidx<-}{
      \code{signature(object = "xcmsSet")}: set \code{groupidx} slot
    }
    \item{groupidx}{
      \code{signature(object = "xcmsSet")}: get \code{groupidx} slot
    }
    \item{\link{groupnames}}{
      \code{signature(object = "xcmsSet")}: get textual names for
      peak groups
    }
    \item{groups<-}{
      \code{signature(object = "xcmsSet")}: set \code{groups} slot
    }
    \item{groups}{
      \code{signature(object = "xcmsSet")}: get \code{groups} slot
    }
    \item{\link{groupval}}{
      \code{signature(object = "xcmsSet")}: get matrix of values
      from peak data with a row for each peak group
    }
    \item{\link{group}}{
      \code{signature(object = "xcmsSet")}: find groups of peaks
      across samples that share similar m/z and retention times
    }
    \item{peaks<-}{
      \code{signature(object = "xcmsSet")}: set \code{peaks} slot
    }
    \item{peaks}{
      \code{signature(object = "xcmsSet")}: get \code{peaks} slot
    }
    \item{\link{plotrt}}{
      \code{signature(object = "xcmsSet")}: plot retention time
      deviation profiles
    }
    \item{profinfo<-}{
      \code{signature(object = "xcmsSet")}: set \code{profinfo} slot
    }
    \item{profinfo}{
      \code{signature(object = "xcmsSet")}: get \code{profinfo} slot
    }
    \item{\link{retcor}}{
      \code{signature(object = "xcmsSet")}: use initial grouping
      of peaks to do nonlinear loess retention time correction
    }
    \item{sampclass<-}{
      \code{signature(object = "xcmsSet")}: \strong{DEPRECATED}. If used, the
      experimental design will be replaced with a data frame with a
      single column matching the supplied factor.
    }
    \item{sampclass}{
      \code{signature(object = "xcmsSet")}: get the interaction of the
      experimental design factors
    }
    \item{phenoData<-}{
      \code{signature(object = "xcmsSet")}: set the \code{phenoData} slot
    }
    \item{phenoData}{
      \code{signature(object = "xcmsSet")}: get the \code{phenoData} slot
    }
    \item{progressCallback<-}{
      \code{signature(object = "xcmsSet")}: set the \code{progressCallback} slot
    }
    \item{progressCallback}{
      \code{signature(object = "xcmsSet")}: get the \code{progressCallback} slot
    }
    \item{sampnames<-}{
      \code{signature(object = "xcmsSet")}: set rownames in the
      \code{phenoData} slot
    }
    \item{\link{sampnames}}{
      \code{signature(object = "xcmsSet")}: get rownames in the
      \code{phenoData} slot
    }
    \item{\link[xcms:split.xcmsSet]{split}}{
      \code{signature("xcmsSet")}: divide into a list of objects
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
  \code{\link{xcmsSet}}
}
\keyword{classes}
