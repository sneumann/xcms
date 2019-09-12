\name{xcmsSet-class}
\docType{class}
\alias{xcmsSet-class}
\alias{[,xcmsSet-method}
\alias{[,xcmsSet,ANY,ANY,ANY-method}
\alias{$<-,xcmsSet-method}
\alias{$,xcmsSet-method}
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
\alias{mslevel,xcmsSet-method}
\alias{mslevel}
\alias{peaks<-,xcmsSet-method}
\alias{peaks<-}
\alias{peaks,xcmsSet-method}
\alias{peaks}
\alias{phenoData}
\alias{phenoData,xcmsSet-method}
\alias{phenoData<-}
\alias{phenoData<-,xcmsSet-method}
\alias{phenoData<-,xcmsSet,ANY-method}
\alias{profinfo<-,xcmsSet-method}
\alias{profinfo<-}
\alias{profinfo,xcmsSet-method}
\alias{profinfo}
\alias{profMethod,xcmsSet-method}
\alias{profStep,xcmsSet-method}
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
\alias{scanrange,xcmsSet-method}
\alias{scanrange}
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
    \item{peaks}{
      \code{matrix} containing peak data.
    }
    \item{filled}{
      A vector with peak indices of peaks which have been added by a
      \code{\link{fillPeaks}} method.
    }
    \item{groups}{
      Matrix containing statistics about peak groups.
    }
    \item{groupidx}{
      List containing indices of peaks in each group.
    }
    \item{phenoData}{
      A \code{data.frame} containing the experimental design factors.
    }
    \item{rt}{
      \code{list} containing two lists, \code{raw} and \code{corrected},
      each containing retention times for every scan of every sample.
    }
    \item{filepaths}{
      Character vector with absolute path name of each NetCDF file.
    }
    \item{profinfo}{
      \code{list} containing the values \code{method} - profile generation
      method, and \code{step} - profile m/z step size and eventual
      additional parameters to the profile function.
    }
    \item{dataCorrection}{
      \code{logical} vector filled if the waters Lock mass correction
      parameter is used.
    }
    \item{polarity}{
      A string ("positive" or "negative" or NULL) describing whether only
      positive or negative scans have been used reading the raw data.
    }
    \item{progressInfo}{
       Progress informations for some xcms functions (for GUI).
    }
    \item{progressCallback}{
       Function to be called, when progressInfo changes (for GUI).
    }
    \item{mslevel}{
      Numeric representing the MS level on which the peak picking was
      performed (by default on MS1). This slot should be accessed
      through its getter method \code{mslevel}.
    }
    \item{scanrange}{
      Numeric of length 2 specifying the scan range (or \code{NULL} for
      the full range). This slot should be accessed through its getter
      method \code{scanrange}. The scan range provided in this slot
      represents the scans to which the whole raw data is subsetted.
    }
    \item{.processHistory}{
      Internal slot to be used to keep track of performed processing
      steps. This slot should not be directly accessed by the user.
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
    \item{\link{getXcmsRaw}}{
      \code{signature(object = "xcmsSet", sampleidx = 1,
	profmethod = profMethod(object), profstep = profStep(object),
	profparam=profinfo(object), mslevel = NULL, scanrange = NULL,
	rt=c("corrected", "raw"), BPPARAM = bpparam())}: read the raw
      data for one or more files in the \code{xcmsSet} and return
      it. The default parameters will apply all settings used in the
      original \code{\link{xcmsSet}} call to generate the \code{xcmsSet}
      object to be applied also to the raw data. Parameter
      \code{sampleidx} allows to specify which raw file(s) should be
      loaded. Argument \code{BPPARAM} allows to setup parallel
      processing.
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
    \item{mslevel}{
      Getter method for the \code{mslevel} slot.
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
    \item{profMethod}{
      \code{signature(object = "xcmsSet")}: extract the method used to
      generate the profile matrix.
    }
    \item{profStep}{
      \code{signature(object = "xcmsSet")}: extract the profile step
      used for the generation of the profile matrix.
    }
    \item{\link{retcor}}{
      \code{signature(object = "xcmsSet")}: use initial grouping
      of peaks to do nonlinear loess retention time correction
    }
    \item{sampclass<-}{
      \code{signature(object = "xcmsSet")}: Replaces the column
      \dQuote{class} in the \code{phenoData} slot. See details for more information.
    }
    \item{sampclass}{
      \code{signature(object = "xcmsSet")}: Returns the content of the
      column \dQuote{class} from the \code{phenoData} slot or, if not
      present, the interaction of the experimental design factors
      (i.e. of the \code{phenoData} \code{data.frame}). See details for
      more information.
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
    \item{scanrange}{
      Getter method for the \code{scanrange} slot. See scanrange slot
      description above for more details.
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
      \code{signature("xcmsSet")}: divide the xcmsSet into a list of
      xcmsSet objects depending on the provided factor. Note that only
      peak data will be preserved, i.e. eventual peak grouping information
      will be lost.
    }
    \item{\code{object$name}, \code{object$name<-value}}{
      Access and set \code{name} column in \code{phenoData}
    }
    \item{\code{object[, i]}}{
      Conducts subsetting of a \code{xcmsSet} instance. Only subsetting
      on columns, i.e. samples, is supported. Subsetting is performed on
      all slots, also on \code{groups} and \code{groupidx}. Parameter
      \code{i} can be an integer vector, a logical vector or a character
      vector of sample names (matching \code{sampnames}).
    }
  }
}
\details{
  The \code{phenoData} slot (and \code{phenoData} parameter in the
  \code{\link{xcmsSet}} function) is intended to contain a \code{data.frame} describing
  all experimental factors, i.e. the samples along with their
  properties. If this \code{data.frame} contains a column named
  \dQuote{class}, this will be returned by the \code{sampclass} method
  and will thus be used by all methods to determine the sample
  grouping/class assignment (e.g. to define the colors in various plots
  or for the \code{\link{group}} method).

  The \code{sampclass<-} method adds or replaces the \dQuote{class}
  column in the \code{phenoData} slot. If a \code{data.frame} is
  submitted to this method, the interaction of its columns will be
  stored into the \dQuote{class} column.

  Also, similar to other classes in Bioconductor, the \code{$} method
  can be used to directly access all columns in the \code{phenoData}
  slot (e.g. use \code{xset$name} on a \code{xcmsSet} object called
  \dQuote{xset} to extract the values from a column named \dQuote{name} in the \code{phenoData} slot).
}
\author{Colin A. Smith, \email{csmith@scripps.edu}, Johannes Rainer \email{johannes.rainer@eurac.edu}}
\seealso{
  \code{\link{xcmsSet}}
}
\keyword{classes}
