\name{xcmsSet}
\alias{xcmsSet}
\title{Constructor for xcmsSet objects which finds peaks in NetCDF/mzXML files}
\description{
  This function handles the construction of xcmsSet objects. It
  finds peaks in batch mode and pre-sorts files from subdirectories
  into different classes suitable for grouping.
}
\usage{
xcmsSet(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
        profmethod = "bin", profparam = list(),
        polarity = NULL, lockMassFreq=FALSE,
	mslevel=NULL, nSlaves=0, progressCallback=NULL,
        scanrange = NULL, BPPARAM = bpparam(),
        stopOnError = TRUE, ...)
}
\arguments{
  \item{files}{path names of the NetCDF/mzXML files to read}
  \item{snames}{sample names. By default the file name without extension
    is used.}
  \item{sclass}{sample classes.}
  \item{phenoData}{\code{data.frame} or \code{AnnotatedDataFrame}
    defining the sample names and classes and other sample related
    properties. If not provided, the argument \code{sclass} or the
    subdirectories in which the samples are stored will be
    used to specify sample grouping.}
  \item{profmethod}{Method to use for profile generation. Supported
    values are \code{"bin"}, \code{"binlin"}, \code{"binlinbase"} and
    \code{"intlin"} (for methods \code{\link{profBin}},
    \code{\link{profBinLin}}, \code{\link{profBinLinBase}} and
    \code{\link{profIntLin}}, respectively). See help on
    \code{\link{profBin}} for a complete list of available methods and
    their supported parameters.}
  \item{profparam}{parameters to use for profile generation.}
  \item{polarity}{filter raw data for positive/negative scans}
	\item{lockMassFreq}{Performs correction for Waters LockMass function}
    \item{mslevel}{perform peak picking on data of given mslevel}
    \item{nSlaves}{\emph{DEPRECATED}, use \code{BPPARAM} argument
      instead.}
  \item{progressCallback}{function to be called, when progressInfo changes (useful for GUIs)}
  \item{scanrange}{scan range to read}
  \item{BPPARAM}{a \code{BiocParallel} parameter object to control how
    and if parallel processing should be performed. Such objects can be
    created by the \code{\link{SerialParam}},
    \code{\link{MulticoreParam}} or
    \code{\link{SnowParam}} functions.
  }
  \item{stopOnError}{Logical specifying whether the feature detection
    call should stop on the first encountered error (the default), or
    whether feature detection is performed in all files regardless
    eventual failures for individual files in which case all errors are
    reported as \code{warnings}.
  }
  \item{\dots}{
    further arguments to the \code{findPeaks} method of the
    \code{xcmsRaw} class
  }
}
\details{
  The default values of the \code{files}, \code{snames}, \code{sclass}, and
  \code{phenoData} arguments cause the function to recursively search
  for readable files. The filename without extention is used for the
  sample name. The subdirectory path is used for the sample class.
  If the files contain both positive and negative spectra, the polarity
  can be selected explicitly. The default (NULL) is to read all scans.

  If \code{phenoData} is provided, it is stored to the \code{phenoData}
  slot of the returned \code{xcmsSet} class. If that \code{data.frame}
  contains a column named \dQuote{class}, its content will be returned
  by the \code{\link{sampclass}} method and thus be used for the
  group/class assignment of the individual files (e.g. for peak grouping
  etc.). For more details see the help of the \code{\link{xcmsSet-class}}.

  The step size (in m/z) to use for profile generation can be submitted
  either using the \code{profparam} argument
  (e.g. \code{profparam=list(step=0.1)}) or by submitting
  \code{step=0.1}. By specifying a value of \code{0} the profile matrix
  generation can be skipped.

  The feature/peak detection algorithm can be specified with the
  \code{method} argument which defaults to the \code{"matchFilter"}
  method (\code{\link{findPeaks.matchedFilter}}). Possible values are
  returned by \code{getOption("BioC")$xcms$findPeaks.methods}.

  The lock mass correction allows for the lock mass scan to be added back in
  with the last working scan. This correction gives better reproducibility
  between sample sets.
}
\value{
  A \code{xcmsSet} object.
}
\author{Colin A. Smith, \email{csmith@scripps.edu}}
\note{
  The arguments \code{profmethod} and \code{profparam} have no influence
  on the feature/peak detection. The step size parameter \code{step} for
  the profile generation in the \code{\link{findPeaks.matchedFilter}}
  peak detection algorithm can be passed using the \code{...}.
}
\seealso{
  \code{\link{xcmsSet-class}},
  \code{\link{findPeaks}},
  \code{\link{profStep}},
  \code{\link{profMethod}},
  \code{\link{profBin}},
  \code{\link{xcmsPapply}}
}
\keyword{file}
