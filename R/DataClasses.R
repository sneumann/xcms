## All class definitions should go in here.
#' @include AllGenerics.R

############################################################
## Class unions
setClassUnion("characterOrNULL", c("character", "NULL"))
setClassUnion("logicalOrNumeric", c("logical", "numeric"))
##setClassUnion("ANYorNULL", c("ANY", "NULL"))


############################################################
## xcmsSet
##
setClass("xcmsSet",
         representation = representation(peaks = "matrix",
                                         groups = "matrix",
                                         groupidx = "list",
                                         filled="numeric",
                                         phenoData = "data.frame",
                                         rt = "list",
                                         filepaths = "character",
                                         profinfo = "list",
                                         dataCorrection="numeric",
                                         polarity = "character",
                                         progressInfo = "list",
                                         progressCallback="function",
                                         mslevel = "numeric",
                                         scanrange = "numeric",
                                         .processHistory = "list"),
         prototype = prototype(peaks = matrix(nrow = 0, ncol = 0),
                               groups = matrix(nrow = 0, ncol = 0),
                               groupidx = list(),
                               filled = integer(0),
                               phenoData = data.frame(),
                               rt = list(),
                               filepaths = character(0),
                               profinfo = vector("list"),
                               dataCorrection=integer(0),
                               polarity = character(0),
                               progressInfo = list(),
                               mslevel = numeric(0),
                               scanrange= numeric(0),
                               progressCallback = function(progress) NULL,
                               .processHistory = list()),
         validity = function(object) {
             msg <- character()
             ## Check if all slots are present.
             slNames <- slotNames(object)
             missingSlots <- character()
             for (i in 1:length(slNames)) {
                 if (!.hasSlot(object, slNames[i]))
                     missingSlots <- c(missingSlots, slNames[i])
             }
             if (length(missingSlots) > 0)
                 msg <- c(msg, paste0("This xcmsSet lacks slot(s): ",
                                      paste(missingSlots, collapse = ","),
                                      ". Please update the object using",
                                      " the 'updateObject' method."))
             ## Check the .processHistory slot.
             if (!any(missingSlots == ".processHistory")) {
                 inh <- unlist(lapply(object@.processHistory,
                                      FUN = function(z) {
                                          return(inherits(z, "ProcessHistory"))
                                      }))
                 if (!all(inh))
                     msg <- c(msg,
                              paste0("Slot '.processHistory' should",
                                     " only contain 'ProcessHistory'",
                                     " objects!"))
             }
             if (length(msg))
                 return(msg)
             return(TRUE)
         })

############################################################
## xcmsEIC
setClass("xcmsEIC",
         representation(eic = "list",
                        mzrange = "matrix",
                        rtrange = "matrix",
                        rt = "character",
                        groupnames = "character"),
         prototype(eic = list(),
                   mzrange = matrix(nrow = 0, ncol = 0),
                   rtrange = matrix(nrow = 0, ncol = 0),
                   rt = character(0),
                   groupnames = character(0)))

############################################################
## xcmsFragments
setClass("xcmsFragments",
         representation(peaks = "matrix",
                        MS2spec = "list",
                        specinfo = "matrix"
                        ##, pipeline = "xcmsRawPipeline"
                        ),
         prototype(peaks = matrix(nrow = 0, ncol = 6),
                   MS2spec=NULL,
                   specinfo=NULL
                   ##, pipeline = new("xcmsRawPipeline")
                   ))

############################################################
## xcmsSource
setClass("xcmsSource", representation("VIRTUAL"))
## If given an xcmsSource object, simply return it unchanged
setMethod("xcmsSource", "xcmsSource", function(object) object)

############################################################
## xcmsFileSource
setClass("xcmsFileSource",
         representation("character"),
         contains="xcmsSource",
         validity=function(object) {
             if (file.exists(object)) TRUE
             else paste("File not found:", object)
         })

############################################################
## xcmsRaw
setClass("xcmsRaw", representation(env = "environment",
                                   tic = "numeric",
                                   scantime = "numeric",
                                   scanindex = "integer",
                                   polarity = "factor",
                                   acquisitionNum = "integer",
                                   profmethod = "character",
                                   profparam = "list",
                                   mzrange = "numeric",
                                   gradient = "matrix",
                                   msnScanindex = "integer",
                                   msnAcquisitionNum = "integer",
                                   msnPrecursorScan = "integer",
                                   msnLevel = "integer",
                                   msnRt = "numeric",
                                   msnPrecursorMz = "numeric",
                                   msnPrecursorIntensity = "numeric",
                                   msnPrecursorCharge = "numeric",
                                   msnCollisionEnergy = "numeric",
                                   filepath = "xcmsSource",
                                   scanrange = "numeric",
                                   mslevel = "numeric"),
         prototype(env = new.env(parent=.GlobalEnv),
                   tic = numeric(0),
                   scantime = numeric(0),
                   scanindex = integer(0),
                   polarity = factor(integer(0)),
                   acquisitionNum = integer(0),
                   profmethod = "bin",
                   profparam = list(),
                   mzrange = numeric(0),
                   gradient = matrix(nrow=0, ncol=0),
                   msnScanindex = integer(0),
                   msnAcquisitionNum = integer(0),
                   msnLevel = integer(0),
                   msnRt = numeric(0),
                   msnPrecursorScan = integer(0),
                   msnPrecursorMz = numeric(0),
                   msnPrecursorIntensity = numeric(0),
                   msnPrecursorCharge = numeric(0),
                   msnCollisionEnergy = numeric(0),
                   scanrange = NULL,
                   mslevel = 1
                   ))

############################################################
## netCdfSource
setClass("netCdfSource", contains="xcmsFileSource")

############################################################
## rampSource
setClass("rampSource", contains="xcmsFileSource")

############################################################
## pwizSource
setClass("pwizSource", contains="xcmsFileSource")

############################################################
## xcmsPeaks
setClass("xcmsPeaks", contains = "matrix")

############################################################
## Processing history type statics
.PROCSTEP.UNKNOWN <- "Unknown"
.PROCSTEP.PEAK.DETECTION <- "Peak detection"
.PROCSTEP.PEAK.GROUPING <- "Peak grouping"
.PROCSTEP.RTIME.CORRECTION <- "Retention time correction"
.PROCSTEP.PEAK.FILLING <- "Missing peak filling"
.PROCSTEP.CALIBRATION <- "Calibration"
.PROCSTEPS <- c(
    .PROCSTEP.UNKNOWN,
    .PROCSTEP.PEAK.DETECTION,
    .PROCSTEP.PEAK.GROUPING,
    .PROCSTEP.RTIME.CORRECTION,
    .PROCSTEP.PEAK.FILLING,
    .PROCSTEP.CALIBRATION
)

############################################################
## ProcessHistory
#' @aliases ProcessHistory
#' 
#' @title Tracking data processing
#'
#' @description Objects of the type \code{ProcessHistory} allow to keep track
#'     of any data processing step in an metabolomics experiment. They are
#'     created by the data processing methods, such as
#'     \code{\link{findChromPeaks}} and added to the corresponding results
#'     objects. Thus, usually, users don't need to create them.
#'
#' @slot type character(1): string defining the type of the processing step.
#'     This string has to match predefined values. Use
#'     \code{\link{processHistoryTypes}} to list them.
#'
#' @slot date character(1): date time stamp when the processing step was started.
#'
#' @slot info character(1): optional additional information.
#'
#' @slot fileIndex integer of length 1 or > 1 to specify on which
#'     samples of the object the processing was performed.
#'
#' @slot error (ANY): used to store eventual calculation errors.
#'
#' @rdname ProcessHistory-class
setClass("ProcessHistory",
         slots = c(
             type = "character",
             date = "character",
             info = "character",
             fileIndex = "integer",
             error = "ANY"
         ),
         contains = "Versioned",
         prototype = prototype(
             type = .PROCSTEP.UNKNOWN,
             date = character(),
             info = character(),
             fileIndex = integer(),  ## This can be of length 1 or > 1.
             error = NULL
         ),
         validity = function(object) {
             msg <- character()
             ## check type:
             if (!any(object@type == .PROCSTEPS))
                 msg <- c(msg, paste0("Got invalid type '", object@type,
                                      "'! Allowd are: ",
                                      paste0("\"", .PROCSTEPS, "\"",
                                             collapse = ", ")))
             if (length(object@type) > 1)
                 msg <- c(msg, paste0("length of 'type' should not be ",
                                      "larger than 1!"))
             if (length(object@date) > 1)
                 msg <- c(msg, paste0("length of 'date' should not be ",
                                      "larger than 1!"))
             if (length(object@info) > 1)
                 msg <- c(msg, paste0("length of 'info' should not be ",
                                      "larger than 1!"))
             if (length(msg))
                 msg
             else
                 TRUE
         }
         )

## BasicParam class
## CentWaveParam
setClass("Param",
         representation = representation("VIRTUAL"),
         contains = c("Versioned"))
setClassUnion("ParamOrNULL", c("Param", "NULL"))

#' @aliases GenericParam
#' 
#' @title Generic parameter class
#'
#' @description The \code{GenericParam} class allows to store generic parameter
#'     information such as the name of the function that was/has to be called
#'     (slot \code{fun}) and its arguments (slot \code{args}). This object is
#'     used to track the process history of the data processings of an
#'     \code{\link{XCMSnExp}} object. This is in contrast to e.g. the
#'     \code{\link{CentWaveParam}} object that is passed to the actual
#'     processing method.
#'
#' @seealso \code{\link{processHistory}} for how to access the process history
#'     of an \code{\link{XCMSnExp}} object.
#'
#' @slot fun \code{character} specifying the function name.
#' 
#' @slot args \code{list} (ideally named) with the arguments to the
#'     function.
#' 
#' @slot .__classVersion__ the version of the class.
#' 
#' @author Johannes Rainer
#' 
#' @rdname GenericParam
#'
#' @examples
#' prm <- GenericParam(fun = "mean")
#'
#' prm <- GenericParam(fun = "mean", args = list(na.rm = TRUE))
setClass("GenericParam",
         slots = c(fun = "character",
                   args = "list"),
         contains = "Param",
         prototype = prototype(
             fun = character(),
             args = list()
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@args) > 0)
                 if (!length(object@fun) > 0)
                     msg <- c(msg, paste0("No function name specified in '@fun'",
                                          " but got '@args'"))
             if (length(object@fun) > 1)
                 msg <- c(msg, paste0("'@fun' has to be of length 1"))
             if (length(msg)) msg
             else TRUE
         }
         )

#' @aliases XProcessHistory
#' 
#' @title Tracking data processing
#'
#' @description The \code{XProcessHistory} extends the \code{ProcessHistory} by
#'     adding a slot \code{param} that allows to store the actual parameter
#'     class of the processing step.
#'
#' @slot param (Param): an object of type \code{Param} (e.g.
#'     \code{\link{CentWaveParam}}) specifying the settings of the processing
#'     step.
#' 
#' @rdname ProcessHistory-class
setClass("XProcessHistory",
         slots = c(
             param = "ParamOrNULL"
         ),
         contains = "ProcessHistory",
         prototype = prototype(
             param = NULL
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@param) > 0)
                 if(!is(object@param, "Param"))
                     msg <- c(msg,
                              paste0("Only objects from type 'Param' ",
                                     "allowed in slot '@param'! I got ",
                                     class(object@param)))
             if (length(msg)) msg
             else TRUE
         })

#' @aliases findChromPeaks
#'
#' @title Chromatographic peak detection methods.
#'
#' @description The \code{findChromPeaks} methods perform the chromatographic
#'     peak detection on LC/GC-MS data and are part of the modernized
#'     \code{xcms} user interface.
#'
#'     The implemented peak detection methods in chromatographic space are:
#'     \describe{
#'     \item{centWave}{chromatographic peak detection using the \emph{centWave}
#'     method. See \code{\link{centWave}} for more details.}
#'
#'     \item{centWave with predicted isotopes}{peak detection using a two-step
#'     centWave-based approach considering also feature isotopes. See
#'     \code{\link{centWaveWithPredIsoROIs}} for more details.}
#'
#'     \item{matchedFilter}{peak detection in chromatographic space. See
#'     \code{\link{matchedFilter}} for more details.}
#'
#'     \item{massifquant}{peak detection using the Kalman filter-based 
#'     method. See \code{\link{massifquant}} for more details.}
#'
#'     \item{MSW}{single-spectrum non-chromatography MS data peak detection.
#'     See \code{\link{MSW}} for more details.}
#'
#'     }
#' 
#' @name chromatographic-peak-detection
#' 
#' @family peak detection methods
#' 
#' @seealso \code{\link{findPeaks}} for the \emph{old} peak detection
#'     methods.
#' 
#'     \code{\link{plotChromPeaks}} to plot identified chromatographic peaks
#'     for one file.
#'
#'     \code{\link{highlightChromPeaks}} to highlight identified chromatographic
#'     peaks in an extracted ion chromatogram plot.
#' 
#' @author Johannes Rainer
NULL
#> NULL

## Main centWave documentation.
#' @title Chromatographic peak detection using the centWave method
#'
#' @aliases centWave
#'
#' @description The centWave algorithm perform peak density and wavelet based
#'     chromatographic peak detection for high resolution LC/MS data in centroid
#'     mode [Tautenhahn 2008].
#'
#' @param ppm \code{numeric(1)} defining the maximal tolerated m/z deviation in
#'     consecutive scans in parts per million (ppm) for the initial ROI
#'     definition.
#' 
#' @param peakwidth \code{numeric(2)} with the expected approximate
#'     peak width in chromatographic space. Given as a range (min, max)
#'     in seconds.
#' 
#' @param snthresh \code{numeric(1)} defining the signal to noise ratio cutoff.
#'
#' @param prefilter \code{numeric(2)}: \code{c(k, I)} specifying the prefilter
#'     step for the first analysis step (ROI detection). Mass traces are only
#'     retained if they contain at least \code{k} peaks with intensity
#'     \code{>= I}.
#' 
#' @param mzCenterFun Name of the function to calculate the m/z center of the
#'     chromatographic peak. Allowed are: \code{"wMean"}: intensity weighted
#'     mean of the peak's m/z values, \code{"mean"}: mean of the peak's m/z
#'     values, \code{"apex"}: use the m/z value at the peak apex,
#'     \code{"wMeanApex3"}: intensity weighted mean of the m/z value at the
#'     peak apex and the m/z values left and right of it and \code{"meanApex3"}:
#'     mean of the m/z value of the peak apex and the m/z values left and right
#'     of it.
#' 
#' @param integrate Integration method. For \code{integrate = 1} peak limits
#'     are found through descent on the mexican hat filtered data, for
#'     \code{integrate = 2} the descent is done on the real data. The latter
#'     method is more accurate but prone to noise, while the former is more
#'     robust, but less exact.
#' 
#' @param mzdiff \code{numeric(1)} representing the minimum difference in m/z
#'     dimension for peaks with overlapping retention times; can be negative to
#'     allow overlap.
#' 
#' @param fitgauss \code{logical(1)} whether or not a Gaussian should be fitted
#'     to each peak.
#' 
#' @param noise \code{numeric(1)} allowing to set a minimum intensity required
#'     for centroids to be considered in the first analysis step (centroids with
#'     intensity \code{< noise} are omitted from ROI detection).
#' 
#' @param verboseColumns \code{logical(1)} whether additional peak meta data
#'     columns should be returned.
#' 
#' @param roiList An optional list of regions-of-interest (ROI) representing
#'     detected mass traces. If ROIs are submitted the first analysis step is
#'     omitted and chromatographic peak detection is performed on the submitted
#'     ROIs. Each ROI is expected to have the following elements specified:
#'     \code{scmin} (start scan index), \code{scmax} (end scan index),
#'     \code{mzmin} (minimum m/z), \code{mzmax} (maximum m/z), \code{length}
#'     (number of scans), \code{intensity} (summed intensity). Each ROI should
#'     be represented by a \code{list} of elements or a single row
#'     \code{data.frame}.
#' 
#' @param firstBaselineCheck \code{logical(1)}. If \code{TRUE} continuous
#'     data within regions of interest is checked to be above the first baseline.
#' 
#' @param roiScales Optional numeric vector with length equal to \code{roiList}
#'     defining the scale for each region of interest in \code{roiList} that
#'     should be used for the centWave-wavelets.
#'
#' @details The centWave algorithm is most suitable for high resolution
#'     LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase
#'     the method identifies \emph{regions of interest} (ROIs) representing
#'     mass traces that are characterized as regions with less than \code{ppm}
#'     m/z deviation in consecutive scans in the LC/MS map. These ROIs are
#'     then subsequently analyzed using continuous wavelet transform (CWT)
#'     to locate chromatographic peaks on different scales. The first analysis
#'     step is skipped, if regions of interest are passed \emph{via} the
#'     \code{param} parameter.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports peak detection on
#'     \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
#'     objects (both defined in the \code{MSnbase} package). All of the settings
#'     to the centWave algorithm can be passed with a \code{CentWaveParam}
#'     object.
#'
#' @family peak detection methods
#' 
#' @seealso The \code{\link{do_findChromPeaks_centWave}} core API function and
#'     \code{\link{findPeaks.centWave}} for the old user interface.
#'
#' @references
#' Ralf Tautenhahn, Christoph B\"{o}ttcher, and Steffen Neumann "Highly
#' sensitive feature detection for high resolution LC/MS" \emph{BMC Bioinformatics}
#' 2008, 9:504
#' 
#' @name findChromPeaks-centWave
#' 
#' @author Ralf Tautenhahn, Johannes Rainer
NULL
#> NULL

#' @description The \code{CentWaveParam} class allows to specify all settings
#'     for a chromatographic peak detection using the centWave method. Instances
#'     should be created with the \code{CentWaveParam} constructor.
#'
#' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname findChromPeaks-centWave
#'
#' @examples
#'
#' ## Create a CentWaveParam object. Note that the noise is set to 10000 to
#' ## speed up the execution of the example - in a real use case the default
#' ## value should be used, or it should be set to a reasonable value.
#' cwp <- CentWaveParam(ppm = 20, noise = 10000)
#' ## Change snthresh parameter
#' snthresh(cwp) <- 25
#' cwp
#'
#' ## Perform the peak detection using centWave on some of the files from the
#' ## faahKO package. Files are read using the readMSData from the MSnbase
#' ## package
#' library(faahKO)
#' library(xcms)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the peak detection using the settings defined above.
#' res <- findChromPeaks(raw_data, param = cwp)
#' head(chromPeaks(res))
setClass("CentWaveParam",
         slots = c(
             ppm = "numeric",
             peakwidth = "numeric",
             snthresh = "numeric",
             prefilter = "numeric",
             mzCenterFun = "character",
             integrate = "integer",
             mzdiff = "numeric",
             fitgauss = "logical",
             noise = "numeric",
             verboseColumns = "logical",
             roiList = "list",
             firstBaselineCheck = "logical",
             roiScales = "numeric"
         ),
         contains = c("Param"),
         prototype = prototype(
             ppm = 25,
             peakwidth = c(20, 50),
             snthresh = 10,
             prefilter = c(3, 100),
             mzCenterFun = "wMean",
             integrate = 1L,
             mzdiff = -0.001,
             fitgauss = FALSE,
             noise = 0,
             verboseColumns = FALSE,
             roiList = list(),
             firstBaselineCheck = TRUE,
             roiScales = numeric()
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@ppm) != 1 | any(object@ppm < 0))
                 msg <- c(msg, paste0("'ppm' has to be positive numeric",
                                      " of length 1."))
             if (length(object@peakwidth) != 2 | any(object@peakwidth < 0))
                 msg <- c(msg, paste0("'peakwidth' has to be a numeric",
                                      " of length 2 with only positive",
                                      " values."))
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- c(msg, paste0("'snthresh' has to be a positive",
                                      " numeric of length 1."))
             if (length(object@prefilter) != 2)
                 msg <- c(msg, paste0("'prefilter' has to be a numeric",
                                      " of length 2."))
             allowed_vals <- c("wMean", "mean", "apex", "wMeanApex3",
                               "meanApex3")
             if (!(object@mzCenterFun) %in% allowed_vals)
                 msg <- c(msg, paste0("'mzCenterFun' has to be one of ",
                                      paste0("'", allowed_vals, "'",
                                             collapse = ", "), "."))
             if (!(object@integrate %in% c(1L, 2L)))
                 msg <- c(msg, paste0("'integrate' has to be either 1",
                                      " or 2."))
             if (length(object@mzdiff) != 1)
                 msg <- c(msg, paste0("'mzdiff' has to be a numeric of",
                                      " length 1."))
             if (length(object@noise) != 1)
                 msg <- c(msg, paste0("'noise' has to be a numeric of",
                                      " length 1."))
             if (length(object@fitgauss) != 1)
                 msg <- c(msg, paste0("'fitgauss' has to be a numeric of",
                                      " length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- c(msg, paste0("'verboseColumns' has to be a ",
                                      "numeric of length 1."))
             if (length(object@firstBaselineCheck) != 1)
                 msg <- c(msg, paste0("'firstBaselineCheck' has to be a",
                                      " numeric of length 1."))
             if (length(object@roiList) > 0) {
                 doHaveExpectedEls <- function(z) {
                     need <- c("scmax", "scmin", "mzmin", "mzmax", "length",
                               "intensity")
                     if (is.null(nrow(z))) {
                         OK <- all(need %in% names(z))
                     } else {
                         OK <- all(need %in% colnames(z))
                     }
                     return(OK)
                 }
                 OKs <- unlist(lapply(object@roiList, doHaveExpectedEls))
                 if (any(!OKs))
                     msg <- c(msg, paste0("'roiList' does not provide ",
                                          "all required fields!"))
             }
             if (length(object@roiScales) > 0) {   
                 if (length(object@roiList) != length(object@roiScales))
                     msg <- c(msg, paste0("'roiScales' has to have the same",
                                          " length than 'roiList'."))
             }
             if (length(msg))
                 msg
             else
                 TRUE
         })

## Main matchedFilter documentation.
#' @title Peak detection in the chromatographic time domain
#'
#' @aliases matchedFilter
#'
#' @description The \emph{matchedFilter} algorithm identifies peaks in the
#'     chromatographic time domain as described in [Smith 2006]. The intensity
#'     values are binned by cutting The LC/MS data into slices (bins) of a mass
#'     unit (\code{binSize} m/z) wide. Within each bin the maximal intensity is
#'     selected. The chromatographic peak detection is then performed in each
#'     bin by extending it based on the \code{steps} parameter to generate
#'     slices comprising bins \code{current_bin - steps +1} to
#'     \code{current_bin + steps - 1}. Each of these slices is then filtered
#'     with matched filtration using a second-derative Gaussian as the model
#'     peak shape. After filtration peaks are detected using a signal-to-ratio
#'     cut-off. For more details and illustrations see [Smith 2006].
#'
#' @param binSize \code{numeric(1)} specifying the width of the
#'     bins/slices in m/z dimension.
#' 
#' @param impute Character string specifying the method to be used for missing
#'     value imputation. Allowed values are \code{"none"} (no linear
#'     interpolation), \code{"lin"} (linear interpolation), \code{"linbase"}
#'     (linear interpolation within a certain bin-neighborhood) and
#'     \code{"intlin"}. See \code{\link{imputeLinInterpol}} for more details.
#' 
#' @param fwhm \code{numeric(1)} specifying the full width at half maximum
#'     of matched filtration gaussian model peak. Only used to calculate the
#'     actual sigma, see below.
#' 
#' @param sigma \code{numeric(1)} specifying the standard deviation (width)
#'     of the matched filtration model peak.
#' 
#' @param max \code{numeric(1)} representing the maximum number of peaks
#'     that are expected/will be identified per slice.
#' 
#' @param snthresh \code{numeric(1)} defining the signal to noise cutoff
#'     to be used in the chromatographic peak detection step.
#' 
#' @param steps \code{numeric(1)} defining the number of bins to be
#'     merged before filtration (i.e. the number of neighboring bins that will
#'     be joined to the slice in which filtration and peak detection will be
#'     performed).
#' 
#' @param mzdiff \code{numeric(1)} defining the minimum difference
#'     in m/z for peaks with overlapping retention times
#' 
#' @param index \code{logical(1)} specifying whether indicies should be
#'     returned instead of values for m/z and retention times.
#'
#' @details The intensities are binned by the provided m/z values within each
#'     spectrum (scan). Binning is performed such that the bins are centered
#'     around the m/z values (i.e. the first bin includes all m/z values between
#'     \code{min(mz) - bin_size/2} and \code{min(mz) + bin_size/2}).
#'
#'     For more details on binning and missing value imputation see
#'     \code{\link{binYonX}} and \code{\link{imputeLinInterpol}} methods.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports chromatographic peak
#'     detection on \code{\link[MSnbase]{MSnExp}} and
#'     \code{\link[MSnbase]{OnDiskMSnExp}} objects (both defined in the
#'     \code{MSnbase} package). All of the settings to the matchedFilter
#'     algorithm can be passed with a \code{MatchedFilterParam} object.
#'
#' @inheritParams imputeLinInterpol
#' 
#' @inheritParams findChromPeaks-centWave
#'
#' @family peak detection methods
#' 
#' @seealso The \code{\link{do_findChromPeaks_matchedFilter}} core API function
#'     and \code{\link{findPeaks.matchedFilter}} for the old user interface.
#'
#' @references
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' \emph{Anal. Chem.} 2006, 78:779-787.
#'
#' @author Colin A Smith, Johannes Rainer
#'
#' @name findChromPeaks-matchedFilter
NULL
#> NULL

#' @description The \code{MatchedFilterParam} class allows to specify all
#'     settings for a chromatographic peak detection using the matchedFilter
#'     method. Instances should be created with the \code{MatchedFilterParam}
#'     constructor.
#'
#' @slot .__classVersion__,binSize,impute,baseValue,distance,fwhm,sigma,max,snthresh,steps,mzdiff,index See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname findChromPeaks-matchedFilter
#'
#' @examples
#'
#' ## Create a MatchedFilterParam object. Note that we use a unnecessarily large
#' ## binSize parameter to reduce the run-time of the example.
#' mfp <- MatchedFilterParam(binSize = 5)
#' ## Change snthresh parameter
#' snthresh(mfp) <- 15
#' mfp
#'
#' ## Perform the peak detection using matchecFilter on the files from the
#' ## faahKO package. Files are read using the readMSData from the MSnbase
#' ## package
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#' ## Perform the chromatographic peak detection using the settings defined
#' ## above. Note that we are also disabling parallel processing in this
#' ## example by registering a "SerialParam"
#' register(SerialParam())
#' res <- findChromPeaks(raw_data, param = mfp)
#' head(chromPeaks(res))
setClass("MatchedFilterParam",
         slots = c(
             binSize = "numeric",
             impute = "character",
             baseValue = "numeric",
             distance = "numeric",
             fwhm = "numeric",
             sigma = "numeric",
             max = "numeric",
             snthresh = "numeric",
             steps = "numeric",
             mzdiff = "numeric",
             index = "logical"
         ),
         contains = c("Param"),
         prototype = prototype(
             binSize = 0.1,
             impute = "none",
             baseValue = numeric(),
             distance = numeric(),
             fwhm = 30,
             sigma = 12.73994,
             max = 5,
             snthresh = 10,
             steps = 2,
             mzdiff = 0.6,
             index = FALSE
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@binSize) != 1 | any(object@binSize < 0))
                 msg <- c(msg, paste0("'binSize' has to be positive",
                                      " numeric of length 1."))
             if (!any(c("none", "lin", "linbase") == object@impute))
                 msg <- c(msg,
                          paste0("Only values 'none', 'lin' and ",
                                 "'linbase' are allowed for'impute'"))
             if (length(object@baseValue) > 1)
                 msg <- c(msg, paste0("'baseValue' has to be a",
                                      " numeric of length 1."))
             if (length(object@distance) > 1)
                 msg <- c(msg, paste0("'distance' has to be a numeric",
                                      " of length 1."))
             if (length(object@fwhm) != 1)
                 msg <- c(msg, paste0("'fwhm' has to be a numeric",
                                      " of length 1."))
             if (length(object@sigma) != 1)
                 msg <- c(msg, paste0("'sigma' has to be a numeric",
                                      " of length 1."))
             if (length(object@max) != 1)
                 msg <- c(msg, paste0("'max' has to be a numeric",
                                      " of length 1."))
             if (length(object@snthresh) != 1)
                 msg <- c(msg, paste0("'snthresh' has to be a numeric",
                                      " of length 1."))
             if (length(object@steps) != 1)
                 msg <- c(msg, paste0("'steps' has to be a numeric",
                                      " of length 1."))
             if (length(object@mzdiff) != 1)
                 msg <- c(msg, paste0("'mzdiff' has to be a numeric",
                                      " of length 1."))
             if (length(object@index) != 1)
                 msg <- c(msg, paste0("'index' has to be a logical",
                                      " of length 1."))
             if (length(msg))
                 msg
             else
                 TRUE
         })


## Main massifquant documentation.
#' @title Chromatographic peak detection using the massifquant method
#'
#' @aliases massifquant
#'
#' @description Massifquant is a Kalman filter (KF)-based chromatographic peak
#'     detection for XC-MS data in centroid mode. The identified peaks
#'     can be further refined with the \emph{centWave} method (see
#'     \code{\link{findChromPeaks-centWave}} for details on centWave)
#'     by specifying \code{withWave = TRUE}.
#'
#' @param peakwidth \code{numeric(2)}. Only the first element is used by
#'     massifquant, which specifices the minimum peak length in time scans.
#'     For \code{withWave = TRUE} the second argument represents the maximum
#'     peak length subject to being greater than the mininum peak length
#'     (see also documentation of \code{\link{do_findChromPeaks_centWave}}).
#' 
#' @param prefilter \code{numeric(2)}. The first argument is only used
#'     if (\code{withWave = TRUE}); see \code{\link{findChromPeaks-centWave}}
#'     for details. The second argument specifies the minimum threshold for the
#'     maximum intensity of a chromatographic peak that must be met.
#' 
#' @param criticalValue \code{numeric(1)}. Suggested values:
#'     (\code{0.1-3.0}). This setting helps determine the the Kalman Filter
#'     prediciton margin of error. A real centroid belonging to a bonafide
#'     peak must fall within the KF prediction margin of error. Much like
#'     in the construction of a confidence interval, \code{criticalVal} loosely
#'     translates to be a multiplier of the standard error of the prediction
#'     reported by the Kalman Filter. If the peak in the XC-MS sample have
#'     a small mass deviance in ppm error, a smaller critical value might be
#'     better and vice versa.
#' 
#' @param consecMissedLimit \code{integer(1)} Suggested values: (\code{1,2,3}).
#'     While a peak is in the proces of being detected by a Kalman Filter, the
#'     Kalman Filter may not find a predicted centroid in every scan. After 1
#'     or more consecutive failed predictions, this setting informs Massifquant
#'     when to stop a Kalman Filter from following a candidate peak.
#' 
#' @param unions \code{integer(1)} set to \code{1} if apply t-test union on
#'     segmentation; set to \code{0} if no t-test to be applied on
#'     chromatographically continous peaks sharing same m/z range.
#'     Explanation: With very few data points, sometimes a Kalman Filter stops
#'     tracking a peak prematurely. Another Kalman Filter is instantiated
#'     and begins following the rest of the signal. Because tracking is done
#'     backwards to forwards, this algorithmic defect leaves a real peak
#'     divided into two segments or more. With this option turned on, the
#'     program identifies segmented peaks and combines them (merges them)
#'     into one with a two sample t-test. The potential danger of this option
#'     is that some truly distinct peaks may be merged.
#' 
#' @param checkBack \code{integer(1)} set to \code{1} if turned on; set to
#'     \code{0} if turned off. The convergence of a Kalman Filter to a peak's
#'     precise m/z mapping is very fast, but sometimes it incorporates erroneous
#'     centroids as part of a peak (especially early on). The \code{scanBack}
#'     option is an attempt to remove the occasional outlier that lies beyond
#'     the converged bounds of the Kalman Filter. The option does not directly
#'     affect identification of a peak because it is a postprocessing measure;
#'     it has not shown to be a extremely useful thus far and the default is set
#'     to being turned off.
#' 
#' @param withWave \code{logical(1)} if \code{TRUE}, the peaks identified first
#'     with Massifquant are subsequently filtered with the second step of the
#'     centWave algorithm, which includes wavelet estimation.
#'
#' @details This algorithm's performance has been tested rigorously
#'     on high resolution LC/{OrbiTrap, TOF}-MS data in centroid mode.
#'     Simultaneous kalman filters identify chromatographic peaks and calculate
#'     their area under the curve. The default parameters are set to operate on
#'     a complex LC-MS Orbitrap sample. Users will find it useful to do some
#'     simple exploratory data analysis to find out where to set a minimum
#'     intensity, and identify how many scans an average peak spans. The
#'     \code{consecMissedLimit} parameter has yielded good performance on
#'     Orbitrap data when set to (\code{2}) and on TOF data it was found best
#'     to be at (\code{1}). This may change as the algorithm has yet to be
#'     tested on many samples. The \code{criticalValue} parameter is perhaps
#'     most dificult to dial in appropriately and visual inspection of peak
#'     identification is the best suggested tool for quick optimization.
#'     The \code{ppm} and \code{checkBack} parameters have shown less influence
#'     than the other parameters and exist to give users flexibility and
#'     better accuracy.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports chromatographic peak
#'     detection on \code{\link[MSnbase]{MSnExp}} and
#'     \code{\link[MSnbase]{OnDiskMSnExp}} objects (both defined in the
#'     \code{MSnbase} package). All of the settings to the massifquant and
#'     centWave algorithm can be passed with a \code{MassifquantParam} object.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @family peak detection methods
#' 
#' @seealso The \code{\link{do_findChromPeaks_massifquant}} core API function
#'     and \code{\link{findPeaks.massifquant}} for the old user interface.
#'
#' @references
#' Conley CJ, Smith R, Torgrip RJ, Taylor RM, Tautenhahn R and Prince JT
#' "Massifquant: open-source Kalman filter-based XC-MS isotope trace feature
#' detection" \emph{Bioinformatics} 2014, 30(18):2636-43.
#'
#' @author Christopher Conley, Johannes Rainer
#'
#' @name findChromPeaks-massifquant
NULL
#> NULL

#' @description The \code{MassifquantParam} class allows to specify all
#'     settings for a chromatographic peak detection using the massifquant
#'     method eventually in combination with the centWave algorithm. Instances
#'     should be created with the \code{MassifquantParam} constructor.
#'
#' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,criticalValue,consecMissedLimit,unions,checkBack,withWave See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname findChromPeaks-massifquant
#'
#' @examples
#'
#' ## Create a MassifquantParam object.
#' mqp <- MassifquantParam()
#' ## Change snthresh parameter
#' snthresh(mqp) <- 30
#' mqp
#'
#' ## Perform the peak detection using massifquant on the files from the
#' ## faahKO package. Files are read using the readMSData from the MSnbase
#' ## package
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#' ## Perform the peak detection using the settings defined above.
#' res <- findChromPeaks(raw_data, param = mqp)
#' head(chromPeaks(res))
setClass("MassifquantParam",
         slots = c(
             ppm = "numeric",
             peakwidth = "numeric",
             snthresh = "numeric",
             prefilter = "numeric",
             mzCenterFun = "character",
             integrate = "integer",
             mzdiff = "numeric",
             fitgauss = "logical",
             noise = "numeric",
             verboseColumns = "logical",
             criticalValue = "numeric",
             consecMissedLimit = "integer",
             unions = "integer",
             checkBack = "integer",
             withWave = "logical"
         ),
         contains = c("Param"),
         prototype = prototype(
             ppm = 25,
             peakwidth = c(20, 50),
             snthresh = 10,
             prefilter = c(3, 100),
             mzCenterFun = "wMean",
             integrate = 1L,
             mzdiff = -0.001,
             fitgauss = FALSE,
             noise = 0,
             verboseColumns = FALSE,
             criticalValue = 1.125,
             consecMissedLimit = 2L,
             unions = 1L,
             checkBack = 0L,
             withWave = FALSE
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@ppm) != 1 | any(object@ppm < 0))
                 msg <- c(msg, paste0("'ppm' has to be positive numeric",
                                      " of length 1."))
             if (length(object@peakwidth) != 2 | any(object@peakwidth < 0))
                 msg <- c(msg, paste0("'peakwidth' has to be a numeric",
                                      " of length 2 with only positive",
                                      " values."))
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- c(msg, paste0("'snthresh' has to be a positive",
                                      " numeric of length 1."))
             if (length(object@prefilter) != 2)
                 msg <- c(msg, paste0("'prefilter' has to be a numeric",
                                      " of length 2."))
             allowed_vals <- c("wMean", "mean", "apex", "wMeanApex3",
                               "meanApex3")
             if (!(object@mzCenterFun) %in% allowed_vals)
                 msg <- c(msg, paste0("'mzCenterFun' has to be one of ",
                                      paste0("'", allowed_vals, "'",
                                             collapse = ", "), "."))
             if (!(object@integrate %in% c(1L, 2L)))
                 msg <- c(msg, paste0("'integrate' has to be either 1",
                                      " or 2."))
             if (length(object@mzdiff) != 1)
                 msg <- c(msg, paste0("'mzdiff' has to be a numeric of",
                                      " length 1."))
             if (length(object@noise) != 1)
                 msg <- c(msg, paste0("'noise' has to be a numeric of",
                                      " length 1."))
             if (length(object@fitgauss) != 1)
                 msg <- c(msg, paste0("'fitgauss' has to be a numeric of",
                                      " length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- c(msg, paste0("'verboseColumns' has to be a ",
                                      "numeric of length 1."))
             if (length(object@criticalValue) != 1)
                 msg <- c(msg, paste0("'criticalValue' has to be a ",
                                      "numeric of length 1."))
             if (length(object@consecMissedLimit) != 1)
                 msg <- c(msg, paste0("'consecMissedLimit' has to be a ",
                                      "numeric of length 1."))
             if (length(object@unions) != 1)
                 msg <- c(msg, paste0("'unions' has to be a ",
                                      "numeric of length 1."))
             if (object@unions != 0 & object@unions != 1)
                 msg <- c(msg, paste0("'unions' has to be either 0 or 1!"))
             if (length(object@checkBack) != 1)
                 msg <- c(msg, paste0("'checkBack' has to be a ",
                                      "numeric of length 1."))
             if (object@checkBack != 0 & object@checkBack != 1)
                 msg <- c(msg, paste0("'checkBack' has to be either 0",
                                      " or 1!"))
             if (length(object@withWave) != 1)
                 msg <- c(msg, paste0("'withWave' has to be a ",
                                      "numeric of length 1."))
             if (length(msg))
                 msg
             else TRUE
         })

## Main MSW documentation.
#' @title Single-spectrum non-chromatography MS data peak detection
#'
#' @aliases MSW
#'
#' @description Perform peak detection in mass spectrometry
#'     direct injection spectrum using a wavelet based algorithm.
#'
#' @details This is a wrapper for the peak picker in Bioconductor's
#'     \code{MassSpecWavelet} package calling
#'     \code{\link[MassSpecWavelet]{peakDetectionCWT}} and
#'     \code{\link[MassSpecWavelet]{tuneInPeakInfo}} functions. See the
#'     \emph{xcmsDirect} vignette for more information.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports peak detection on
#'     \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
#'     objects (both defined in the \code{MSnbase} package). All of the settings
#'     to the algorithm can be passed with a \code{MSWParam} object.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @family peak detection methods
#' 
#' @seealso The \code{\link{do_findPeaks_MSW}} core API function
#'     and \code{\link{findPeaks.MSW}} for the old user interface.
#'
#' @author Joachim Kutzera, Steffen Neumann, Johannes Rainer
#'
#' @name findPeaks-MSW
NULL
#> NULL

#' @description The \code{MSWParam} class allows to specify all
#'     settings for a peak detection using the MSW method. Instances should be
#'     created with the \code{MSWParam} constructor.
#'
#' @slot .__classVersion__,snthresh,verboseColumns,scales,nearbyPeak,peakScaleRange,ampTh,minNoiseLevel,ridgeLength,peakThr,tuneIn,addParams See corresponding parameter above. \code{.__classVersion__} stores the version from the class. Slots values
#' should exclusively be accessed \emph{via} the corresponding getter and
#' setter methods listed above.
#'
#' @rdname findPeaks-MSW
#'
#' @examples
#'
#' ## Create a MSWParam object
#' mp <- MSWParam()
#' ## Change snthresh parameter
#' snthresh(mp) <- 15
#' mp
#'
#' ## Loading a small subset of direct injection, single spectrum files
#' library(msdata)
#' fticrf <- list.files(system.file("fticr", package = "msdata"),
#'                     recursive = TRUE, full.names = TRUE)
#' fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
#'
#' ## Perform the MSW peak detection on these:
#' p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
#'              SNR.method = "data.mean", winSize.noise = 500)
#' fticr <- findChromPeaks(fticr, param = p)
#'
#' head(chromPeaks(fticr))
setClass("MSWParam",
         slots = c(
             snthresh = "numeric",
             verboseColumns = "logical",
             ## params from the peakDetectionCWT
             scales = "numeric",
             nearbyPeak = "logical",
             peakScaleRange = "numeric",
             ampTh = "numeric",
             minNoiseLevel = "numeric",
             ridgeLength = "numeric",
             peakThr = "numeric",
             tuneIn = "logical",
             addParams = "list"
         ),
         contains = c("Param"),
         prototype = prototype(
             snthresh = 3,
             verboseColumns = FALSE,
             scales = c(1, seq(2, 30, 2), seq(32, 64, 4)),
             nearbyPeak = TRUE,
             peakScaleRange = 5,
             ampTh = 0.01,
             minNoiseLevel = (0.01 / 3),
             ridgeLength = 24,
             peakThr = numeric(),
             tuneIn = FALSE,
             addParams = list()
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- c(msg, paste0("'snthresh' has to be a positive",
                                      " numeric of length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- c(msg, paste0("'verboseColumns' has to be a ",
                                      "numeric of length 1."))
             if (length(object@nearbyPeak) != 1)
                 msg <- c(msg, paste0("'nearbyPeak' has to be a ",
                                      "logical of length 1."))
             if (length(object@peakScaleRange) != 1 |
                 any(object@peakScaleRange < 0))
                 msg <- c(msg, paste0("'peakScaleRange' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@ampTh) != 1 | any(object@ampTh < 0))
                 msg <- c(msg, paste0("'ampTh' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@minNoiseLevel) != 1 |
                 any(object@minNoiseLevel < 0))
                 msg <- c(msg, paste0("'minNoiseLevel' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@ridgeLength) != 1 |
                 any(object@ridgeLength < 0))
                 msg <- c(msg, paste0("'ridgeLength' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@peakThr) > 1)
                 msg <- c(msg, paste0("'peakThr' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@tuneIn) != 1)
                 msg <- c(msg, paste0("'tuneIn' has to be a ",
                                      "logical of length 1."))
             if (length(msg))
                 msg
             else TRUE
         })

#' @title Two-step centWave peak detection considering also isotopes
#'
#' @aliases centWaveWithPredIsoROIs
#'
#' @description This method performs a two-step centWave-based chromatographic
#'     peak detection: in a first centWave run peaks are identified for which
#'     then the location of their potential isotopes in the mz-retention time is
#'     predicted. A second centWave run is then performed on these
#'     \emph{regions of interest} (ROIs). The final list of chromatographic
#'     peaks comprises all non-overlapping peaks from both centWave runs.
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @param maxCharge \code{integer(1)} defining the maximal isotope charge.
#'     Isotopes will be defined for charges \code{1:maxCharge}.
#'
#' @param maxIso \code{integer(1)} defining the number of isotope peaks that
#'     should be predicted for each peak identified in the first centWave run.
#'
#' @param mzIntervalExtension \code{logical(1)} whether the mz range for the
#'     predicted isotope ROIs should be extended to increase detection of low
#'     intensity peaks.
#'
#' @param snthreshIsoROIs \code{numeric(1)} defining the signal to noise ratio
#'     cutoff to be used in the second centWave run to identify peaks for
#'     predicted isotope ROIs.
#'
#' @param polarity \code{character(1)} specifying the polarity of the data.
#'     Currently not used, but has to be \code{"positive"}, \code{"negative"} or
#'     \code{"unknown"} if provided.
#'
#' @details See \code{\link{centWave}} for details on the centWave method.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports chromatographic peak
#'     detection on \code{\link[MSnbase]{MSnExp}} and
#'     \code{\link[MSnbase]{OnDiskMSnExp}} objects (both defined in the
#'     \code{MSnbase} package). All of the settings to the algorithm can be
#'     passed with a \code{CentWavePredIsoParam} object.
#'
#' @family peak detection methods
#' 
#' @seealso The \code{\link{do_findChromPeaks_centWaveWithPredIsoROIs}} core
#'     API function and \code{\link{findPeaks.centWave}} for the old user
#'     interface. \code{\link{CentWaveParam}} for the class the
#'     \code{CentWavePredIsoParam} extends.
#'
#' @name findChromPeaks-centWaveWithPredIsoROIs
#' 
#' @author Hendrik Treutler, Johannes Rainer
NULL
#> NULL

#' @description The \code{CentWavePredIsoParam} class allows to specify all
#'     settings for the two-step centWave-based peak detection considering also
#'     predicted isotopes of peaks identified in the first centWave run.
#'     Instances should be created with the \code{CentWavePredIsoParam}
#'     constructor. See also the documentation of the
#'     \code{\link{CentWaveParam}} for all methods and arguments this class
#'     inherits.
#'
#' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales,snthreshIsoROIs,maxCharge,maxIso,mzIntervalExtension,polarity See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname findChromPeaks-centWaveWithPredIsoROIs
#'
#' @examples
#'
#' ## Create a param object
#' p <- CentWavePredIsoParam(maxCharge = 4)
#' ## Change snthresh parameter
#' snthresh(p) <- 25
#' p
#'
setClass("CentWavePredIsoParam",
         slots = c(
             snthreshIsoROIs = "numeric",
             maxCharge = "integer",
             maxIso = "integer",
             mzIntervalExtension = "logical",
             polarity = "character"
         ),
         contains = c("CentWaveParam"),
         prototype = prototype(
             snthreshIsoROIs = 6.25,
             maxCharge = 3L,
             maxIso = 5L,
             mzIntervalExtension = TRUE,
             polarity = "unknown"
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@snthreshIsoROIs) != 1 |
                 any(object@snthreshIsoROIs < 0))
                 msg <- c(msg, paste0("'snthreshIsoROIs' has to be a ",
                                      "positive numeric of length 1."))
             if (length(object@maxCharge) != 1 | any(object@maxCharge < 0))
                 msg <- c(msg, paste0("'maxCharge' has to be a ",
                                      "positive integer of length 1."))
             if (length(object@maxIso) != 1 | any(object@maxIso < 0))
                 msg <- c(msg, paste0("'maxIso' has to be a ",
                                      "positive integer of length 1."))
             if (length(object@mzIntervalExtension) != 1)
                 msg <- c(msg, paste0("'mzIntervalExtension' has to be a",
                                      " logical of length 1."))
             if (length(object@polarity) != 1)
                 msg <- c(msg, paste0("'polarity' has to be a",
                                      " character of length 1."))
             if (!(object@polarity %in% c("positive", "negative", "unknown")))
                 msg <- c(msg, paste0("'polarity' has to be either ",
                                      "'positive', 'negative' or ",
                                      "'unknown'!"))
             if (length(msg))
                 msg
             else TRUE
         })


## General groupChromPeaks method.
#' @title Correspondence: Chromatographic peak grouping methods.
#'
#' @description The \code{groupChromPeaks} method(s) perform the correspondence,
#'     i.e. the grouping of chromatographic peaks within and between samples.
#'     These methods are part of the modernized \code{xcms} user interface.
#'     The resulting peak groups are referred to as (mz-rt) features and can be
#'     accessed \emph{via} the \code{\link{featureDefinitions}} method on the
#'     result object.
#'
#'     The implemented peak grouping methods are:
#'     \describe{
#' 
#'     \item{density}{peak grouping based on time dimension peak densities.
#'     See \code{\link{groupChromPeaks-density}} for more details.}
#'
#'     \item{mzClust}{high resolution peak grouping for single spectra (direct
#'     infusion) MS data. See \code{\link{groupChromPeaks-mzClust}} for more
#'     details.}
#'
#'     \item{nearest}{chromatographic peak grouping based on their proximity in
#'     the mz-rt space. See \code{\link{groupChromPeaks-nearest}} for more
#'     details.}
#' 
#' }
#' @name groupChromPeaks
#' 
#' @family peak grouping methods
#' 
#' @seealso \code{\link{group}} for the \emph{old} peak grouping methods.
#'     \code{\link{featureDefinitions}} and
#'     \code{\link{featureValues,XCMSnExp-method}} for methods to access peak
#'     grouping results.
#' 
#' @author Johannes Rainer
NULL
#> NULL

#' @title Peak grouping based on time dimension peak densities
#'
#' @description This method performs performs correspondence (chromatographic
#'     peak grouping) based on the density (distribution) of identified peaks
#'     along the retention time axis within slices of overlapping mz ranges.
#'     All peaks (from the same or from different samples) being close on the
#'     retention time axis are grouped into a feature (\emph{peak group}).
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{group}} methods. All of the settings to the algorithm
#'     can be passed with a \code{PeakDensityParam} object.
#'
#' @param sampleGroups A vector of the same length than samples defining the
#'     sample group assignments (i.e. which samples belong to which sample
#'     group).
#'
#' @param bw \code{numeric(1)} defining the bandwidth (standard deviation ot the
#'     smoothing kernel) to be used. This argument is passed to the
#'     \code{\link{density}} method.
#'
#' @param minFraction \code{numeric(1)} defining the minimum fraction of samples
#'     in at least one sample group in which the peaks have to be present to be
#'     considered as a peak group (feature).
#'
#' @param minSamples \code{numeric(1)} with the minimum number of samples in at
#'     least one sample group in which the peaks have to be detected to be
#'     considered a peak group (feature).
#'
#' @param binSize \code{numeric(1)} defining the size of the overlapping slices
#'     in mz dimension.
#'
#' @param maxFeatures \code{numeric(1)} with the maximum number of peak groups
#'     to be identified in a single mz slice.
#' 
#' @family peak grouping methods
#' 
#' @seealso The \code{\link{do_groupChromPeaks_density}} core
#'     API function and \code{\link{group.density}} for the old user interface.
#' 
#' @seealso \code{\link{plotChromPeakDensity}} to plot peak densities and
#'     evaluate different algorithm settings.
#'     \code{\link{featureDefinitions}} and
#'     \code{\link{featureValues,XCMSnExp-method}} for methods to access the
#'     features (i.e. the peak grouping results).
#'
#' @name groupChromPeaks-density
#' 
#' @author Colin Smith, Johannes Rainer
#'
#' @references
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' \emph{Anal. Chem.} 2006, 78:779-787.
NULL
#> NULL

#' @description The \code{PeakDensityParam} class allows to specify all
#'     settings for the peak grouping based on peak densities along the time
#'     dimension. Instances should be created with the \code{PeakDensityParam}
#'     constructor.
#'
#' @slot .__classVersion__,sampleGroups,bw,minFraction,minSamples,binSize,maxFeatures See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname groupChromPeaks-density
#'
#' @examples
#'
#' ## Create a PeakDensityParam object
#' p <- PeakDensityParam(binSize = 0.05)
#' ## Change hte minSamples slot
#' minSamples(p) <- 3
#' p
#'
#' ##############################
#' ## Chromatographic peak detection and grouping.
#' ##
#' ## Below we perform first a peak detection (using the matchedFilter
#' ## method) on some of the test files from the faahKO package followed by
#' ## a peak grouping using the density method.
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' 
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the chromatographic peak detection using the matchedFilter method.
#' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
#' res <- findChromPeaks(raw_data, param = mfp)
#'
#' head(chromPeaks(res))
#' ## The number of peaks identified per sample:
#' table(chromPeaks(res)[, "sample"])
#'
#' ## Performing the chromatographic peak grouping
#' fdp <- PeakDensityParam()
#' res <- groupChromPeaks(res, fdp)
#'
#' ## The definition of the features (peak groups):
#' featureDefinitions(res)
#'
#' ## Using the featureValues method to extract a matrix with the intensities of
#' ## the features per sample.
#' head(featureValues(res, value = "into"))
#' 
#' ## The process history:
#' processHistory(res)
setClass("PeakDensityParam",
         slots = c(sampleGroups = "ANY",
                   bw = "numeric",
                   minFraction = "numeric",
                   minSamples = "numeric",
                   binSize = "numeric",
                   maxFeatures = "numeric"),
         contains = "Param",
         prototype = prototype(
             sampleGroups = numeric(),
             bw = 30,
             minFraction = 0.5,
             minSamples = 1,
             binSize = 0.25,
             maxFeatures = 50),
         validity = function(object) {
             msg <- character()
             if (length(object@bw) > 1 | any(object@bw < 0))
                 msg <- c(msg, paste0("'bw' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@minFraction) > 1 | any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- c(msg, paste0("'minFraction' has to be a ",
                                      "single positive number between ",
                                      "0 and 1!"))
             if (length(object@minSamples) > 1 | any(object@minSamples < 0))
                 msg <- c(msg, paste0("'minSamples' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@binSize) > 1 | any(object@binSize < 0))
                 msg <- c(msg, paste0("'binSize' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@maxFeatures) > 1 | any(object@maxFeatures < 0))
                 msg <- c(msg, paste0("'maxFeatures' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(msg))
                 return(msg)
             else
                 return(TRUE)
         })

## Main group.mzClust documentation.
#' @title High resolution peak grouping for single spectra samples
#'
#' @description This method performs high resolution correspondence for single
#'     spectra samples.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{group}} methods. All of the settings to the algorithm
#'     can be passed with a \code{MzClustParam} object.
#'
#' @inheritParams groupChromPeaks-density
#'
#' @param ppm \code{numeric(1)} representing the relative mz error for the
#'     clustering/grouping (in parts per million).
#' 
#' @param absMz \code{numeric(1)} representing the absolute mz error for the
#'     clustering.
#' 
#' @family peak grouping methods
#' 
#' @seealso The \code{\link{do_groupPeaks_mzClust}} core API function and
#'     \code{\link{group.mzClust}} for the old user interface.
#'     \code{\link{featureDefinitions}} and
#'     \code{\link{featureValues,XCMSnExp-method}} for methods to access peak
#'     grouping results (i.e. the features).
#'
#' @name groupChromPeaks-mzClust
#'
#' @references Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill
#' and David F. Grant\cr \emph{Alignment of high resolution mass spectra:
#' development of a heuristic approach for metabolomics}.\cr Metabolomics,
#' Vol. 2, No. 2, 75-83 (2006)
NULL
#> NULL

#' @description The \code{MzClustParam} class allows to specify all
#'     settings for the peak grouping based on the \emph{mzClust} algorithm.
#'     Instances should be created with the \code{MzClustParam} constructor.
#'
#' @slot .__classVersion__,sampleGroups,ppm,absMz,minFraction,minSamples See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname groupChromPeaks-mzClust
#'
#' @examples
#'
#' ## Loading a small subset of direct injection, single spectrum files
#' library(msdata)
#' fticrf <- list.files(system.file("fticr", package = "msdata"),
#'                     recursive = TRUE, full.names = TRUE)
#' fticr <- readMSData(fticrf[1:2], msLevel. = 1, mode = "onDisk")
#'
#' ## Perform the MSW peak detection on these:
#' p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
#'              SNR.method = "data.mean", winSize.noise = 500)
#' fticr <- findChromPeaks(fticr, param = p)
#'
#' head(chromPeaks(fticr))
#'
#' ## Now create the MzClustParam parameter object: we're assuming here that
#' ## both samples are from the same sample group.
#' p <- MzClustParam(sampleGroups = c(1, 1))
#'
#' fticr <- groupChromPeaks(fticr, param = p)
#'
#' ## Get the definition of the features.
#' featureDefinitions(fticr)
setClass("MzClustParam",
         slots = c(sampleGroups = "ANY",
                   ppm = "numeric",
                   absMz = "numeric",
                   minFraction = "numeric",
                   minSamples = "numeric"),
         contains = "Param",
         prototype = prototype(
             sampleGroups = numeric(),
             ppm = 20,
             absMz = 0,
             minFraction = 0.5,
             minSamples = 1),
         validity = function(object) {
             msg <- character()
             if (length(object@ppm) > 1 | any(object@ppm < 0))
                 msg <- c(msg, paste0("'ppm' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@absMz) > 1 | any(object@absMz < 0))
                 msg <- c(msg, paste0("'absMz' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@minFraction) > 1 | any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- c(msg, paste0("'minFraction' has to be a ",
                                      "single positive number between ",
                                      "0 and 1!"))
             if (length(object@minSamples) > 1 | any(object@minSamples < 0))
                 msg <- c(msg, paste0("'minSamples' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(msg))
                 msg
             else
                 TRUE
         })

## Main group.nearest documentation.
#' @title Peak grouping based on proximity in the mz-rt space
#'
#' @description This method is inspired by the grouping algorithm of mzMine
#'     [Katajamaa 2006] and performs correspondence based on proximity of peaks
#'     in the space spanned by retention time and mz values.
#'     The method creates first a \emph{master peak list} consisting of all
#'     chromatographic peaks from the sample in which most peaks were
#'     identified, and starting from that, calculates distances to peaks from
#'     the sample with the next most number of peaks. If peaks are closer than
#'     the defined threshold they are grouped together.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{group}} methods. All of the settings to the algorithm
#'     can be passed with a \code{NearestPeaksParam} object.
#'
#' @inheritParams groupChromPeaks-density
#'
#' @param mzVsRtBalance \code{numeric(1)} representing the factor by which mz
#'     values are multiplied before calculating the (euclician) distance between
#'     two peaks.
#'
#' @param absMz \code{numeric(1)} maximum tolerated distance for mz values.
#'
#' @param absRt \code{numeric(1)} maximum tolerated distance for rt values.
#'
#' @param kNN \code{numeric(1)} representing the number of nearest neighbors
#'     to check.
#' 
#' @family peak grouping methods
#' 
#' @seealso The \code{\link{do_groupChromPeaks_nearest}} core
#'     API function and \code{\link{group.nearest}} for the old user interface.
#'     \code{\link{featureDefinitions}} and
#'     \code{\link{featureValues,XCMSnExp-method}} for methods to access
#'     peak grouping results (i.e. the features).
#'
#' @name groupChromPeaks-nearest
#'
#' @references Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for
#' processing and visualization of mass spectrometry based molecular profile
#' data. \emph{Bioinformatics} 2006, 22:634-636. 
NULL
#> NULL

#' @description The \code{NearestPeaksParam} class allows to specify all
#'     settings for the peak grouping based on the \emph{nearest} algorithm.
#'     Instances should be created with the \code{NearestPeaksParam} constructor.
#'
#' @slot .__classVersion__,sampleGroups,mzVsRtBalance,absMz,absRt,kNN See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname groupChromPeaks-nearest
#'
#' @examples
#'
#' ## Create a NearestPeaksParam object
#' p <- NearestPeaksParam(kNN = 3)
#' p
#'
#' ##############################
#' ## Chromatographic peak detection and grouping.
#' ##
#' ## Below we perform first a chromatographic peak detection (using the
#' ## matchedFilter method) on some of the test files from the faahKO package
#' ## followed by a peaks grouping using the "nearest" method.
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' 
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the peak detection using the matchedFilter method.
#' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
#' res <- findChromPeaks(raw_data, param = mfp)
#'
#' head(chromPeaks(res))
#' ## The number of peaks identified per sample:
#' table(chromPeaks(res)[, "sample"])
#'
#' ## Performing the peak grouping
#' p <- NearestPeaksParam()
#' res <- groupChromPeaks(res, param = p)
#'
#' ## The results from the peak grouping:
#' featureDefinitions(res)
#'
#' ## Using the featureValues method to extract a matrix with the intensities of
#' ## the features per sample.
#' head(featureValues(res, value = "into"))
#'
#' ## The process history:
#' processHistory(res)
setClass("NearestPeaksParam",
         slots = c(sampleGroups = "ANY",
                   mzVsRtBalance = "numeric",
                   absMz = "numeric",
                   absRt = "numeric",
                   kNN = "numeric"),
         contains = "Param",
         prototype = prototype(
             sampleGroups = numeric(),
             mzVsRtBalance = 10,
             absMz = 0.2,
             absRt = 15,
             kNN = 10),
         validity = function(object) {
             msg <- character()
             if (length(object@mzVsRtBalance) > 1 |
                 any(object@mzVsRtBalance < 0))
                 msg <- c(msg, paste0("'mzVsRtBalance' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@absMz) > 1 | any(object@absMz < 0))
                 msg <- c(msg, paste0("'absMz' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@absRt) > 1 | any(object@absRt < 0))
                 msg <- c(msg, paste0("'absRt' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@kNN) > 1 | any(object@kNN < 0))
                 msg <- c(msg, paste0("'kNN' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(msg))
                 msg
             else TRUE
         })



#' @title Alignment: Retention time correction methods.
#'
#' @description The \code{adjustRtime} method(s) perform retention time
#'     correction (alignment) between chromatograms of different samples. These
#'     methods are part of the modernized \code{xcms} user interface.
#'
#'     The implemented retention time adjustment methods are:
#'     \describe{
#'     \item{peakGroups}{retention time correction based on aligment of
#'     features (peak groups) present in most/all samples.
#'     See \code{\link{adjustRtime-peakGroups}} for more details.}
#'
#'     \item{obiwarp}{alignment based on the complete mz-rt data. This method
#'     does not require any identified peaks or defined features. See
#'     \code{\link{adjustRtime-obiwarp}} for more details.}
#'     }
#' @name adjustRtime
#' 
#' @family retention time correction methods
#' 
#' @seealso \code{\link{retcor}} for the \emph{old} retention time correction 
#'     methods.
#'     \code{\link{plotAdjustedRtime}} for visualization of alignment results.
#' 
#' @author Johannes Rainer
NULL
#> NULL

## Main retcor.peakgroups documentation.
#' @title Retention time correction based on alignment of house keeping peak
#' groups
#'
#' @description This method performs retention time adjustment based on the
#'     alignment of chromatographic peak groups present in all/most samples
#'     (hence corresponding to house keeping compounds). First the retention
#'     time deviation of these peak groups is described by fitting either a
#'     polynomial (\code{smooth = "loess"}) or a linear (
#'     \code{smooth = "linear"}) model to the data points. These models are
#'     subsequently used to adjust the retention time of each spectrum in
#'     each sample.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{group}} methods. All of the settings to the alignment
#'     algorithm can be passed with a \code{PeakGroupsParam} object.
#'
#'     The matrix with the (raw) retention times of the peak groups used
#'     in the alignment is added to the \code{peakGroupsMatrix} slot of the
#'     \code{PeakGroupsParam} object that is stored into the corresponding
#'     \emph{process history step} (see \code{\link{processHistory}} for how
#'     to access the process history).
#'
#' @param minFraction \code{numeric(1)} between 0 and 1 defining the minimum
#'     required fraction of samples in which peaks for the peak group were
#'     identified. Peak groups passing this criteria will aligned across
#'     samples and retention times of individual spectra will be adjusted
#'     based on this alignment. For \code{minFraction = 1} the peak group
#'     has to contain peaks in all samples of the experiment.
#' 
#' @param extraPeaks \code{numeric(1)} defining the maximal number of
#'     additional peaks for all samples to be assigned to a peak group (i.e.
#'     feature) for retention time correction. For a data set with 6 samples,
#'     \code{extraPeaks = 1} uses all peak groups with a total peak count
#'     \code{<= 6 + 1}. The total peak count is the total number of peaks being
#'     assigned to a peak group and considers also multiple peaks within a
#'     sample being assigned to the group.
#'
#' @param smooth character defining the function to be used, to interpolate
#'     corrected retention times for all peak groups. Either \code{"loess"} or
#'     \code{"linear"}.
#'
#' @param span \code{numeric(1)} defining the degree of smoothing (if
#'     \code{smooth = "loess"}). This parameter is passed to the internal call
#'     to \code{\link{loess}}.
#'
#' @param family character defining the method to be used for loess smoothing.
#'     Allowed values are \code{"gaussian"} and \code{"symmetric"}.See
#'     \code{\link{loess}} for more information.
#'
#' @param peakGroupsMatrix optional \code{matrix} of (raw) retention times for
#'     the peak groups on which the alignment should be performed. Each column
#'     represents a sample, each row a feature/peak group. Such a matrix is
#'     for example returned by the \code{\link{adjustRtimePeakGroups}} method.
#' 
#' @family retention time correction methods
#' 
#' @seealso The \code{\link{do_adjustRtime_peakGroups}} core
#'     API function and \code{\link{retcor.peakgroups}} for the old user
#'     interface.
#'     \code{\link{plotAdjustedRtime}} for visualization of alignment results.
#' 
#' @name adjustRtime-peakGroups
#'
#' @author Colin Smith, Johannes Rainer
#' 
#' @references
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' \emph{Anal. Chem.} 2006, 78:779-787.
NULL
#> NULL

#' @description The \code{PeakGroupsParam} class allows to specify all
#'     settings for the retention time adjustment based on \emph{house keeping}
#'     peak groups present in most samples.
#'     Instances should be created with the \code{PeakGroupsParam} constructor.
#'
#' @slot .__classVersion__,minFraction,extraPeaks,smooth,span,family,peakGroupsMatrix See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname adjustRtime-peakGroups
#'
#' @examples
#' ##############################
#' ## Chromatographic peak detection and grouping.
#' ##
#' ## Below we perform first a peak detection (using the matchedFilter
#' ## method) on some of the test files from the faahKO package followed by
#' ## a peak grouping.
#' library(faahKO)
#' library(xcms)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' 
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform the peak detection using the matchedFilter method.
#' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
#' res <- findChromPeaks(raw_data, param = mfp)
#'
#' head(chromPeaks(res))
#' ## The number of peaks identified per sample:
#' table(chromPeaks(res)[, "sample"])
#'
#' ## Performing the peak grouping using the "peak density" method.
#' p <- PeakDensityParam(sampleGroups = c(1, 1))
#' res <- groupChromPeaks(res, param = p)
#'
#' ## Perform the retention time adjustment using peak groups found in both
#' ## files.
#' fgp <- PeakGroupsParam(minFraction = 1)
#'
#' ## Before running the alignment we can evaluate which features (peak groups)
#' ## would be used based on the specified parameters.
#' pkGrps <- adjustRtimePeakGroups(res, param = fgp)
#'
#' ## We can also plot these to evaluate if the peak groups span a large portion
#' ## of the retention time range.
#' plot(x = pkGrps[, 1], y = rep(1, nrow(pkGrps)), xlim = range(rtime(res)),
#'     ylim = c(1, 2), xlab = "rt", ylab = "", yaxt = "n")
#' points(x = pkGrps[, 2], y = rep(2, nrow(pkGrps)))
#' segments(x0 = pkGrps[, 1], x1 = pkGrps[, 2],
#'     y0 = rep(1, nrow(pkGrps)), y1 = rep(2, nrow(pkGrps)))
#' grid()
#' axis(side = 2, at = c(1, 2), labels = colnames(pkGrps))
#'
#' ## Next we perform the alignment.
#' res <- adjustRtime(res, param = fgp)
#'
#' ## Any grouping information was dropped
#' hasFeatures(res)
#'
#' ## Plot the raw against the adjusted retention times.
#' plot(rtime(raw_data), rtime(res), pch = 16, cex = 0.25, col = fromFile(res))
#'
#' ## Adjusterd retention times can be accessed using
#' ## rtime(object, adjusted = TRUE) and adjustedRtime
#' all.equal(rtime(res), adjustedRtime(res))
#'
#' ## To get the raw, unadjusted retention times:
#' all.equal(rtime(res, adjusted = FALSE), rtime(raw_data))
#'
#' ## To extract the retention times grouped by sample/file:
#' rts <- rtime(res, bySample = TRUE)
setClass("PeakGroupsParam",
         slots = c(minFraction = "numeric",
                   extraPeaks = "numeric",
                   smooth = "character",
                   span = "numeric",
                   family = "character",
                   peakGroupsMatrix = "matrix"),
         contains = "Param",
         prototype = prototype(
             minFraction = 0.9,
             extraPeaks = 1,
             smooth = "loess",
             span = 0.2,
             family = "gaussian",
             peakGroupsMatrix = matrix(ncol = 0, nrow = 0)
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@minFraction) > 1 |
                 any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- c(msg, paste0("'minFraction' has to be a single",
                                      " number between 0 and 1!"))
             if (length(object@extraPeaks) > 1 |
                 any(object@extraPeaks < 0))
                 msg <- c(msg, paste0("'extraPeaks' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@span) > 1 | any(object@span < 0))
                 msg <- c(msg, paste0("'span' has to be a ",
                                      "positive numeric of length 1!"))
             if (length(object@smooth) > 1 |
                 !all(object@smooth %in% c("loess", "linear")))
                 msg <- c(msg, paste0("'smooth' has to be either \"",
                                      "loess\" or \"linear\"!"))
             if (length(object@family) > 1 |
                 !all(object@family %in% c("gaussian", "symmetric")))
                 msg <- c(msg, paste0("'family' has to be either \"",
                                      "gaussian\" or \"symmetric\"!"))
             if (length(msg))
                 msg
             else TRUE
         })

#' @title Align retention times across samples using Obiwarp
#'
#' @description This method performs retention time adjustment using the
#'     Obiwarp method [Prince 2006]. It is based on the code at
#'     \url{http://obi-warp.sourceforge.net} but supports alignment of multiple
#'     samples by aligning each against a \emph{center} sample. The alignment is
#'     performed directly on the \code{\link{profile-matrix}} and can hence be
#'     performed independently of the peak detection or peak grouping.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{retcor}} methods. All of the settings to the alignment
#'     algorithm can be passed with a \code{ObiwarpParam} object.
#' 
#' @param binSize \code{numeric(1)} defining the bin size (in mz dimension)
#'     to be used for the \emph{profile matrix} generation. See \code{step}
#'     parameter in \code{\link{profile-matrix}} documentation for more details.
#'
#' @param centerSample \code{integer(1)} defining the index of the center sample
#'     in the experiment. It defaults to
#'     \code{floor(median(1:length(fileNames(object))))}.
#'
#' @param response \code{numeric(1)} defining the \emph{responsiveness} of
#'     warping with \code{response = 0} giving linear warping on start and end
#'     points and \code{response = 100} warping using all bijective anchors.
#'
#' @param distFun character defining the distance function to be used. Allowed
#'     values are \code{"cor"} (Pearson's correlation), \code{"cor_opt"}
#'     (calculate only 10\% diagonal band of distance matrix; better runtime),
#'     \code{"cov"} (covariance), \code{"prd"} (product) and \code{"euc"}
#'     (Euclidian distance). The default value is \code{distFun = "cor_opt"}.
#'
#' @param gapInit \code{numeric(1)} defining the penalty for gap opening. The
#'     default value for \code{gapInit} depends on the value of \code{distFun}:
#'     for \code{distFun = "cor"} and \code{distFun = "cor_opt"} it is
#'     \code{0.3}, for \code{distFun = "cov"} and \code{distFun = "prd"}
#'     \code{0.0} and for \code{distFun = "euc"} \code{0.9}.
#'
#' @param gapExtend \code{numeric(1)} defining the penalty for gap enlargement.
#'     The default value for \code{gapExtend} depends on the value of
#'     \code{distFun}, for \code{distFun = "cor"} and
#'     \code{distFun = "cor_opt"} it is \code{2.4}, for \code{distFun = "cov"}
#'     \code{11.7}, for \code{distFun = "euc"} \code{1.8} and for
#'     \code{distFun = "prd"} {7.8}.
#'
#' @param factorDiag \code{numeric(1)} defining the local weight applied to
#'     diagonal moves in the alignment.
#'
#' @param factorGap \code{numeric(1)} defining the local weight for gap moves
#'     in the alignment.
#'
#' @param localAlignment \code{logical(1)} whether a local alignment should be
#'     performed instead of the default global alignment.
#'
#' @param initPenalty \code{numeric(1)} defining the penalty for initiating an
#'     alignment (for local alignment only).
#' 
#' @family retention time correction methods
#' 
#' @seealso \code{\link{retcor.obiwarp}} for the old user interface.
#'     \code{\link{plotAdjustedRtime}} for visualization of alignment results.
#'
#' @name adjustRtime-obiwarp
#'
#' @author Colin Smith, Johannes Rainer
#' 
#' @references
#' John T. Prince and Edward M. Marcotte. "Chromatographic Alignment of
#' ESI-LC-MS Proteomics Data Sets by Ordered Bijective Interpolated Warping"
#' \emph{Anal. Chem.} 2006, 78(17):6140-6152.

NULL
#> NULL

#' @description The \code{ObiwarpParam} class allows to specify all
#'     settings for the retention time adjustment based on the \emph{obiwarp}
#'     method. Class Instances should be created using the
#'     \code{ObiwarpParam} constructor.
#'
#' @slot .__classVersion__,binSize,centerSample,response,distFun,gapInit,gapExtend,factorDiag,factorGap,localAlignment,initPenalty See corresponding parameter above. \code{.__classVersion__} stores
#' the version from the class. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname adjustRtime-obiwarp
#'
#' @examples
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' 
#' ## Reading 2 of the KO samples
#' raw_data <- readMSData(fls[1:2], mode = "onDisk")
#'
#' ## Perform retention time correction on the OnDiskMSnExp:
#' res <- adjustRtime(raw_data, param = ObiwarpParam())
#' 
#' ## As a result we get a numeric vector with the adjusted retention times for
#' ## all spectra.
#' head(res)
#'
#' ## We can split this by file to get the adjusted retention times for each
#' ## file
#' resL <- split(res, fromFile(raw_data))
#'
#' ##############################
#' ## Perform retention time correction on an XCMSnExp:
#' ##
#' ## Perform first the chromatographic peak detection using the matchedFilter
#' ## method.
#' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
#' res <- findChromPeaks(raw_data, param = mfp)
#'
#' ## Performing the retention time adjustment using obiwarp.
#' res_2 <- adjustRtime(res, param = ObiwarpParam())
#'
#' head(rtime(res_2))
#' head(rtime(raw_data))
#'
#' ## Also the retention times of the detected peaks were adjusted.
#' tail(chromPeaks(res))
#' tail(chromPeaks(res_2))
setClass("ObiwarpParam",
         slots = c(binSize = "numeric",
                   centerSample = "integer",
                   response = "integer",
                   distFun = "character",
                   gapInit = "numeric",
                   gapExtend = "numeric",
                   factorDiag = "numeric",
                   factorGap = "numeric",
                   localAlignment = "logical",
                   initPenalty = "numeric"),
         contains = "Param",
         prototype = prototype(
             binSize = 1,
             centerSample = integer(),
             response = 1L,
             distFun = "cor_opt",
             gapInit = numeric(),
             gapExtend = numeric(),
             factorDiag = 2,
             factorGap = 1,
             localAlignment = FALSE,
             initPenalty = 0),
         validity = function(object) {
             msg <- character()
             if (length(object@binSize) > 1 |
                 any(object@binSize < 0))
                 msg <- c(msg, paste0("'binSize' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@centerSample) > 1 |
                 any(object@centerSample < 0))
                 msg <- c(msg, paste0("'centerSample' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@response) > 1 |
                 any(object@response < 0) |
                 any(object@response > 100))
                 msg <- c(msg, paste0("'response' has to be a single ",
                                      " integer from 1 to 100!"))
             if (length(object@distFun) > 1 |
                 any(!(object@distFun %in% c("cor", "cor_opt", "cov", "euc",
                                             "prd"))))
                 msg <- c(msg, paste0("'distFun' has to be one of \"cor\"",
                                      ", \"cor_opt\", \"cov\", \"euc\"",
                                      " or \"prd\"!"))
             if (length(object@gapInit) > 1 | any(object@gapInit < 0))
                 msg <- c(msg, paste0("'gapInit' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@gapExtend) > 1 | any(object@gapExtend < 0))
                 msg <- c(msg, paste0("'gapExtend' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@factorDiag) > 1 | any(object@factorDiag < 0))
                 msg <- c(msg, paste0("'factorDiag' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@factorGap) > 1 | any(object@factorGap < 0))
                 msg <- c(msg, paste0("'factorGap' has to be a positive",
                                      " numeric of length 1!"))
             if (length(object@localAlignment) > 1)
                 msg <- c(msg, paste0("'localAlignment' has to be a ",
                                      "logical of length 1!"))
             if (length(object@initPenalty) > 1 | any(object@initPenalty < 0))
                 msg <- c(msg, paste0("'initPenalty' has to be a positive",
                                      " numeric of length 1!"))
             if (length(msg))
                 msg
             else TRUE
         })

#' @description The \code{FillChromPeaksParam} object encapsules all settings for
#' the signal integration for missing peaks.
#' 
#' @slot .__classVersion__,expandMz,expandRt,ppm See corresponding parameter above. \code{.__classVersion__} stores the version of the class.
#' 
#' @rdname fillChromPeaks
setClass("FillChromPeaksParam",
         slots = c(expandMz = "numeric",
                   expandRt = "numeric",
                   ppm = "numeric"),
         contains = "Param",
         prototype = prototype(
             expandMz = 0,
             expandRt = 0,
             ppm = 0
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@expandMz) > 1 | any(object@expandMz < -1))
                 msg <- c(msg, "'expandMz' has to be > -1 and of length 1")
             if (length(object@expandRt) > 1 | any(object@expandRt < -1))
                 msg <- c(msg, "'expandRt' has to be > -1 and of length 1")
             if (length(object@ppm) > 1 | any(object@ppm < 0))
                 msg <- c(msg, paste0("'ppm' has to be a positive",
                                      " numeric of length 1!"))
             if (length(msg))
                 msg
             else TRUE
         }
         )

#' @aliases MsFeatureData
#'
#' @title Data container storing xcms preprocessing results
#'
#' @description The \code{MsFeatureData} class is designed to encapsule all
#'     data related to the preprocessing of metabolomics data using the
#'     \code{xcms} package, i.e. it contains a \code{matrix} with the
#'     chromatographic peaks identified by the peak detection, a
#'     \code{DataFrame} with the definition on grouped chromatographic peaks
#'     across samples and a \code{list} with the adjusted retention times per
#'     sample.
#'
#' @rdname XCMSnExp-class
setClass("MsFeatureData", contains = c("environment", "Versioned"),
         prototype = prototype(.xData = new.env(parent = emptyenv())))

.REQ_PEAKS_COLS <- c("mz", "mzmin", "mzmax", "rt", "rtmin",
                     "rtmax", "into", "sample")
.REQ_PEAKG_COLS <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                     "peakidx")

#' @aliases XCMSnExp
#' 
#' @title Data container storing xcms preprocessing results
#'
#' @description The \code{XCMSnExp} object is designed to contain all results
#'     from metabolomics data preprocessing (chromatographic peak detection,
#'     peak grouping (correspondence) and retention time correction). The
#'     corresponding elements in the \code{msFeatureData} slot are
#'     \code{"chromPeaks"} (a \code{matrix}), \code{"featureDefinitions"}
#'     (a \code{DataFrame}) and \code{"adjustedRtime"} (a \code{list} of
#'     numeric vectors). Note that these should not be accessed directly but
#'     rather \emph{via} their accessor methods.
#'     Along with the results, the object contains the processing history that
#'     allow to track each processing step along with the used settings. The
#'     object also directly extends the \code{\link[MSnbase]{OnDiskMSnExp}}
#'     object hence allowing easy access to the full data on which the peak
#'     detection was performed.
#'
#'     Objects from this class should not be created directly, they are
#'     returned as result from the \code{\link{findChromPeaks}} method.
#'
#'     \code{XCMSnExp} objects can be coerced into \code{\linkS4class{xcmsSet}}
#'     objects using the \code{as} method.
#'
#' @note The \code{"chromPeaks"} element in the \code{msFeatureData} slot is
#'     equivalent to the \code{@peaks} slot of the \code{xcmsSet} object, the
#'     \code{"featureDefinitions"} contains information from the \code{@groups}
#'     and \code{@groupidx} slots from an \code{xcmsSet} object.
#'
#' @slot .processHistory \code{list} with \code{XProcessHistory} objects
#'     tracking all individual analysis steps that have been performed.
#'
#' @slot msFeatureData \code{MsFeatureData} class extending \code{environment}
#'     and containing the results from a chromatographic peak detection (element
#'     \code{"chromPeaks"}), peak grouping (element \code{"featureDefinitions"})
#'     and retention time correction (element \code{"adjustedRtime"}) steps.
#'
#' @param object For \code{adjustedRtime}, \code{featureDefinitions},
#'     \code{chromPeaks}, \code{hasAdjustedRtime}, \code{hasFeatures} and
#'     \code{hasChromPeaks} either a \code{MsFeatureData} or a \code{XCMSnExp}
#'     object, for all other methods a \code{XCMSnExp} object.
#'
#' @param value For \code{adjustedRtime<-}: a \code{list} (length equal to the
#'     number of samples) with numeric vectors representing the adjusted
#'     retention times per scan.
#'
#'     For \code{featureDefinitions<-}: a \code{DataFrame} with peak
#'     grouping information. See return value for the \code{featureDefinitions}
#'     method for the expected format.
#'
#'     For \code{chromPeaks<-}: a \code{matrix} with information on
#'     detected peaks. See return value for the \code{chromPeaks} method for the
#'     expected format.
#'
#' @author Johannes Rainer
#'
#' @seealso \code{\linkS4class{xcmsSet}} for the old implementation.
#'     \code{\link[MSnbase]{OnDiskMSnExp}}, \code{\link[MSnbase]{MSnExp}}
#'     and \code{\link[MSnbase]{pSet}} for a complete list of inherited methods.
#'
#'     \code{\link{findChromPeaks}} for available peak detection methods
#'     returning a \code{XCMSnExp} object as a result.
#' 
#'     \code{\link{groupChromPeaks}} for available peak grouping
#'     methods and \code{\link{featureDefinitions}} for the method to extract
#'     the feature definitions representing the peak grouping results.
#'     \code{\link{adjustRtime}} for retention time adjustment methods.
#'
#'     \code{\link[MSnbase]{chromatogram}} to extract MS data as
#'     \code{\link[MSnbase]{Chromatogram}} objects.
#' 
#'     \code{\link{extractMsData}} for the method to extract MS data as
#'     \code{data.frame}s.
#' 
#' @rdname XCMSnExp-class
#'
#' @examples
#'
#' ## Loading the data from 2 files of the faahKO package.
#' library(faahKO)
#' od <- readMSData(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
#'                    system.file("cdf/KO/ko16.CDF", package = "faahKO")),
#'                  mode = "onDisk")
#' ## Now we perform a chromatographic peak detection on this data set using the
#' ## matched filter method. We are tuning the settings such that it performs
#' ## faster.
#' mfp <- MatchedFilterParam(binSize = 6)
#' xod <- findChromPeaks(od, param = mfp)
#'
#' ## The results from the peak detection are now stored in the XCMSnExp
#' ## object
#' xod
#'
#' ## The detected peaks can be accessed with the chromPeaks method.
#' head(chromPeaks(xod))
#'
#' ## The settings of the chromatographic peak detection can be accessed with
#' ## the processHistory method
#' processHistory(xod)
#'
#' ## Also the parameter class for the peak detection can be accessed
#' processParam(processHistory(xod)[[1]])
#'
#' ## The XCMSnExp inherits all methods from the pSet and OnDiskMSnExp classes
#' ## defined in Bioconductor's MSnbase package. To access the (raw) retention
#' ## time for each spectrum we can use the rtime method. Setting bySample = TRUE
#' ## would cause the retention times to be grouped by sample
#' head(rtime(xod))
#'
#' ## Similarly it is possible to extract the mz values or the intensity values
#' ## using the mz and intensity method, respectively, also with the option to
#' ## return the results grouped by sample instead of the default, which is
#' ## grouped by spectrum. Finally, to extract all of the data we can use the
#' ## spectra method which returns Spectrum objects containing all raw data.
#' ## Note that all these methods read the information from the original input
#' ## files and subsequently apply eventual data processing steps to them.
#' mzs <- mz(xod, bySample = TRUE)
#' length(mzs)
#' lengths(mzs)
#'
#' ## The full data could also be read using the spectra data, which returns
#' ## a list of Spectrum object containing the mz, intensity and rt values.
#' ## spctr <- spectra(xod)
#' ## To get all spectra of the first file we can split them by file
#' ## head(split(spctr, fromFile(xod))[[1]])
#'
#' ############
#' ## Filtering
#' ##
#' ## XCMSnExp objects can be filtered by file, retention time, mz values or
#' ## MS level. For some of these filter preprocessing results (mostly
#' ## retention time correction and peak grouping results) will be dropped.
#' ## Below we filter the XCMSnExp object by file to extract the results for
#' ## only the second file.
#' xod_2 <- filterFile(xod, file = 2)
#' xod_2
#'
#' ## Now the objects contains only the idenfified peaks for the second file
#' head(chromPeaks(xod_2))
#'
#' head(chromPeaks(xod)[chromPeaks(xod)[, "sample"] == 2, ])
#'
#' ##########
#' ## Coercing to an xcmsSet object
#' ##
#' ## We can also coerce the XCMSnExp object into an xcmsSet object:
#' xs <- as(xod, "xcmsSet")
#' head(peaks(xs))
setClass("XCMSnExp",
         slots = c(
             .processHistory = "list",
             msFeatureData = "MsFeatureData"
         ),
         prototype = prototype(
             .processHistory = list(),
             msFeatureData = new("MsFeatureData")
         ),
         contains = c("OnDiskMSnExp"),
         validity = function(object) {
             msg <- character()
             if (length(object@.processHistory) > 0) {
                 isOK <- unlist(lapply(object@.processHistory, function(z) {
                     return(inherits(z, "ProcessHistory"))
                 }))
                 if (!all(isOK))
                     msg <- c(msg, paste0("Only 'ProcessHistory' ",
                                          "objects are allowed in slot ",
                                          ".processHistory!"))
             }
             ## TODO @jo add checks:
             ## 1) call validMsFeatureData
             msg <- c(msg, validateMsFeatureData(object@msFeatureData))
             if (length(msg)) return(msg)
             ## 2) peaks[, "sample"] is within 1:number of samples
             if (any(ls(object@msFeatureData) == "chromPeaks")) {
                 if (!all(object@msFeatureData$chromPeaks[, "sample"] %in%
                          1:length(fileNames(object))))
                     msg <- c(msg, paste0("The number of available ",
                                          "samples does not match with ",
                                          "the sample assignment of ",
                                          "peaks in the 'chromPeaks' ",
                                          "element of the msFeatureData ",
                                          "slot!"))
             }
             ## 3) Check that the length of the adjustedRtime matches!
             if (any(ls(object@msFeatureData) == "adjustedRtime")) {
                 rt <- rtime(object, bySample = TRUE)
                 if (length(rt) != length(object@msFeatureData$adjustedRtime)) {
                     msg <- c(msg, paste0("The number of numeric vectors",
                                          " in the 'adjustedRtime' element",
                                          " of the msFeatureData slot does",
                                          " not match the number of",
                                          " samples!"))
                 } else {
                     if (any(lengths(rt) !=
                             lengths(object@msFeatureData$adjustedRtime)))
                         msg <- c(msg,
                                  paste0("The lengths of the numeric ",
                                         "vectors in the 'adjustedRtime'",
                                         " element of the msFeatureData ",
                                         "slot does not match the number",
                                         " of scans per sample!"))
                 }
             }
             ## 3) If we've got peaks, check that we have also a related
             ##    processing history step.
             if (length(msg))
                 msg
             else TRUE
         }
)

#' @aliases mz,CalibrantMassParam
#'
#' @title Calibrant mass based calibration of chromatgraphic peaks
#'
#' @description Calibrate peaks using mz values of known masses/calibrants.
#'     mz values of identified peaks are adjusted based on peaks that are close
#'     to the provided mz values. See details below for more information.
#'
#' @param mz a `numeric` or `list` of `numeric` vectors with reference mz
#'     values. If a `numeric` vector is provided, this is used for each sample
#'     in the `XCMSnExp` object. If a `list` is provided, it's length has to be
#'     equal to the number of samples in the experiment.
#'
#' @param mzabs `numeric(1)` the absolute error/deviation for matching peaks to
#'     calibrants (in Da).
#'
#' @param mzppm `numeric(1)` the relative error for matching peaks to calibrants
#'     in ppm (parts per million).
#'
#' @param neighbors `integer(1)` with the maximal number of peaks within the
#'     permitted distance to the calibrants that are considered. Among these the
#'     mz value of the peak with the largest intensity is used in the
#'     calibration function estimation.
#'
#' @param method `character(1)` defining the method that should be used to
#'     estimate the calibration function. Can be `"shift"`, `"linear"` (default)
#'     or `"edgeshift"`.
#' 
#' @details The method does first identify peaks that are close to the provided
#'     mz values and, given that there difference to the calibrants is smaller
#'     than the user provided cut off (based on arguments `mzabs` and `mzppm`),
#'     their mz values are replaced with the provided mz values. The mz values
#'     of all other peaks are either globally shifted (for `method = "shift"`
#'     or estimated by a linear model through all calibrants.
#'     Peaks are considered close to a calibrant mz if the difference between
#'     the calibrant and its mz is `<= mzabs + mz * mzppm /1e6`.
#' 
#' **Adjustment methods**: adjustment function/factor is estimated using
#' the difference between calibrant and peak mz values only for peaks
#' that are close enough to the calibrants. The availabel methods are:
#' * `shift`: shifts the m/z of each peak by a global factor which
#'   corresponds to the average difference between peak mz and calibrant mz.
#' * `linear`: fits a linear model throught the differences between
#'   calibrant and peak mz values and adjusts the mz values of all peaks
#'   using this.
#' * `edgeshift`: performs same adjustment as `linear` for peaks that are
#'   within the mz range of the calibrants and shift outside of it.
#'
#' For more information, details and examples refer to the
#' *xcms-direct-injection* vignette.
#' 
#' @note `CalibrantMassParam` classes don't have exported getter or setter
#'     methods.
#'
#' @return For `CalibrantMassParam`: a `CalibrantMassParam` instance.
#'     For `calibrate`: an [XCMSnExp] object with chromatographic peaks being
#'     calibrated. **Be aware** that the actual raw mz values are not (yet)
#'     calibrated, but **only** the identified chromatographic peaks.
#' 
#' @author Joachim Bargsten, Johannes Rainer
#'
#' @md
#' 
#' @rdname calibrate-calibrant-mass
setClass("CalibrantMassParam",
         slots = c(
             mz = "list",
             mzabs = "numeric",
             mzppm = "numeric",
             neighbors = "integer",
             method = "character"
         ),
         contains = c("Param"),
         prototype = prototype(
             mz = list(),
             mzabs = 0.0001,
             mzppm = 5,
             neighbors = 3L,
             method = "linear"
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@mz)) {
                 is_num <- vapply(object@mz, FUN = is.numeric,
                                  FUN.VALUE = logical(1), USE.NAMES = FALSE)
                 if (any(!is_num))
                     msg <- c(msg, paste0("'mz' has to be a list of numeric",
                                          " vectors"))
                 is_unsorted <- vapply(object@mz, FUN = is.unsorted,
                                       FUN.VALUE = logical(1),
                                       USE.NAMES = FALSE)
                 if (any(is_unsorted))
                     msg <- c(msg, paste0("the mz values in 'mz' have to be ",
                                          "increasingly ordered"))
             }
             if (length(object@mzppm) != 1 | any(object@mzppm < 0))
                 msg <- c(msg, paste0("'mzppm' has to be positive numeric",
                                      " of length 1."))
             if (length(object@mzabs) != 1 | any(object@mzabs < 0))
                 msg <- c(msg, paste0("'mzabs' has to be positive numeric",
                                      " of length 1."))
             if (length(object@neighbors) != 1 | any(object@neighbors <= 0))
                 msg <- c(msg, paste0("'neighbors' has to be positive integer",
                                      " of length 1."))
             if (length(object@method) != 1)
                 msg <- c(msg, paste0("'method' has to be of length 1."))
             if (!all(object@method %in% c("linear", "shift", "edgeshift")))
                 msg <- c(msg, paste0("'method' should be one of 'linear'",
                                      ", 'shift' or 'edgeshift'."))
             if (length(msg))
                 msg
             else
                 TRUE
         })
