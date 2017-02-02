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
             msg <- validMsg(NULL, NULL)
             ## Check if all slots are present.
             slNames <- slotNames(object)
             missingSlots <- character()
             for (i in 1:length(slNames)) {
                 if (!.hasSlot(object, slNames[i]))
                     missingSlots <- c(missingSlots, slNames[i])
             }
             if (length(missingSlots) > 0)
                 msg <- validMsg(msg, paste0("This xcmsSet lacks slot(s): ",
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
                     msg <- validMsg(msg,
                                     paste0("Slot '.processHistory' should",
                                            " only contain 'ProcessHistory'",
                                            " objects!"))
             }
             if (!is.null(msg))
                 return(msg)
             return(TRUE)
         }
         )

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
setClass("netCdfSource",
         ## representation(cdf="integer"),
         contains="xcmsFileSource",
         ## validity=function(object) {
         ##     if (!is.null(attr(object@cdf, "errortext"))) {
         ##         mzR:::netCDFClose(object@cdf)
         ##         attr(object@cdf, "errortext")
         ##     } else TRUE
         ## }
         )

############################################################
## rampSource
setClass("rampSource",
         ## representation(rampid="integer"),
         contains="xcmsFileSource",
         ## validity=function(object) {
         ##     if (object@rampid < 0) {
         ##         mzR:::rampClose(object@rampid)
         ##         paste("Could not open mzML/mzXML/mzData file:", object)
         ##     } else TRUE
         ## }
         )

############################################################
## xcmsPeaks
setClass("xcmsPeaks", contains = "matrix")

############################################################
## Processing history type statics
.PROCSTEP.UNKNOWN <- "Unknown"
.PROCSTEP.FEATURE.DETECTION <- "Feature detection"
.PROCSTEP.FEATURE.ALIGNMENT <- "Feature alignment"
.PROCSTEP.RTIME.CORRECTION <- "Retention time correction"
.PROCSTEPS <- c(
    .PROCSTEP.UNKNOWN,
    .PROCSTEP.FEATURE.DETECTION,
    .PROCSTEP.FEATURE.ALIGNMENT,
    .PROCSTEP.RTIME.CORRECTION
)

############################################################
## ProcessHistory
##' @aliases ProcessHistory
##' @title Tracking data processing
##'
##' @description Objects of the type \code{ProcessHistory} allow to keep track
##' of any data processing step in an metabolomics experiment. They are created
##' by the data processing methods, such as \code{\link{detectFeatures}} and
##' added to the corresponding results objects. Thus, usually, users don't need
##' to create them.
##'
##' @slot type character(1): string defining the type of the processing step.
##' This string has to match predefined values defined in the internal variable
##' \code{.PROCSTEPS}.
##'
##' @slot date character(1): date time stamp when the processing step was started.
##' @slot info character(1): optional additional information.
##' @slot fileIndex integer of length 1 or > 1 to specify on which
##' samples of the object the processing was performed.
##' @slot error (ANY): used to store eventual calculation errors.
##' @rdname ProcessHistory-class
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
##             new("Versioned", versions = c(ProcessHistory = "0.0.2"))
         ),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             ## check type:
             if (!any(object@type == .PROCSTEPS))
                 msg <- validMsg(msg, paste0("Got invalid type '", object@type,
                                             "'! Allowd are: ",
                                             paste0("\"", .PROCSTEPS, "\"",
                                                    collapse = ", ")))
             if (length(object@type) > 1)
                 msg <- validMsg(msg, paste0("length of 'type' should not be ",
                                             "larger than 1!"))
             if (length(object@date) > 1)
                 msg <- validMsg(msg, paste0("length of 'date' should not be ",
                                             "larger than 1!"))
             if (length(object@info) > 1)
                 msg <- validMsg(msg, paste0("length of 'info' should not be ",
                                             "larger than 1!"))
             if (is.null(msg)) TRUE
             else msg
         }
         )

## BasicParam class
## CentWaveParam
setClass("Param",
         representation = representation("VIRTUAL"),
         contains = c("Versioned"))
setClassUnion("ParamOrNULL", c("Param", "NULL"))

##' @aliases XProcessHistory
##' @title Tracking data processing
##'
##' @description The \code{XProcessHistory} extends the \code{ProcessHistory} by
##' adding a slot \code{param} that allows to store the actual parameter class
##' of the processing step.
##'
##' @slot param (Param): an object of type \code{Param} (e.g.
##' \code{\link{CentWaveParam}}) specifying the settings of the processing step.
##' @rdname ProcessHistory-class
setClass("XProcessHistory",
         slots = c(
             param = "ParamOrNULL"
         ),
         contains = "ProcessHistory",
         prototype = prototype(
             param = NULL
         ),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@param) > 0)
                 if(!is(object@param, "Param"))
                     msg <- validMsg(msg,
                                     paste0("Only objects from type 'Param' ",
                                            "allowed in slot '@param'! I got ",
                                            class(object@param)))
             if (is.null(msg)) TRUE
             else msg
         })



## General detectFeatures method.
##' @title Feature detection methods.
##'
##' @description The \code{detectFeature} methods are part of the modernized
##' \code{xcms} user interface.
##'
##' The implemented feature detection methods are:
##' \describe{
##' \item{centWave}{feature detection using the \emph{centWave} method.
##' See \code{\link{centWave}} for more details.}
##'
##' \item{centWave with predicted isotopes}{feature detection using a two-step
##' centWave-based approach considering also feature isotopes. See
##' \code{\link{centWaveWithPredIsoROIs}} for more details.}
##'
##' \item{matchedFilter}{peak detection in chromatographic space. See
##' \code{\link{matchedFilter}} for more details.}
##'
##' \item{massifquant}{peak detection using the Kalman filter-based feature
##' method. See \code{\link{massifquant}} for more details.}
##'
##' \item{MSW}{single-spectrum non-chromatography MS data feature detection.
##' See \code{\link{MSW}} for more details.}
##'
##' }
##' @name detectFeatures
##' @family feature detection methods
##' @seealso \code{\link{findPeaks}} for the \emph{old} feature detection
##' methods.
##' @author Johannes Rainer
NULL
#> NULL

## Main centWave documentation.
##' @title Feature detection using the centWave method
##'
##' @aliases centWave
##'
##' @description The centWave algorithm perform peak density and wavelet based
##' feature detection for high resolution LC/MS data in centroid
##' mode [Tautenhahn 2008].
##'
##' @param ppm Maximal tolerated m/z deviation in consecutive scans in parts
##' per million (ppm).
##' @param peakwidth numeric(2) with the expected approximate
##' feature/peak width in chromatographic space. Given as a range (min, max)
##' in seconds.
##' @param snthresh numeric(1) defining the signal to noise ratio cutoff.
##' @param prefilter numeric(2): \code{c(k, I)} specifying the prefilter
##' step for the first analysis step (ROI detection). Mass traces are only
##' retained if they contain at least \code{k} peaks with intensity \code{>= I}.
##' @param mzCenterFun Name of the function to calculate the m/z center of the
##' feature. Allowed are: \code{"wMean"}: intensity weighted mean of the feature's
##' m/z values, \code{"mean"}: mean of the feature's m/z values, \code{"apex"}:
##' use the m/z value at the peak apex, \code{"wMeanApex3"}: intensity weighted
##' mean of the m/z value at the peak apex and the m/z values left and right of
##' it and \code{"meanApex3"}: mean of the m/z value of the peak apex and the
##' m/z values left and right of it.
##' @param integrate Integration method. For \code{integrate = 1} peak limits
##' are found through descent on the mexican hat filtered data, for
##' \code{integrate = 2} the descent is done on the real data. The latter method
##' is more accurate but prone to noise, while the former is more robust, but
##' less exact.
##' @param mzdiff Numeric representing the minimum difference in m/z dimension
##' for peaks with overlapping retention times; can be negatove to allow overlap.
##' @param fitgauss Logical whether or not a Gaussian should be fitted to each
##' peak.
##' @param noise numeric(1) allowing to set a minimum intensity required
##' for centroids to be considered in the first analysis step (centroids with
##' intensity \code{< noise} are omitted from ROI detection).
##' @param verboseColumns Logical whether additional feature meta data columns
##' should be returned.
##' @param roiList An optional list of regions-of-interest (ROI) representing
##' detected mass traces. If ROIs are submitted the first analysis step is
##' omitted and feature detection is performed on the submitted ROIs. Each
##' ROI is expected to have the following elements specified:
##' \code{scmin} (start scan index), \code{scmax} (end scan index),
##' \code{mzmin} (minimum m/z), \code{mzmax} (maximum m/z), \code{length}
##' (number of scans), \code{intensity} (summed intensity). Each ROI should be
##' represented by a \code{list} of elements or a single row \code{data.frame}.
##' @param firstBaselineCheck logical(1). If \code{TRUE} continuous
##' data within regions of interest is checked to be above the first baseline.
##' @param roiScales Optional numeric vector with length equal to \code{roiList}
##' defining the scale for each region of interest in \code{roiList} that should
##' be used for the centWave-wavelets.
##'
##' @details The centWave algorithm is most suitable for high resolution
##' LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase the
##' method identifies \emph{regions of interest} (ROIs) representing mass traces
##' that are characterized as regions with less than \code{ppm} m/z deviation in
##' consecutive scans in the LC/MS map. These ROIs are then subsequently
##' analyzed using continuous wavelet transform (CWT) to locate chromatographic
##' peaks on different scales. The first analysis step is skipped, if regions
##' of interest are passed \emph{via} the \code{param} parameter.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{findPeaks}} methods. It supports feature detection on
##' \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
##' objects (both defined in the \code{MSnbase} package). All of the settings
##' to the centWave algorithm can be passed with a \code{CentWaveParam} object.
##'
##' @family feature detection methods
##' @seealso The \code{\link{do_detectFeatures_centWave}} core API function and
##' \code{\link{findPeaks.centWave}} for the old user interface.
##'
##' @references
##' Ralf Tautenhahn, Christoph B\"{o}ttcher, and Steffen Neumann "Highly
##' sensitive feature detection for high resolution LC/MS" \emph{BMC Bioinformatics}
##' 2008, 9:504
##' @name featureDetection-centWave
##' @author Ralf Tautenhahn, Johannes Rainer
NULL
#> NULL

##' @description The \code{CentWaveParam} class allows to specify all settings for
##' a feature detection using the centWave method. Instances should be created
##' with the \code{CentWaveParam} constructor.
##'
##' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname featureDetection-centWave
##'
##' @examples
##'
##' ## Create a CentWaveParam object. Note that the noise is set to 10000 to
##' ## speed up the execution of the example - in a real use case the default
##' ## value should be used, or it should be set to a reasonable value.
##' cwp <- CentWaveParam(ppm = 20, noise = 10000)
##' ## Change snthresh parameter
##' snthresh(cwp) <- 25
##' cwp
##'
##' ## Perform the feature detection using centWave on some of the files from the
##' ## faahKO package. Files are read using the readMSData2 from the MSnbase
##' ## package
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' raw_data <- readMSData2(fls[1:2])
##'
##' ## Perform the feature detection using the settings defined above. We're
##' ## returning the results as an xcmsSet object.
##' res <- detectFeatures(raw_data, param = cwp, return.type = "xcmsSet")
##' head(peaks(res))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@ppm) != 1 | any(object@ppm < 0))
                 msg <- validMsg(msg, paste0("'ppm' has to be positive numeric",
                                             " of length 1."))
             if (length(object@peakwidth) != 2 | any(object@peakwidth < 0))
                 msg <- validMsg(msg, paste0("'peakwidth' has to be a numeric",
                                             " of length 2 with only positive",
                                             " values."))
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- validMsg(msg, paste0("'snthresh' has to be a positive",
                                             " numeric of length 1."))
             if (length(object@prefilter) != 2)
                 msg <- validMsg(msg, paste0("'prefilter' has to be a numeric",
                                             " of length 2."))
             allowed_vals <- c("wMean", "mean", "apex", "wMeanApex3",
                               "meanApex3")
             if (!(object@mzCenterFun) %in% allowed_vals)
                 msg <- validMsg(msg, paste0("'mzCenterFun' has to be one of ",
                                             paste0("'", allowed_vals, "'",
                                             collapse = ", "), "."))
             if (!(object@integrate %in% c(1L, 2L)))
                 msg <- validMsg(msg, paste0("'integrate' has to be either 1",
                                             " or 2."))
             if (length(object@mzdiff) != 1)
                 msg <- validMsg(msg, paste0("'mzdiff' has to be a numeric of",
                                             " length 1."))
             if (length(object@noise) != 1)
                 msg <- validMsg(msg, paste0("'noise' has to be a numeric of",
                                             " length 1."))
             if (length(object@fitgauss) != 1)
                 msg <- validMsg(msg, paste0("'fitgauss' has to be a numeric of",
                                             " length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- validMsg(msg, paste0("'verboseColumns' has to be a ",
                                             "numeric of length 1."))
             if (length(object@firstBaselineCheck) != 1)
                 msg <- validMsg(msg, paste0("'firstBaselineCheck' has to be a",
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
                     msg <- validMsg(msg, paste0("'roiList' does not provide ",
                                                 "all required fields!"))
             }
             if (length(object@roiList) > 0 &
                 length(object@roiList) != length(object@roiScales))
                 msg <- validMsg(msg, paste0("'roiScales' has to have the same",
                                             " length than 'roiList'."))
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })

## Main matchedFilter documentation.
##' @title Peak detection in the chromatographic time domain
##'
##' @aliases matchedFilter
##'
##' @description The \emph{matchedFilter} algorithm identifies features in the
##' chromatographic time domain as described in [Smith 2006]. The intensity
##' values are binned by cutting The LC/MS data into slices (bins) of a mass unit
##' (\code{binSize} m/z) wide. Within each bin the maximal intensity is selected.
##' The feature detection is then performed in each bin by extending it based on
##' the \code{steps} parameter to generate slices comprising bins
##' \code{current_bin - steps +1} to \code{current_bin + steps - 1}. Each of
##' these slices is then filtered with matched filtration using a second-derative
##' Gaussian as the model feature/peak shape. After filtration features are
##' detected using a signal-to-ration cut-off. For more details and
##' illustrations see [Smith 2006].
##'
##' @param binSize numeric(1) specifying the width of the
##' bins/slices in m/z dimension.
##' @param impute Character string specifying the method to be used for missing
##' value imputation. Allowed values are \code{"none"} (no linear interpolation),
##' \code{"lin"} (linear interpolation), \code{"linbase"} (linear interpolation
##' within a certain bin-neighborhood) and \code{"intlin"}. See
##' \code{\link{imputeLinInterpol}} for more details.
##' @param fwhm numeric(1) specifying the full width at half maximum
##' of matched filtration gaussian model peak. Only used to calculate the actual
##' sigma, see below.
##' @param sigma numeric(1) specifying the standard deviation (width)
##' of the matched filtration model peak.
##' @param max numeric(1) representing the maximum number of peaks
##' that are expected/will be identified per slice.
##' @param snthresh numeric(1) defining the signal to noise cutoff
##' to be used in the feature detection step.
##' @param steps numeric(1) defining the number of bins to be
##' merged before filtration (i.e. the number of neighboring bins that will be
##' joined to the slice in which filtration and peak detection will be
##' performed).
##' @param mzdiff numeric(1) defining the minimum difference
##' in m/z for peaks with overlapping retention times
##' @param index Logical specifying whether indicies should be returned instead
##' of values for m/z and retention times.
##'
##' @details The intensities are binned by the provided m/z values within each
##' spectrum (scan). Binning is performed such that the bins are centered around
##' the m/z values (i.e. the first bin includes all m/z values between
##' \code{min(mz) - bin_size/2} and \code{min(mz) + bin_size/2}).
##'
##' For more details on binning and missing value imputation see
##' \code{\link{binYonX}} and \code{\link{imputeLinInterpol}} methods.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{findPeaks}} methods. It supports feature detection on
##' \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
##' objects (both defined in the \code{MSnbase} package). All of the settings
##' to the matchedFilter algorithm can be passed with a
##' \code{MatchedFilterParam} object.
##'
##' @inheritParams imputeLinInterpol
##' @inheritParams featureDetection-centWave
##'
##' @family feature detection methods
##' @seealso The \code{\link{do_detectFeatures_matchedFilter}} core API function
##' and \code{\link{findPeaks.matchedFilter}} for the old user interface.
##'
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
##' @author Colin A Smith, Johannes Rainer
##'
##' @name featureDetection-matchedFilter
NULL
#> NULL

##' @description The \code{MatchedFilterParam} class allows to specify all
##' settings for a feature detection using the matchedFilter method. Instances
##' should be created with the \code{MatchedFilterParam} constructor.
##'
##' @slot .__classVersion__,binSize,impute,baseValue,distance,fwhm,sigma,max,snthresh,steps,mzdiff,index See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname featureDetection-matchedFilter
##'
##' @examples
##'
##' ## Create a MatchedFilterParam object
##' mfp <- MatchedFilterParam(binSize = 0.5)
##' ## Change snthresh parameter
##' snthresh(mfp) <- 15
##' mfp
##'
##' ## Perform the feature detection using matchecFilter on the files from the
##' ## faahKO package. Files are read using the readMSData2 from the MSnbase
##' ## package
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' raw_data <- readMSData2(fls)
##' ## Perform the feature detection using the settings defined above. We're
##' ## returning the results as an xcmsSet object. Note that we are also
##' ## disabling parallel processing in this example by registering a "SerialParam"
##' register(SerialParam())
##' res <- detectFeatures(raw_data, param = mfp, return.type = "xcmsSet")
##' head(peaks(res))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@binSize) != 1 | any(object@binSize < 0))
                 msg <- validMsg(msg, paste0("'binSize' has to be positive",
                                             " numeric of length 1."))
             if (!any(c("none", "lin", "linbase") == object@impute))
                 msg <- validMsg(msg,
                                 paste0("Only values 'none', 'lin' and ",
                                        "'linbase' are allowed for'impute'"))
             if (length(object@baseValue) > 1)
                 msg <- validMsg(msg, paste0("'baseValue' has to be a",
                                             " numeric of length 1."))
             if (length(object@distance) > 1)
                 msg <- validMsg(msg, paste0("'distance' has to be a numeric",
                                             " of length 1."))
             if (length(object@fwhm) != 1)
                 msg <- validMsg(msg, paste0("'fwhm' has to be a numeric",
                                             " of length 1."))
             if (length(object@sigma) != 1)
                 msg <- validMsg(msg, paste0("'sigma' has to be a numeric",
                                             " of length 1."))
             if (length(object@max) != 1)
                 msg <- validMsg(msg, paste0("'max' has to be a numeric",
                                             " of length 1."))
             if (length(object@snthresh) != 1)
                 msg <- validMsg(msg, paste0("'snthresh' has to be a numeric",
                                             " of length 1."))
             if (length(object@steps) != 1)
                 msg <- validMsg(msg, paste0("'steps' has to be a numeric",
                                             " of length 1."))
             if (length(object@mzdiff) != 1)
                 msg <- validMsg(msg, paste0("'mzdiff' has to be a numeric",
                                             " of length 1."))
             if (length(object@index) != 1)
                 msg <- validMsg(msg, paste0("'index' has to be a logical",
                                             " of length 1."))
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })


## Main massifquant documentation.
##' @title Feature detection using the massifquant method
##'
##' @aliases massifquant
##'
##' @description Massifquant is a Kalman filter (KF)-based feature
##' detection for XC-MS data in centroid mode. The identified features
##' can be further refined with the \emph{centWave} method (see
##' \code{\link{do_detectFeatures_centWave}} for details on centWave)
##' by specifying \code{withWave = TRUE}.
##'
##' @param peakwidth numeric(2). Only the first element is used by
##' massifquant, which specifices the minimum feature length in time scans.
##' For \code{withWave = TRUE} the second argument represents the maximum
##' feature length subject to being greater than the mininum feature length
##' (see also documentation of \code{\link{do_detectFeatures_centWave}}).
##' @param prefilter numeric(2). The first argument is only used
##' if (\code{withWave = TRUE}); see \code{\link{do_detectFeatures_centWave}}
##' for details. The second argument specifies the minimum threshold for the
##' maximum intensity of a feature that must be met.
##' @param criticalValue numeric(1). Suggested values:
##' (\code{0.1-3.0}). This setting helps determine the the Kalman Filter
##' prediciton margin of error. A real centroid belonging to a bonafide
##' feature must fall within the KF prediction margin of error. Much like
##' in the construction of a confidence interval, \code{criticalVal} loosely
##' translates to be a multiplier of the standard error of the prediction
##' reported by the Kalman Filter. If the features in the XC-MS sample have
##' a small mass deviance in ppm error, a smaller critical value might be
##' better and vice versa.
##' @param consecMissedLimit Integer: Suggested values: (\code{1,2,3}). While
##' a feature is in the proces of being detected by a Kalman Filter, the
##' Kalman Filter may not find a predicted centroid in every scan. After 1
##' or more consecutive failed predictions, this setting informs Massifquant
##' when to stop a Kalman Filter from following a candidate feature.
##' @param unions Integer: set to \code{1} if apply t-test union on
##' segmentation; set to \code{0} if no t-test to be applied on
##' chromatographically continous features sharing same m/z range.
##' Explanation: With very few data points, sometimes a Kalman Filter stops
##' tracking a feature prematurely. Another Kalman Filter is instantiated
##' and begins following the rest of the signal. Because tracking is done
##' backwards to forwards, this algorithmic defect leaves a real feature
##' divided into two segments or more. With this option turned on, the
##' program identifies segmented features and combines them (merges them)
##' into one with a two sample t-test. The potential danger of this option
##' is that some truly distinct features may be merged.
##' @param checkBack Integer: set to \code{1} if turned on; set to \code{0}
##' if turned off. The convergence of a Kalman Filter to a feature's precise
##' m/z mapping is very fast, but sometimes it incorporates erroneous centroids
##' as part of a feature (especially early on). The \code{scanBack} option is an
##' attempt to remove the occasional outlier that lies beyond the converged
##' bounds of the Kalman Filter. The option does not directly affect
##' identification of a feature because it is a postprocessing measure; it
##' has not shown to be a extremely useful thus far and the default is set
##' to being turned off.
##' @param withWave Logical: if \code{TRUE}, the features identified first
##' with Massifquant are subsequently filtered with the second step of the
##' centWave algorithm, which includes wavelet estimation.
##'
##' @details This algorithm's performance has been tested rigorously
##' on high resolution LC/{OrbiTrap, TOF}-MS data in centroid mode.
##' Simultaneous kalman filters identify features and calculate their
##' area under the curve. The default parameters are set to operate on
##' a complex LC-MS Orbitrap sample. Users will find it useful to do some
##' simple exploratory data analysis to find out where to set a minimum
##' intensity, and identify how many scans an average feature spans. The
##' \code{consecMissedLimit} parameter has yielded good performance on
##' Orbitrap data when set to (\code{2}) and on TOF data it was found best
##' to be at (\code{1}). This may change as the algorithm has yet to be
##' tested on many samples. The \code{criticalValue} parameter is perhaps
##' most dificult to dial in appropriately and visual inspection of peak
##' identification is the best suggested tool for quick optimization.
##' The \code{ppm} and \code{checkBack} parameters have shown less influence
##' than the other parameters and exist to give users flexibility and
##' better accuracy.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{findPeaks}} methods. It supports feature detection on
##' \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
##' objects (both defined in the \code{MSnbase} package). All of the settings
##' to the massifquant and centWave algorithm can be passed with a
##' \code{MassifquantParam} object.
##'
##' @inheritParams featureDetection-centWave
##'
##' @family feature detection methods
##' @seealso The \code{\link{do_detectFeatures_massifquant}} core API function
##' and \code{\link{findPeaks.massifquant}} for the old user interface.
##'
##' @references
##' Conley CJ, Smith R, Torgrip RJ, Taylor RM, Tautenhahn R and Prince JT
##' "Massifquant: open-source Kalman filter-based XC-MS isotope trace feature
##' detection" \emph{Bioinformatics} 2014, 30(18):2636-43.
##' @author Christopher Conley, Johannes Rainer
##'
##' @name featureDetection-massifquant
NULL
#> NULL

##' @description The \code{MassifquantParam} class allows to specify all
##' settings for a feature detection using the massifquant method eventually in
##' combination with the centWave algorithm. Instances should be created with
##' the \code{MassifquantParam} constructor.
##'
##' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,criticalValue,consecMissedLimit,unions,checkBack,withWave See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname featureDetection-massifquant
##'
##' @examples
##'
##' ## Create a MassifquantParam object.
##' mqp <- MassifquantParam()
##' ## Change snthresh parameter
##' snthresh(mqp) <- 30
##' mqp
##'
##' ## Perform the feature detection using massifquant on the files from the
##' ## faahKO package. Files are read using the readMSData2 from the MSnbase
##' ## package
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' raw_data <- readMSData2(fls[1:2])
##' ## Perform the feature detection using the settings defined above. We're
##' ## returning the results as an xcmsSet object.
##' res <- detectFeatures(raw_data, param = mqp, return.type = "xcmsSet")
##' head(peaks(res))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@ppm) != 1 | any(object@ppm < 0))
                 msg <- validMsg(msg, paste0("'ppm' has to be positive numeric",
                                             " of length 1."))
             if (length(object@peakwidth) != 2 | any(object@peakwidth < 0))
                 msg <- validMsg(msg, paste0("'peakwidth' has to be a numeric",
                                             " of length 2 with only positive",
                                             " values."))
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- validMsg(msg, paste0("'snthresh' has to be a positive",
                                             " numeric of length 1."))
             if (length(object@prefilter) != 2)
                 msg <- validMsg(msg, paste0("'prefilter' has to be a numeric",
                                             " of length 2."))
             allowed_vals <- c("wMean", "mean", "apex", "wMeanApex3",
                               "meanApex3")
             if (!(object@mzCenterFun) %in% allowed_vals)
                 msg <- validMsg(msg, paste0("'mzCenterFun' has to be one of ",
                                             paste0("'", allowed_vals, "'",
                                             collapse = ", "), "."))
             if (!(object@integrate %in% c(1L, 2L)))
                 msg <- validMsg(msg, paste0("'integrate' has to be either 1",
                                             " or 2."))
             if (length(object@mzdiff) != 1)
                 msg <- validMsg(msg, paste0("'mzdiff' has to be a numeric of",
                                             " length 1."))
             if (length(object@noise) != 1)
                 msg <- validMsg(msg, paste0("'noise' has to be a numeric of",
                                             " length 1."))
             if (length(object@fitgauss) != 1)
                 msg <- validMsg(msg, paste0("'fitgauss' has to be a numeric of",
                                             " length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- validMsg(msg, paste0("'verboseColumns' has to be a ",
                                             "numeric of length 1."))
             if (length(object@criticalValue) != 1)
                 msg <- validMsg(msg, paste0("'criticalValue' has to be a ",
                                             "numeric of length 1."))
             if (length(object@consecMissedLimit) != 1)
                 msg <- validMsg(msg, paste0("'consecMissedLimit' has to be a ",
                                             "numeric of length 1."))
             if (length(object@unions) != 1)
                 msg <- validMsg(msg, paste0("'unions' has to be a ",
                                             "numeric of length 1."))
             if (object@unions != 0 & object@unions != 1)
                 msg <- validMsg(msg, paste0("'unions' has to be either 0 or 1!"))
             if (length(object@checkBack) != 1)
                 msg <- validMsg(msg, paste0("'checkBack' has to be a ",
                                             "numeric of length 1."))
             if (object@checkBack != 0 & object@checkBack != 1)
                 msg <- validMsg(msg, paste0("'checkBack' has to be either 0",
                                             " or 1!"))
             if (length(object@withWave) != 1)
                 msg <- validMsg(msg, paste0("'withWave' has to be a ",
                                             "numeric of length 1."))
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })

## Main MSW documentation.
##' @title Single-spectrum non-chromatography MS data feature detection
##'
##' @aliases MSW
##'
##' @description Perform feature detection in mass spectrometry
##' direct injection spectrum using a wavelet based algorithm.
##'
##' @details This is a wrapper for the peak picker in Bioconductor's
##' \code{MassSpecWavelet} package calling
##' \code{\link[MassSpecWavelet]{peakDetectionCWT}} and
##' \code{\link[MassSpecWavelet]{tuneInPeakInfo}} functions. See the
##' \emph{xcmsDirect} vignette for more information.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{findPeaks}} methods. It supports feature detection on
##' \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
##' objects (both defined in the \code{MSnbase} package). All of the settings
##' to the massifquant and centWave algorithm can be passed with a
##' \code{MassifquantParam} object.
##'
##' @inheritParams featureDetection-centWave
##'
##' @family feature detection methods
##' @seealso The \code{\link{do_detectFeatures_MSW}} core API function
##' and \code{\link{findPeaks.MSW}} for the old user interface.
##'
##' @author Joachim Kutzera, Steffen Neumann, Johannes Rainer
##'
##' @name featureDetection-MSW
NULL
#> NULL

##' @description The \code{MSWParam} class allows to specify all
##' settings for a feature detection using the MSW method. Instances should be
##' created with the \code{MSWParam} constructor.
##'
##' @slot .__classVersion__,snthresh,verboseColumns,scales,nearbyPeak,peakScaleRange,ampTh,minNoiseLevel,ridgeLength,peakThr,tuneIn,addParams See corresponding parameter above. \code{.__classVersion__} stores the version from the class. Slots values
##' should exclusively be accessed \emph{via} the corresponding getter and
##' setter methods listed above.
##'
##' @rdname featureDetection-MSW
##'
##' @examples
##'
##' ## Create a MSWParam object
##' mp <- MSWParam()
##' ## Change snthresh parameter
##' snthresh(mp) <- 15
##' mp
##'
##' ## Loading a small subset of direct injection, single spectrum files
##' library(msdata)
##' fticrf <- list.files(system.file("fticr", package = "msdata"),
##'                     recursive = TRUE, full.names = TRUE)
##' fticr <- readMSData2(fticrf[1:2], msLevel. = 1)
##'
##' ## Perform the MSW feature detection on these:
##' p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
##'              SNR.method = "data.mean", winSize.noise = 500)
##' fticr <- detectFeatures(fticr, param = p)
##'
##' head(features(fticr))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@snthresh) != 1 | any(object@snthresh < 0))
                 msg <- validMsg(msg, paste0("'snthresh' has to be a positive",
                                             " numeric of length 1."))
             if (length(object@verboseColumns) != 1)
                 msg <- validMsg(msg, paste0("'verboseColumns' has to be a ",
                                             "numeric of length 1."))
             if (length(object@nearbyPeak) != 1)
                 msg <- validMsg(msg, paste0("'nearbyPeak' has to be a ",
                                             "logical of length 1."))
             if (length(object@peakScaleRange) != 1 |
                 any(object@peakScaleRange < 0))
                 msg <- validMsg(msg, paste0("'peakScaleRange' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@ampTh) != 1 | any(object@ampTh < 0))
                 msg <- validMsg(msg, paste0("'ampTh' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@minNoiseLevel) != 1 |
                 any(object@minNoiseLevel < 0))
                 msg <- validMsg(msg, paste0("'minNoiseLevel' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@ridgeLength) != 1 |
                 any(object@ridgeLength < 0))
                 msg <- validMsg(msg, paste0("'ridgeLength' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@peakThr) > 1)
                 msg <- validMsg(msg, paste0("'peakThr' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@tuneIn) != 1)
                 msg <- validMsg(msg, paste0("'tuneIn' has to be a ",
                                             "logical of length 1."))
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })

##' @title Two-step centWave feature detection considering also feature isotopes
##'
##' @aliases centWaveWithPredIsoROIs
##'
##' @description This method performs a two-step centWave-based feature
##' detection: in a first centWave run features are identified for which then
##' the location of their potential isotopes in the mz-retention time is
##' predicted. A second centWave run is then performed on these
##' \emph{regions of interest} (ROIs). The final list of features comprises all
##' non-overlapping features from both centWave runs.
##'
##' @inheritParams featureDetection-centWave
##'
##' @param maxCharge integer(1) defining the maximal isotope charge. Isotopes
##' will be defined for charges \code{1:maxCharge}.
##'
##' @param maxIso integer(1) defining the number of isotope peaks that should be
##' predicted for each feature identified in the first centWave run.
##'
##' @param mzIntervalExtension logical(1) whether the mz range for the predicted
##' isotope ROIs should be extended to increase detection of low intensity peaks.
##'
##' @param snthreshIsoROIs numeric(1) defining the signal to noise ratio cutoff
##' to be used in the second centWave run to identify features for predicted
##' isotope ROIs.
##'
##' @param polarity character(1) specifying the polarity of the data. Currently
##' not used, but has to be \code{"positive"}, \code{"negative"} or
##' \code{"unknown"} if provided.
##'
##' @details See \code{\link{centWave}} for details on the centWave method.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{findPeaks}} methods. It supports feature detection on
##' \code{\link[MSnbase]{MSnExp}} and \code{\link[MSnbase]{OnDiskMSnExp}}
##' objects (both defined in the \code{MSnbase} package). All of the settings
##' to the centWave algorithm can be passed with a \code{CentWaveParam} object.
##'
##' @family feature detection methods
##' @seealso The \code{\link{do_detectFeatures_centWaveWithPredIsoROIs}} core
##' API function and \code{\link{findPeaks.centWave}} for the old user interface.
##' \code{\link{CentWaveParam}} for the class the \code{CentWavePredIsoParam}
##' extends.
##'
##' @name featureDetection-centWaveWithPredIsoROIs
##' @author Hendrik Treutler, Johannes Rainer
NULL
#> NULL

##' @description The \code{CentWavePredIsoParam} class allows to specify all
##' settings for the two-step centWave-based feature detection considering also
##' predicted isotopes of features identified in the first centWave run.
##' Instances should be created with the \code{CentWavePredIsoParam} constructor.
##' See also the documentation of the \code{\link{CentWaveParam}} for all methods
##' and arguments this class inherits.
##'
##' @slot .__classVersion__,ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales,snthreshIsoROIs,maxCharge,maxIso,mzIntervalExtension,polarity See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname featureDetection-centWaveWithPredIsoROIs
##'
##' @examples
##'
##' ## Create a param object
##' p <- CentWavePredIsoParam(maxCharge = 4)
##' ## Change snthresh parameter
##' snthresh(p) <- 25
##' p
##'
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
             msg <- validMsg(NULL, NULL)
             if (length(object@snthreshIsoROIs) != 1 |
                 any(object@snthreshIsoROIs < 0))
                 msg <- validMsg(msg, paste0("'snthreshIsoROIs' has to be a ",
                                             "positive numeric of length 1."))
             if (length(object@maxCharge) != 1 | any(object@maxCharge < 0))
                 msg <- validMsg(msg, paste0("'maxCharge' has to be a ",
                                             "positive integer of length 1."))
             if (length(object@maxIso) != 1 | any(object@maxIso < 0))
                 msg <- validMsg(msg, paste0("'maxIso' has to be a ",
                                             "positive integer of length 1."))
             if (length(object@mzIntervalExtension) != 1)
                 msg <- validMsg(msg, paste0("'mzIntervalExtension' has to be a",
                                             " logical of length 1."))
             if (length(object@polarity) != 1)
                 msg <- validMsg(msg, paste0("'polarity' has to be a",
                                             " character of length 1."))
             if (!(object@polarity %in% c("positive", "negative", "unknown")))
                 msg <- validMsg(msg, paste0("'polarity' has to be either ",
                                             "'positive', 'negative' or ",
                                             "'unknown'!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })


## General groupFeatures method.
##' @title Feature alignment methods.
##'
##' @description The \code{groupFeatures} method(s) perform alignment of features
##' within and between samples. These methods are part of the modernized
##' \code{xcms} user interface.
##'
##' The implemented feature alignment methods are:
##' \describe{
##' \item{density}{feature alignment based on time dimension feature densities.
##' See \code{\link{groupFeatures-density}} for more details.}
##'
##' \item{mzClust}{high resolution feature alignment for single spectra (direct
##' infusion) MS data. See \code{\link{groupFeatures-mzClust}} for more details.}
##'
##' \item{nearest}{feature alignment based on their proximity in the mz-rt space.
##' See \code{\link{groupFeatures-nearest}} for more details.}
##' 
##' }
##' @name groupFeatures
##' @family feature alignment methods
##' @seealso \code{\link{group}} for the \emph{old} feature alignment
##' methods.
##' @seealso \code{\link{featureGroups}} and \code{\link{groupval,XCMSnExp-method}}
##' for methods to access feature grouping results.
##' 
##' @author Johannes Rainer
NULL
#> NULL

## Main group.density documentation.
##' @title Feature alignment based on time dimension feature densities
##'
##' @description This method performs performs feature alignment based on the
##' density (distribution) of identified features along the retention time axis
##' within slices of overlapping mz ranges. All features (from the same or from
##' different samples) being close on the retention time axis are grouped into
##' a \emph{feature group}.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{group}} methods. All of the settings to the alignment algorithm
##' can be passed with a \code{FeatureDensityParam} object.
##'
##' @param sampleGroups A vector of the same length than samples defining the
##' sample group assignments.
##'
##' @param bw numeric(1) defining the bandwidth (standard deviation ot the
##' smoothing kernel) to be used. This argument is passed to the
##' \code{\link{density}} method.
##'
##' @param minFraction numeric(1) defining the minimum fraction of samples in at
##' least one sample group in which the features have to be present to be
##' considered as a feature group.
##'
##' @param minSamples numeric(1) with the minimum number of samples in at least
##' one sample group in which the features have to be detected to be considered
##' as a feature group.
##'
##' @param binSize numeric(1) defining the size of the overlapping slices in mz
##' dimension.
##'
##' @param maxFeatures numeric(1) with the maximum number of feature groups to
##' be identified in a single mz slice.
##' 
##' @family feature alignment methods
##' @seealso The \code{\link{do_groupFeatures_density}} core
##' API function and \code{\link{group.density}} for the old user interface.
##' @seealso \code{\link{featureGroups}} and \code{\link{groupval,XCMSnExp-method}}
##' for methods to access feature grouping results.
##'
##' @name groupFeatures-density
##' 
##' @author Colin Smith, Johannes Rainer
##'
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
NULL
#> NULL

##' @description The \code{FeatureDensityParam} class allows to specify all
##' settings for the feature alignment based on feature densities along the time
##' dimension. Instances should be created with the \code{FeatureDensityParam}
##' constructor.
##'
##' @slot .__classVersion__,sampleGroups,bw,minFraction,minSamples,binSize,maxFeatures See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname groupFeatures-density
##'
##' @examples
##'
##' ## Create a FeatureDensityParam object
##' p <- FeatureDensityParam(binSize = 0.05)
##' ## Change hte minSamples slot
##' minSamples(p) <- 3
##' p
##'
##' ##############################
##' ## feature detection and alignment.
##' ##
##' ## Below we perform first a feature detection (using the matchedFilter
##' ## method) on some of the test files from the faahKO package followed by
##' ## a feature alignment using the density method.
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' 
##' ## Reading 2 of the KO samples
##' raw_data <- readMSData2(fls[1:2])
##'
##' ## Perform the feature detection using the matchedFilter method.
##' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
##' res <- detectFeatures(raw_data, param = mfp)
##'
##' head(features(res))
##' ## The number of features identified per sample:
##' table(features(res)[, "sample"])
##'
##' ## Performing the feature alignment
##' fdp <- FeatureDensityParam()
##' res <- groupFeatures(res, fdp)
##'
##' ## The results from the feature alignment:
##' featureGroups(res)
##'
##' ## Using the groupval method to extract a matrix with the intensities of
##' ## the feature groups per sample.
##' head(groupval(res, value = "into"))
##' 
##' ## The process history:
##' processHistory(res)
setClass("FeatureDensityParam",
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
             msg <- validMsg(NULL, NULL)
             if (length(object@bw) > 1 | any(object@bw < 0))
                 msg <- validMsg(msg, paste0("'bw' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@minFraction) > 1 | any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- validMsg(msg, paste0("'minFraction' has to be a ",
                                             "single positive number between ",
                                             "0 and 1!"))
             if (length(object@minSamples) > 1 | any(object@minSamples < 0))
                 msg <- validMsg(msg, paste0("'minSamples' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@binSize) > 1 | any(object@binSize < 0))
                 msg <- validMsg(msg, paste0("'binSize' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@maxFeatures) > 1 | any(object@maxFeatures < 0))
                 msg <- validMsg(msg, paste0("'maxFeatures' has to be a ",
                                             "positive numeric of length 1!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })

## Main group.mzClust documentation.
##' @title High resolution feature alignment for single spectra samples
##'
##' @description This method performs high resolution alignment for single
##' spectra samples.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{group}} methods. All of the settings to the alignment algorithm
##' can be passed with a \code{MzClustParam} object.
##'
##' @inheritParams groupFeatures-density
##'
##' @param ppm numeric(1) representing the relative mz error for the
##' clustering/grouping (in parts per million).
##' 
##' @param absMz numeric(1) representing the absolute mz error for the clustering.
##' 
##' @family feature alignment methods
##' @seealso The \code{\link{do_groupFeatures_mzClust}} core
##' API function and \code{\link{group.mzClust}} for the old user interface.
##' @seealso \code{\link{featureGroups}} and \code{\link{groupval,XCMSnExp-method}}
##' for methods to access feature grouping results.
##'
##' @name groupFeatures-mzClust
##'
##' @references Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill
##' and David F. Grant\cr \emph{Alignment of high resolution mass spectra:
##' development of a heuristic approach for metabolomics}.\cr Metabolomics,
##' Vol. 2, No. 2, 75-83 (2006)
NULL
#> NULL

##' @description The \code{MzClustParam} class allows to specify all
##' settings for the feature alignment based on the \emph{mzClust} algorithm.
##' Instances should be created with the \code{MzClustParam} constructor.
##'
##' @slot .__classVersion__,sampleGroups,ppm,absMz,minFraction,minSamples See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname groupFeatures-mzClust
##'
##' @examples
##'
##' ## Loading a small subset of direct injection, single spectrum files
##' library(msdata)
##' fticrf <- list.files(system.file("fticr", package = "msdata"),
##'                     recursive = TRUE, full.names = TRUE)
##' fticr <- readMSData2(fticrf[1:2], msLevel. = 1)
##'
##' ## Perform the MSW feature detection on these:
##' p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
##'              SNR.method = "data.mean", winSize.noise = 500)
##' fticr <- detectFeatures(fticr, param = p)
##'
##' head(features(fticr))
##'
##' ## Now create the MzClustParam parameter object: we're assuming here that
##' ## both samples are from the same sample group.
##' p <- MzClustParam(sampleGroups = c(1, 1))
##'
##' fticr <- groupFeatures(fticr, param = p)
##'
##' featureGroups(fticr)
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
             msg <- validMsg(NULL, NULL)
             if (length(object@ppm) > 1 | any(object@ppm < 0))
                 msg <- validMsg(msg, paste0("'ppm' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@absMz) > 1 | any(object@absMz < 0))
                 msg <- validMsg(msg, paste0("'absMz' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@minFraction) > 1 | any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- validMsg(msg, paste0("'minFraction' has to be a ",
                                             "single positive number between ",
                                             "0 and 1!"))
             if (length(object@minSamples) > 1 | any(object@minSamples < 0))
                 msg <- validMsg(msg, paste0("'minSamples' has to be a ",
                                             "positive numeric of length 1!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })

## Main group.nearest documentation.
##' @title Feature alignment based on proximity in the mz-rt space
##'
##' @description This method is inspired by the alignment algorithm of mzMine
##' [Katajamaa 2006] and performs alignment based on proximity of features in the
##' space spanned by retention time and mz values.
##' The method creates first a \emph{master feature list} consisting of all
##' features from the sample in which most features were identified, and
##' starting from that, calculates distances to features from the sample with the
##' next most features. If features are closer than the defined threshold they
##' are grouped together.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{group}} methods. All of the settings to the alignment algorithm
##' can be passed with a \code{NearestFeaturesParam} object.
##'
##' @inheritParams groupFeatures-density
##'
##' @param mzVsRtBalance numeric(1) representing the factor by which mz values are
##' multiplied before calculating the (euclician) distance between two features.
##'
##' @param absMz numeric(1) maximum tolerated distance for mz values.
##'
##' @param absRt numeric(1) maximum tolerated distance for rt values.
##'
##' @param kNN numeric(1) representing the number of nearest neighbors to check.
##' 
##' @family feature alignment methods
##' 
##' @seealso The \code{\link{do_groupFeatures_nearest}} core
##' API function and \code{\link{group.nearest}} for the old user interface.
##' @seealso \code{\link{featureGroups}} and \code{\link{groupval,XCMSnExp-method}}
##' for methods to access feature grouping results.
##'
##' @name groupFeatures-nearest
##'
##' @references Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for
##' processing and visualization of mass spectrometry based molecular profile
##' data. \emph{Bioinformatics} 2006, 22:634-636. 
NULL
#> NULL

##' @description The \code{NearestFeaturesParam} class allows to specify all
##' settings for the feature alignment based on the \emph{nearest} algorithm.
##' Instances should be created with the \code{NearestFeaturesParam} constructor.
##'
##' @slot .__classVersion__,sampleGroups,mzVsRtBalance,absMz,absRt,kNN See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname groupFeatures-nearest
##'
##' @examples
##'
##' ## Create a NearestFeaturesParam object
##' p <- NearestFeaturesParam(kNN = 3)
##' p
##'
##' ##############################
##' ## feature detection and alignment.
##' ##
##' ## Below we perform first a feature detection (using the matchedFilter
##' ## method) on some of the test files from the faahKO package followed by
##' ## a feature alignment using the "nearest" method.
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' 
##' ## Reading 2 of the KO samples
##' raw_data <- readMSData2(fls[1:2])
##'
##' ## Perform the feature detection using the matchedFilter method.
##' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
##' res <- detectFeatures(raw_data, param = mfp)
##'
##' head(features(res))
##' ## The number of features identified per sample:
##' table(features(res)[, "sample"])
##'
##' ## Performing the feature alignment
##' p <- NearestFeaturesParam()
##' res <- groupFeatures(res, param = p)
##'
##' ## The results from the feature alignment:
##' featureGroups(res)
##'
##' ## Using the groupval method to extract a matrix with the intensities of
##' ## the feature groups per sample.
##' head(groupval(res, value = "into"))
##'
##' ## The process history:
##' processHistory(res)
setClass("NearestFeaturesParam",
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
             msg <- validMsg(NULL, NULL)
             if (length(object@mzVsRtBalance) > 1 |
                 any(object@mzVsRtBalance < 0))
                 msg <- validMsg(msg, paste0("'mzVsRtBalance' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@absMz) > 1 | any(object@absMz < 0))
                 msg <- validMsg(msg, paste0("'absMz' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@absRt) > 1 | any(object@absRt < 0))
                 msg <- validMsg(msg, paste0("'absRt' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@kNN) > 1 | any(object@kNN < 0))
                 msg <- validMsg(msg, paste0("'kNN' has to be a ",
                                             "positive numeric of length 1!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })



##' @title Retention time correction methods.
##'
##' @description The \code{adjustRtime} method(s) perform retention time
##' correction between chromatograms of different samples. These methods are
##' part of the modernized \code{xcms} user interface.
##'
##' The implemented feature alignment methods are:
##' \describe{
##' \item{featureGroups}{retention time correction based on aligment of feature
##' groups present in most/all samples.
##' See \code{\link{adjustRtime-featureGroups}} for more details.}
##' 
##' }
##' @name adjustRtime
##' @family retention time correction methods
##' @seealso \code{\link{retcor}} for the \emph{old} retention time correction 
##' methods.
##' @author Johannes Rainer
NULL
#> NULL

## Main retcor.peakgroups documentation.
##' @title Retention time correction based on alignment of house keeping feature
##' groups
##'
##' @description This method performs retention time adjustment based on the
##' alignment of feature groups present in all/most samples (hence corresponding
##' to house keeping compounds). First the retention time deviation of these
##' feature groups is described by fitting either a polynomial
##' (\code{smooth = "loess"}) or a linear (\code{smooth = "linear"}) model to the
##' data points. These models are subsequently used to adjust the retention time
##' of each spectrum in each sample.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{group}} methods. All of the settings to the alignment algorithm
##' can be passed with a \code{FeatureGroupsParam} object.
##'
##' @param minFraction numeric(1) between 0 and 1 defining the minimum required
##' fraction of samples in which features for the feature group were identified.
##' Feature groups passing this criteria will aligned across samples and retention
##' times of individual spectra will be adjusted based on this alignment. For
##' \code{minFraction = 1} the feature group has to contain features in all
##' samples of the experiment.
##' 
##' @param extraFeatures numeric(1) defining the maximal number of additional
##' features for all samples to be assigned to a feature group for retention time
##' correction. For a data set with 6 samples, \code{extraFeatures = 1} uses all
##' feature groups with a total feature count \code{<= 6 + 1}. The total feature
##' count is the total number of features being assigned to a feature group and
##' considers also multiple features within a sample being assigned to the group.
##'
##' @param smooth character defining the function to be used, to interpolate
##' corrected retention times for all feature groups. Either \code{"loess"} or
##' \code{"linear"}.
##'
##' @param span numeric(1) defining the degree of smoothing (if
##' \code{smooth = "loess"}). This parameter is passed to the internal call
##' to \code{\link{loess}}.
##'
##' @param family character defining the method to be used for loess smoothing.
##' Allowed values are \code{"gaussian"} and \code{"symmetric"}.See
##' \code{\link{loess}} for more information.
##' 
##' @family retention time correction methods
##' 
##' @seealso The \code{\link{do_adjustRtime_featureGroups}} core
##' API function and \code{\link{retcor.peakgroups}} for the old user interface.
##'
##' @name adjustRtime-featureGroups
##'
##' @author Colin Smith, Johannes Rainer
##' 
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
NULL
#> NULL

##' @description The \code{FeatureGroupsParam} class allows to specify all
##' settings for the retention time adjustment based on \emph{house keeping}
##' feature groups present in most samples.
##' Instances should be created with the \code{FeatureGroupsParam} constructor.
##'
##' @slot .__classVersion__,minFraction,extraFeatures,smooth,span,family See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname adjustRtime-featureGroups
##'
##' @examples
##' ##############################
##' ## feature detection and alignment.
##' ##
##' ## Below we perform first a feature detection (using the matchedFilter
##' ## method) on some of the test files from the faahKO package followed by
##' ## a feature alignment.
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' 
##' ## Reading 2 of the KO samples
##' raw_data <- readMSData2(fls[1:2])
##'
##' ## Perform the feature detection using the matchedFilter method.
##' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
##' res <- detectFeatures(raw_data, param = mfp)
##'
##' head(features(res))
##' ## The number of features identified per sample:
##' table(features(res)[, "sample"])
##'
##' ## Performing the feature alignment using the "feature density" method.
##' p <- FeatureDensityParam(sampleGroups = c(1, 1))
##' res <- groupFeatures(res, param = p)
##'
##' ## Perform the retention time adjustment using feature groups found in both
##' ## files.
##' fgp <- FeatureGroupsParam(minFraction = 1)
##' res <- adjustRtime(res, param = fgp)
##'
##' ## Any grouping information was dropped
##' hasAlignedFeatures(res)
##'
##' ## Plot the raw against the adjusted retention times.
##' plot(rtime(raw_data), rtime(res), pch = 16, cex = 0.25, col = fromFile(res))
##'
##' ## Adjusterd retention times can be accessed using
##' ## rtime(object, adjusted = TRUE) and adjustedRtime
##' all.equal(rtime(res), adjustedRtime(res))
##'
##' ## To get the raw, unadjusted retention times:
##' all.equal(rtime(res, adjusted = FALSE), rtime(raw_data))
##'
##' ## To extract the retention times grouped by sample/file:
##' rts <- rtime(res, bySample = TRUE)
setClass("FeatureGroupsParam",
         slots = c(minFraction = "numeric",
                   extraFeatures = "numeric",
                   smooth = "character",
                   span = "numeric",
                   family = "character"),
         contains = "Param",
         prototype = prototype(
             minFraction = 0.9,
             extraFeatures = 1,
             smooth = "loess",
             span = 0.2,
             family = "gaussian"
         ),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             if (length(object@minFraction) > 1 |
                 any(object@minFraction < 0) |
                 any(object@minFraction > 1))
                 msg <- validMsg(msg, paste0("'minFraction' has to be a single",
                                             " number between 0 and 1!"))
             if (length(object@extraFeatures) > 1 |
                 any(object@extraFeatures < 0))
                 msg <- validMsg(msg, paste0("'extraFeatures' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@span) > 1 | any(object@span < 0))
                 msg <- validMsg(msg, paste0("'span' has to be a ",
                                             "positive numeric of length 1!"))
             if (length(object@smooth) > 1 |
                 !all(object@smooth %in% c("loess", "linear")))
                 msg <- validMsg(msg, paste0("'smooth' has to be either \"",
                                             "loess\" or \"linear\"!"))
             if (length(object@family) > 1 |
                 !all(object@family %in% c("gaussian", "symmetric")))
                 msg <- validMsg(msg, paste0("'family' has to be either \"",
                                             "gaussian\" or \"symmetric\"!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })

##' @title Align retention times across samples using Obiwarp
##'
##' @description This method performs retention time adjustment using the
##' Obiwarp method [Prince 2006]. It is based on the code at
##' \url{http://obi-warp.sourceforge.net} but supports alignment of multiple
##' samples by aligning each against a \emph{center} sample. The alignment is
##' performed directly on the \code{\link{profile-matrix}} and can hence be
##' performed independently of the feature detection or feature grouping.
##'
##' @note These methods and classes are part of the updated and modernized
##' \code{xcms} user interface which will eventually replace the
##' \code{\link{retcor}} methods. All of the settings to the alignment algorithm
##' can be passed with a \code{ObiwarpParam} object.
##' 
##' @param binSize numeric(1) defining the bin size (in mz dimension) to be used
##' for the \emph{profile matrix} generation. See \code{step} parameter in
##' \code{\link{profile-matrix}} documentation for more details.
##'
##' @param centerSample integer(1) defining the index of the center sample in the
##' experiment. It defaults to \code{floor(median(1:length(fileNames(object))))}.
##'
##' @param response numeric(1) defining the \emph{responsiveness} of warping with
##' \code{response = 0} giving linear warping on start and end points and
##' \code{response = 100} warping using all bijective anchors.
##'
##' @param distFun character defining the distance function to be used. Allowed
##' values are \code{"cor"} (Pearson's correlation), \code{"cor_opt"} (calculate
##' only 10\% diagonal band of distance matrix; better runtime), \code{"cov"}
##' (covariance), \code{"prd"} (product) and \code{"euc"} (Euclidian distance).
##' The default value is \code{distFun = "cor_opt"}.
##'
##' @param gapInit numeric(1) defining the penalty for gap opening. The default
##' value for \code{gapInit} depends on the value of \code{distFun}: for
##' \code{distFun = "cor"} and \code{distFun = "cor_opt"} it is \code{0.3}, for
##' \code{distFun = "cov"} and \code{distFun = "prd"} \code{0.0} and for
##' \code{distFun = "euc"} \code{0.9}.
##'
##' @param gapExtend numeric(1) defining the penalty for gap enlargement. The
##' default value for \code{gapExtend} depends on the value of \code{distFun},
##' for \code{distFun = "cor"} and \code{distFun = "cor_opt"} it is \code{2.4},
##' for \code{distFun = "cov"} \code{11.7}, for \code{distFun = "euc"} \code{1.8}
##' and for \code{distFun = "prd"} {7.8}.
##'
##' @param factorDiag numeric(1) defining the local weight applied to diagonal
##' moves in the alignment.
##'
##' @param factorGap numeric(1) defining the local weight for gap moves in the
##' alignment.
##'
##' @param localAlignment logical(1) whether a local alignment should be
##' performed instead of the default global alignment.
##'
##' @param initPenalty numeric(1) defining the penalty for initiating an
##' alignment (for local alignment only).
##' 
##' @family retention time correction methods
##' 
##' @seealso \code{\link{retcor.obiwarp}} for the old user interface.
##'
##' @name adjustRtime-obiwarp
##'
##' @author Colin Smith, Johannes Rainer
##' 
##' @references
##' John T. Prince and Edward M. Marcotte. "Chromatographic Alignment of
##' ESI-LC-MS Proteomics Data Sets by Ordered Bijective Interpolated Warping"
##' \emph{Anal. Chem.} 2006, 78(17):6140-6152.

NULL
#> NULL

##' @description The \code{ObiwarpParam} class allows to specify all
##' settings for the retention time adjustment based on the \emph{obiwarp}
##' method. Class Instances should be created using the
##' \code{ObiwarpParam} constructor.
##'
##' @slot .__classVersion__,binSize,centerSample,response,distFun,gapInit,gapExtend,factorDiag,factorGap,localAlignment,initPenalty See corresponding parameter above. \code{.__classVersion__} stores
##' the version from the class. Slots values should exclusively be accessed
##' \emph{via} the corresponding getter and setter methods listed above.
##'
##' @rdname adjustRtime-obiwarp
##'
##' @examples
##' library(faahKO)
##' library(MSnbase)
##' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
##'            full.names = TRUE)
##' 
##' ## Reading 2 of the KO samples
##' raw_data <- readMSData2(fls[1:2])
##'
##' ## Perform retention time correction on the OnDiskMSnExp:
##' res <- adjustRtime(raw_data, param = ObiwarpParam())
##' 
##' ## As a result we get a numeric vector with the adjusted retention times for
##' ## all spectra.
##' head(res)
##'
##' ## We can split this by file to get the adjusted retention times for each
##' ## file
##' resL <- split(res, fromFile(raw_data))
##'
##' ##############################
##' ## Perform retention time correction on an XCMSnExp:
##' ##
##' ## Perform first the feature detection using the matchedFilter method.
##' mfp <- MatchedFilterParam(snthresh = 20, binSize = 1)
##' res <- detectFeatures(raw_data, param = mfp)
##'
##' ## Performing the retention time adjustment using obiwarp.
##' res_2 <- adjustRtime(res, param = ObiwarpParam())
##'
##' head(rtime(res_2))
##' head(rtime(raw_data))
##'
##' ## Also the retention times of the detected features were adjusted.
##' tail(features(res))
##' tail(features(res_2))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@binSize) > 1 |
                 any(object@binSize < 0))
                 msg <- validMsg(msg, paste0("'binSize' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@centerSample) > 1 |
                 any(object@centerSample < 0))
                 msg <- validMsg(msg, paste0("'centerSample' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@response) > 1 |
                 any(object@response < 0) |
                 any(object@response > 100))
                 msg <- validMsg(msg, paste0("'response' has to be a single ",
                                             " integer from 1 to 100!"))
             if (length(object@distFun) > 1 |
                 any(!(object@distFun %in% c("cor", "cor_opt", "cov", "euc",
                                             "prd"))))
                 msg <- validMsg(msg, paste0("'distFun' has to be one of \"cor\"",
                                             ", \"cor_opt\", \"cov\", \"euc\"",
                                             " or \"prd\"!"))
             if (length(object@gapInit) > 1 | any(object@gapInit < 0))
                 msg <- validMsg(msg, paste0("'gapInit' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@gapExtend) > 1 | any(object@gapExtend < 0))
                 msg <- validMsg(msg, paste0("'gapExtend' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@factorDiag) > 1 | any(object@factorDiag < 0))
                 msg <- validMsg(msg, paste0("'factorDiag' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@factorGap) > 1 | any(object@factorGap < 0))
                 msg <- validMsg(msg, paste0("'factorGap' has to be a positive",
                                             " numeric of length 1!"))
             if (length(object@localAlignment) > 1)
                 msg <- validMsg(msg, paste0("'localAlignment' has to be a ",
                                             "logical of length 1!"))
             if (length(object@initPenalty) > 1 | any(object@initPenalty < 0))
                 msg <- validMsg(msg, paste0("'initPenalty' has to be a positive",
                                             " numeric of length 1!"))
             if (is.null(msg))
                 return(TRUE)
             else
                 return(msg)
         })



##' @aliases MsFeatureData
##' @title Data container storing xcms preprocessing results
##'
##' @description The \code{MsFeatureData} class is designed to encapsule all
##' data related to the preprocessing of metabolomics data using the \code{xcms}
##' package, i.e. it contains a \code{matrix} with the features identified by the
##' feature detection, a \code{DataFrame} with the information on aligned
##' features across samples and a \code{list} with the adjusted retention times
##' per sample.
##'
##' @rdname XCMSnExp-class
setClass("MsFeatureData", contains = c("environment", "Versioned"),
         prototype = prototype(.xData = new.env(parent = emptyenv())))

.XCMS_REQ_FEATS_COLS <- c("mz", "mzmin", "mzmax", "rt", "rtmin",
                          "rtmax", "into", "sample")
.XCMS_REQ_FEATG_COLS <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                          "featureidx")

##' @aliases XCMSnExp
##' @title Data container storing xcms preprocessing results
##'
##' @description The \code{XCMSnExp} object is designed to contain all results
##' from metabolomics data preprocessing (feature detection, feature alignment
##' and retention time correction). The corresponding elements in the
##' \code{msFeatureData} slot are \code{"features"} (a \code{matrix}),
##' \code{"featureGroups"} (a \code{DataFrame}) and \code{"adjustedRtime"} (a
##' \code{list} of numeric vectors). Note that these should not be accessed
##' directly but rather \emph{via} their accessor methods. Along with the results,
##' the object contains the processing history that allow to track each
##' processing step along with the used settings. The object also directly
##' extends the \code{\link[MSnbase]{OnDiskMSnExp}} object hence allowing easy
##' access to the full data on which the feature detection was performed.
##'
##' Objects from this class should not be created directly, they are returned as
##' result from the \code{\link{detectFeatures}} method.
##'
##' \code{XCMSnExp} objects can be coerced into \code{\linkS4class{xcmsSet}}
##' objects using the \code{as} method.
##'
##' @note The \code{"features"} element in the \code{msFeatureData} slot is
##' equivalent to the \code{@peaks} slot of the \code{xcmsSet} object, the
##' \code{"featureGroups"} contains information from the \code{}
##'
##' @slot .processHistory \code{list} with \code{XProcessHistory} objects
##' tracking all individual analysis steps that have been performed.
##'
##' @slot msFeatureData \code{MsFeatureData} class extending \code{environment}
##' and containing the results from a feature detection (element
##' \code{"features"}), feature alignment (element \code{"featureGroups"}) and
##' retention time correction (element \code{""}) steps.
##'
##' @param object For \code{adjustedRtime}, \code{featureGroups},
##' \code{features}, \code{hasAdjustedRtime}, \code{hasAlignedFeatures} and
##' \code{hasDetectedFeatures} either a \code{MsFeatureData} or a \code{XCMSnExp}
##' object, for all other methods a \code{XCMSnExp} object.
##'
##' @param value For \code{adjustedRtime<-}: a \code{list} (length equal to the
##' number of samples) with numeric vectors representing the adjusted retention
##' times per scan.
##'
##' For \code{featureGroups<-}: a \code{DataFrame} with feature
##' alignment information. See return value for the \code{featureGroups} method
##' for the expected format.
##'
##' For \code{features<-}: a \code{matrix} with information on
##' detected features. See return value for the \code{features} method for the
##' expected format.
##'
##' @author Johannes Rainer
##'
##' @seealso \code{\linkS4class{xcmsSet}} for the old implementation.
##' @seealso \code{\link[MSnbase]{OnDiskMSnExp}}, \code{\link[MSnbase]{MSnExp}}
##' and \code{\link[MSnbase]{pSet}} for a complete list of inherited methods.
##' @seealso \code{\link{detectFeatures}} for available feature detection methods
##' returning a \code{XCMSnExp} object as a result.
##' @seealso \code{\link{groupFeatures}} for available feature grouping methods
##' and \code{\link{groupval,XCMSnExp-method}} for the method to extract feature
##' grouping results.
##' @seealso \code{\link{adjustRtime}} for retention time adjustment methods.
##'
##' @rdname XCMSnExp-class
##'
##' @examples
##'
##' ## Loading the data from 2 files of the faahKO package.
##' library(faahKO)
##' od <- readMSData2(c(system.file("cdf/KO/ko15.CDF", package = "faahKO"),
##'                     system.file("cdf/KO/ko16.CDF", package = "faahKO")))
##' ## Now we perform a feature detection on this data set using the
##' ## matched filter method. We are tuning the settings such that it performs
##' ## faster.
##' mfp <- MatchedFilterParam(binSize = 4)
##' xod <- detectFeatures(od, param = mfp)
##'
##' ## The results from the feature detection are now stored in the XCMSnExp
##' ## object
##' xod
##'
##' ## The detected features can be accessed with the features method.
##' head(features(xod))
##'
##' ## The settings of the feature detection can be accessed with the
##' ## processHistory method
##' processHistory(xod)
##'
##' ## Also the parameter class for the feature detection can be accessed
##' processParam(processHistory(xod)[[1]])
##'
##' ## The XCMSnExp inherits all methods from the pSet and OnDiskMSnExp classes
##' ## defined in Bioconductor's MSnbase package. To access the (raw) retention
##' ## time for each spectrum we can use the rtime method. Setting bySample = TRUE
##' ## would cause the retention times to be grouped by sample
##' head(rtime(xod))
##'
##' ## Similarly it is possible to extract the mz values or the intensity values
##' ## using the mz and intensity method, respectively, also with the option to
##' ## return the results grouped by sample instead of the default, which is
##' ## grouped by spectrum. Finally, to extract all of the data we can use the
##' ## spectra method which returns Spectrum objects containing all raw data.
##' ## Note that all these methods read the information from the original input
##' ## files and subsequently apply eventual data processing steps to them.
##' head(mz(xod, bySample = TRUE))
##'
##' ## Reading all data
##' spctr <- spectra(xod)
##' ## To get all spectra of the first file we can split them by file
##' head(split(spctr, fromFile(xod))[[1]])
##'
##' ############
##' ## Filtering
##' ##
##' ## XCMSnExp objects can be filtered by file, retention time, mz values or
##' ## MS level. For some of these filter preprocessing results (mostly
##' ## retention time correction and feature alignment results) will be dropped.
##' ## Below we filter the XCMSnExp object by file to extract the results for
##' ## only the second file.
##' xod_2 <- filterFile(xod, file = 2)
##' xod_2
##'
##' ## Now the objects contains only the idenfified features for the second file
##' head(features(xod_2))
##'
##' head(features(xod)[features(xod)[, "sample"] == 2, ])
##'
##' ##########
##' ## Coercing to an xcmsSet object
##' ##
##' ## We can also coerce the XCMSnExp object into an xcmsSet object:
##' xs <- as(xod, "xcmsSet")
##' head(peaks(xs))
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
             msg <- validMsg(NULL, NULL)
             if (length(object@.processHistory) > 0) {
                 isOK <- unlist(lapply(object@.processHistory, function(z) {
                     return(inherits(z, "ProcessHistory"))
                 }))
                 if (!all(isOK))
                     msg <- validMsg(msg, paste0("Only 'ProcessHistory' ",
                                                 "objects are allowed in slot ",
                                                 ".processHistory!"))
             }
             ## TODO @jo add checks:
             ## 1) call validMsFeatureData
             msg <- validMsg(msg, validateMsFeatureData(object@msFeatureData))
             if (!is.null(msg)) return(msg)
             ## 2) features[, "sample"] is within 1:number of samples
             if (any(ls(object@msFeatureData) == "features")) {
                 if (!all(object@msFeatureData$features[, "sample"] %in%
                          1:length(fileNames(object))))
                     msg <- validMsg(msg, paste0("The number of available ",
                                                 "samples does not match with ",
                                                 "the sample assignment of ",
                                                 "features in the 'features' ",
                                                 "element of the msFeatureData ",
                                                 "slot!"))
             }
             ## 3) Check that the length of the adjustedRtime matches!
             if (any(ls(object@msFeatureData) == "adjustedRtime")) {
                 rt <- rtime(object, bySample = TRUE)
                 if (length(rt) != length(object@msFeatureData$adjustedRtime)) {
                     msg <- validMsg(msg, paste0("The number of numeric vectors",
                                                 " in the 'adjustedRtime' element",
                                                 " of the msFeatureData slot does",
                                                 " not match the number of",
                                                 " samples!"))
                 } else {
                     if (any(lengths(rt) !=
                             lengths(object@msFeatureData$adjustedRtime)))
                         msg <- validMsg(msg,
                                         paste0("The lengths of the numeric ",
                                                "vectors in the 'adjustedRtime'",
                                                " element of the msFeatureData ",
                                                "slot does not match the number",
                                                " of scans per sample!"))
                 }
             }
             ## 3) If we've got features, check that we have also a related
             ##    processing history step.
             if (is.null(msg))
                 return(TRUE)
             else return(msg)
         }
)

