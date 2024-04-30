## All class definitions should go in here.
#' @include AllGenerics.R functions-XChromatogram.R functions-XChromatograms.R

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
.PROCSTEP.PEAK.REFINEMENT <- "Peak refinement"
.PROCSTEP.PEAK.GROUPING <- "Peak grouping"
.PROCSTEP.RTIME.CORRECTION <- "Retention time correction"
.PROCSTEP.PEAK.FILLING <- "Missing peak filling"
.PROCSTEP.CALIBRATION <- "Calibration"
.PROCSTEP.FEATURE.GROUPING <- "Feature grouping"
.PROCSTEP.FEATURE.FILTERING <- "Feature filtering"
.PROCSTEPS <- c(
    .PROCSTEP.UNKNOWN,
    .PROCSTEP.PEAK.DETECTION,
    .PROCSTEP.PEAK.REFINEMENT,
    .PROCSTEP.PEAK.GROUPING,
    .PROCSTEP.RTIME.CORRECTION,
    .PROCSTEP.PEAK.FILLING,
    .PROCSTEP.CALIBRATION,
    .PROCSTEP.FEATURE.GROUPING,
    .PROCSTEP.FEATURE.FILTERING
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
setClassUnion("ParamOrNULL", c("Param", "NULL"))

#' @aliases GenericParam Param class:Param Param-class
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
#' @slot msLevel: \code{integer} definining the MS level(s) on which the
#'     analysis was performed.
#'
#' @rdname ProcessHistory-class
setClass("XProcessHistory",
         slots = c(
             param = "ParamOrNULL",
             msLevel = "integer"
         ),
         contains = "ProcessHistory",
         prototype = prototype(
             param = NULL,
             msLevel = NA_integer_
         ),
         validity = function(object) {
             msg <- character()
             if (length(object@param) > 0)
                 if (!is(object@param, "Param"))
                     msg <- c(msg,
                              paste0("Only objects from type 'Param' ",
                                     "allowed in slot '@param'! I got ",
                                     class(object@param)))
             if (!is.na(msLevel(object)))
                 if (msLevel(object) < 0)
                     msg <- c(msg, "msLevel has to be a positive integer")
             if (length(msg)) msg
             else TRUE
         })

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
#'     dimension required for peaks with overlapping retention times; can be
#'     negative to allow overlap. During peak post-processing, peaks
#'     defined to be overlapping are reduced to the one peak with the largest
#'     signal.
#'
#' @param fitgauss \code{logical(1)} whether or not a Gaussian should be fitted
#'     to each peak. This affects mostly the retention time position of the
#'     peak.
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
#'     In detail, a first rough estimate of the noise is calculated and peak
#'     detection is performed only in regions in which multiple sequential
#'     signals are higher than this first estimated baseline/noise level.
#'
#' @param roiScales Optional numeric vector with length equal to \code{roiList}
#'     defining the scale for each region of interest in \code{roiList} that
#'     should be used for the centWave-wavelets.
#'
#' @param extendLengthMSW Option to force centWave to use all scales when
#' running centWave rather than truncating with the EIC length. Uses the "open"
#' method to extend the EIC to a integer base-2 length prior to being passed to
#' \code{convolve} rather than the default "reflect" method. See
#' https://github.com/sneumann/xcms/issues/445 for more information.
#'
#' @param verboseBetaColumns Option to calculate two additional metrics of peak
#' quality via comparison to an idealized bell curve. Adds \code{beta_cor} and
#' \code{beta_snr} to the \code{chromPeaks} output, corresponding to a Pearson
#' correlation coefficient to a bell curve with several degrees of skew as well
#' as an estimate of signal-to-noise using the residuals from the best-fitting
#' bell curve. See https://github.com/sneumann/xcms/pull/685 and
#' https://doi.org/10.1186/s12859-023-05533-4 for more information.
#'
#' @param x The parameter object.
#'
#' @details
#'
#' The centWave algorithm is most suitable for high resolution
#' LC/\{TOF,OrbiTrap,FTICR\}-MS data in centroid mode. In the first phase
#' the method identifies \emph{regions of interest} (ROIs) representing
#' mass traces that are characterized as regions with less than \code{ppm}
#' m/z deviation in consecutive scans in the LC/MS map. In detail, starting
#' with a single m/z, a ROI is extended if a m/z can be found in the next scan
#' (spectrum) for which the difference to the mean m/z of the ROI is smaller
#' than the user defined \code{ppm} of the m/z. The mean m/z of the ROI is then
#' updated considering also the newly included m/z value.
#'
#' These ROIs are then, after some cleanup, analyzed using continuous wavelet
#' transform (CWT) to locate chromatographic peaks on different scales.
#' The first analysis step is skipped, if regions of interest are passed
#' \emph{via} the \code{param} parameter.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports peak detection on
#'     \code{\link{OnDiskMSnExp}} objects (defined in the \code{MSnbase}
#'     package). All of the settings to the centWave algorithm can be passed
#'     with a \code{CentWaveParam} object.
#'
#' @family peak detection methods
#'
#' @seealso
#'
#' The \code{\link{do_findChromPeaks_centWave}} core API function and
#' \code{\link{findPeaks.centWave}} for the old user interface.
#'
#' \code{\link{peaksWithCentWave}} for functions to perform centWave peak
#' detection in purely chromatographic data.
#'
#' @references
#' Ralf Tautenhahn, Christoph Böttcher, and Steffen Neumann "Highly
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
#' @slot ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales,extendLengthMSW,verboseBetaColumns See corresponding parameter above. Slots values should exclusively be accessed
#' \emph{via} the corresponding getter and setter methods listed above.
#'
#' @rdname findChromPeaks-centWave
#'
#' @examples
#'
#' ## Create a CentWaveParam object. Note that the noise is set to 10000 to
#' ## speed up the execution of the example - in a real use case the default
#' ## value should be used, or it should be set to a reasonable value.
#' cwp <- CentWaveParam(ppm = 20, noise = 10000, prefilter = c(3, 10000))
#' ## Change snthresh parameter
#' snthresh(cwp) <- 25
#' cwp
#'
#' ## Perform the peak detection using centWave on some of the files from the
#' ## faahKO package. Files are read using the `readMsExperiment` function
#' ## from the MsExperiment package
#' library(faahKO)
#' library(xcms)
#' library(MsExperiment)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMsExperiment(fls[1])
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
             roiScales = "numeric",
             extendLengthMSW = "logical",
             verboseBetaColumns = "logical"
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
             roiScales = numeric(),
             extendLengthMSW = FALSE,
             verboseBetaColumns = FALSE
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
#'     detection on
#'     \code{\link{OnDiskMSnExp}} objects (defined in the
#'     \code{MSnbase} package). All of the settings to the matchedFilter
#'     algorithm can be passed with a \code{MatchedFilterParam} object.
#'
#' @inheritParams imputeLinInterpol
#'
#' @inheritParams findChromPeaks-centWave
#'
#' @family peak detection methods
#'
#' @seealso
#'
#' The \code{\link{do_findChromPeaks_matchedFilter}} core API function
#' and \code{\link{findPeaks.matchedFilter}} for the old user interface.
#'
#' \code{\link{peaksWithMatchedFilter}} for functions to perform matchedFilter
#' peak detection in purely chromatographic data.
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
#' @slot binSize,impute,baseValue,distance,fwhm,sigma,max,snthresh,steps,mzdiff,index
#'     See corresponding parameter above. Slots values should exclusively
#'     be accessed \emph{via} the corresponding getter and setter methods listed
#'     above.
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
#' raw_data <- readMSData(fls[1], mode = "onDisk")
#' ## Perform the chromatographic peak detection using the settings defined
#' ## above. Note that we are also disabling parallel processing in this
#' ## example by registering a "SerialParam"
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
#'     on high resolution LC/(OrbiTrap, TOF)-MS data in centroid mode.
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
#'     detection on
#'     \code{\link{OnDiskMSnExp}} objects (defined in the
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
#' @slot ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,criticalValue,consecMissedLimit,unions,checkBack,withWave
#'      See corresponding parameter above. Slots values should
#'      exclusively be accessed \emph{via} the corresponding getter and setter
#'      methods listed above.
#'
#' @rdname findChromPeaks-massifquant
#'
#' @examples
#'
#' ## Create a MassifquantParam object.
#' mqp <- MassifquantParam()
#' ## Change snthresh prefilter parameters
#' snthresh(mqp) <- 30
#' prefilter(mqp) <- c(6, 10000)
#' mqp
#'
#' ## Perform the peak detection using massifquant on the files from the
#' ## faahKO package. Files are read using the readMSData from the MSnbase
#' ## package
#' library(faahKO)
#' library(MSnbase)
#' fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
#'            full.names = TRUE)
#' raw_data <- readMSData(fls[1], mode = "onDisk")
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
#'     \code{\link{peakDetectionCWT}} and
#'     \code{\link{tuneInPeakInfo}} functions. See the
#'     \emph{xcmsDirect} vignette for more information.
#'
#' @note These methods and classes are part of the updated and modernized
#'     \code{xcms} user interface which will eventually replace the
#'     \code{\link{findPeaks}} methods. It supports peak detection on
#'     \code{\link{OnDiskMSnExp}}
#'     objects (defined in the \code{MSnbase} package). All of the settings
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
#' @slot snthresh,verboseColumns,scales,nearbyPeak,peakScaleRange,ampTh,minNoiseLevel,ridgeLength,peakThr,tuneIn,addParams
#'      See corresponding parameter above.
#'
#' @rdname findPeaks-MSW
#'
#' @examples
#'
#' library(MSnbase)
#' ## Create a MSWParam object
#' mp <- MSWParam()
#' ## Change snthresh parameter
#' snthresh(mp) <- 15
#' mp
#'
#' ## Loading a small subset of direct injection, single spectrum files
#' library(msdata)
#' fticrf <- list.files(system.file("fticr-mzML", package = "msdata"),
#'                     recursive = TRUE, full.names = TRUE)
#' fticr <- readMSData(fticrf[1], msLevel. = 1, mode = "onDisk")
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
#'     detection on
#'     \code{\link{OnDiskMSnExp}} objects (defined in the
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
#' @slot ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselineCheck,roiScales,extendLengthMSW,verboseBetaColumns,snthreshIsoROIs,maxCharge,maxIso,mzIntervalExtension,polarity
#'      See corresponding parameter above.
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

setClass("PeakDensityParam",
         slots = c(sampleGroups = "ANY",
                   bw = "numeric",
                   minFraction = "numeric",
                   minSamples = "numeric",
                   binSize = "numeric",
                   maxFeatures = "numeric",
                   ppm = "numeric"),
         contains = "Param",
         prototype = prototype(
             sampleGroups = numeric(),
             bw = 30,
             minFraction = 0.5,
             minSamples = 1,
             binSize = 0.25,
             ppm = 0,
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

setClass("PeakGroupsParam",
         slots = c(minFraction = "numeric",
                   extraPeaks = "numeric",
                   smooth = "character",
                   span = "numeric",
                   family = "character",
                   peakGroupsMatrix = "matrix",
                   subset = "integer",
                   subsetAdjust = "character"),
         contains = "Param",
         prototype = prototype(
             minFraction = 0.9,
             extraPeaks = 1,
             smooth = "loess",
             span = 0.2,
             family = "gaussian",
             peakGroupsMatrix = matrix(ncol = 0, nrow = 0),
             subset = integer(),
             subsetAdjust = "average"
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

setClass("LamaParama",
         slots = c(lamas = "matrix",
                   method = "character",
                   span = "numeric",
                   outlierTolerance = "numeric",
                   zeroWeight = "numeric",
                   ppm = "numeric",
                   tolerance = "numeric",
                   toleranceRt = "numeric",
                   bs = "character",
                   rtMap = "list",
                   nChromPeaks = "numeric"),
         contains = "Param",
         prototype = prototype(
             lamas = matrix(ncol = 2, nrow = 0),
             method = "loess",
             span = 0.5,
             outlierTolerance = 3,
             zeroWeight = 10,
             ppm = 20,
             tolerance = 0,
             toleranceRt = 20,
             bs = "tp",
             rtMap = list(),
             nChromPeaks = numeric()),
         validity = function(object) {
             msg <- NULL
             if (!nrow(object@lamas))
                 msg <- c(msg, paste0("'lamas' cannot be empty"))
             else {
             }
             if (length(object@method) > 1 |
                 !all(object@method %in% c("gam", "loess")))
                 msg <- c(msg, paste0("'method' has to be either \"",
                                      "gam\" or \"loess\"!"))
             msg
         })


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
                   initPenalty = "numeric",
                   subset = "integer",
                   subsetAdjust = "character",
                   rtimeDifferenceThreshold = "numeric"),
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
             initPenalty = 0,
             subset = integer(),
             subsetAdjust = "average",
             rtimeDifferenceThreshold = 5),
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

#' @slot expandMz,expandRt,ppm,fixedMz,fixedRt See corresponding parameter
#'      above.
#'
#' @rdname fillChromPeaks
setClass("FillChromPeaksParam",
         slots = c(expandMz = "numeric",
                   expandRt = "numeric",
                   ppm = "numeric",
                   fixedMz = "numeric",
                   fixedRt = "numeric"),
         contains = "Param",
         prototype = prototype(
             expandMz = 0,
             expandRt = 0,
             ppm = 0,
             fixedMz = 0,
             fixedRt = 0
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
             if (length(object@fixedMz) > 1)
                 msg <- c(msg, "'fixedMz' has to be a numeric of length 1")
             if (length(object@fixedRt) > 1)
                 msg <- c(msg, "'fixedRt' has to be a numeric of length 1")
             if (length(msg))
                 msg
             else TRUE
         }
         )

#' @rdname fillChromPeaks
#'
#' @slot rtmin,rtmax,mzmin,mzmax See corresponding parameter above.
setClass("ChromPeakAreaParam",
         slots = c(rtmin = "function",
                   rtmax = "function",
                   mzmin = "function",
                   mzmax = "function"),
         contains = "Param"
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
#' @noRd
#'
#' @rdname XCMSnExp-class
setClass("MsFeatureData", contains = c("environment"),
         prototype = prototype(.xData = new.env(parent = emptyenv())))

.REQ_PEAKS_COLS <- c("mz", "mzmin", "mzmax", "rt", "rtmin",
                     "rtmax", "into", "sample")
.REQ_PEAKG_COLS <- c("mzmed", "mzmin", "mzmax", "rtmed", "rtmin", "rtmax",
                     "peakidx")

#' @aliases XCMSnExp updateObject,XCMSnExp-method
#'
#' @title Data container storing xcms preprocessing results
#'
#' @description
#'
#' The \code{XCMSnExp} object is a container for the results of a G/LC-MS
#' data preprocessing that comprises chromatographic peak detection, alignment
#' and correspondence. These results can be accessed with the \code{chromPeaks},
#' \code{adjustedRtime} and \code{featureDefinitions} functions; see below
#' (after the Usage, Arguments, Value and Slots sections) for more details).
#' Along with the results, the object contains the processing history that
#' allows to track each processing step along with the used settings. This
#' can be extracted with the \code{\link{processHistory}} method.
#' \code{XCMSnExp} objects, by directly extending the
#' \code{\link{OnDiskMSnExp}} object from the \code{MSnbase} package, inherit
#' all of its functionality and allows thus an easy access to the full raw
#' data at any stage of an analysis.
#' To support interaction with packages requiring the \emph{old} objects,
#' \code{XCMSnExp} objects can be coerced into \code{\linkS4class{xcmsSet}}
#' objects using the \code{as} method (see examples below). All
#' preprocessing results will be passed along to the resulting
#' \code{xcmsSet} object.
#'
#' General functions for \code{XCMSnExp} objects are (see further below for
#' specific function to handle chromatographic peak data, alignment and
#' correspondence results):
#'
#' @section Chromatographic peak data:
#'
#' Chromatographic peak data is added to an \code{XCMSnExp} object by the
#' \code{\link{findChromPeaks}} function. Functions to access chromatographic
#' peak data are:
#'
#' \itemize{
#' \item \code{hasChromPeaks} whether chromatographic peak data is available,
#' see below for help of the function.
#'
#' \item \code{chromPeaks} access chromatographic peaks (see below for help).
#'
#' \item \code{dropChromPeaks} remove chromatographic peaks (see below for
#' help).
#'
#' \item \code{dropFilledChromPeaks} remove filled-in peaks (see below for
#' help).
#'
#' \item \code{\link{fillChromPeaks}} fill-in missing peaks (see respective
#' help page).
#'
#' \item \code{\link{plotChromPeaks}} plot identified peaks for a file (see
#' respective help page).
#'
#' \item \code{\link{plotChromPeakImage}} plot distribution of peaks along the
#' retention time axis (see respective help page).
#'
#' \item \code{\link{highlightChromPeaks}} add chromatographic peaks to an
#' existing plot of a \code{\link{Chromatogram}} (see respective help page).
#'
#' }
#'
#'
#' @section Adjusted retention times:
#'
#' Adjusted retention times are stored in an \code{XCMSnExp} object besides the
#' original, raw, retention times, allowing to switch between raw and adjusted
#' times. It is also possible to replace the raw retention times with the
#' adjusted ones with the \code{\link{applyAdjustedRtime}}. The adjusted
#' retention times are added to an \code{XCMSnExp} by the
#' \code{\link{adjustRtime}} function. All functions related to the access of
#'  adjusted retention times are:
#'
#' \itemize{
#'
#' \item \code{hasAdjustedRtime} whether adjusted retention times are available
#' (see below for help).
#'
#' \item \code{dropAdjustedRtime} remove adjusted retention times (see below
#' for help).
#'
#' \item \code{\link{applyAdjustedRtime}} replace the raw retention times with
#' the adjusted ones (see respective help page).
#'
#' \item \code{\link{plotAdjustedRtime}} plot differences between adjusted and
#' raw retention times (see respective help page).
#'
#' }
#'
#'
#' @section Correspondence results, features:
#'
#' The correspondence analysis (\code{\link{groupChromPeaks}}) adds the feature
#' definitions to an \code{XCMSnExp} object. All functions related to these are
#' listed below:
#'
#' \itemize{
#'
#' \item \code{hasFeatures} whether correspondence results are available (see
#' below for help).
#'
#' \item \code{featureDefinitions} access the definitions of the features (see
#' below for help).
#'
#' \item \code{dropFeatureDefinitions} remove correspondence results (see below
#' for help).
#'
#' \item \code{\link{featureValues}} access values for features (see respective
#' help page).
#'
#' \item \code{\link{featureSummary}} perform a simple summary of the defined
#' features (see respective help page).
#'
#' \item \code{\link{overlappingFeatures}} identify features that are
#' overlapping or close in the m/z - rt space (see respective help page).
#'
#' \item \code{\link{quantify}} extract feature intensities and put them, along
#' with feature definitions and phenodata information, into a
#' \code{\link{SummarizedExperiment}}. See help page for details.
#' }
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
#'     This object should not be manipulated directly.
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
#'
#' @author Johannes Rainer
#'
#' @seealso \code{\linkS4class{xcmsSet}} for the old implementation.
#'     \code{\link{OnDiskMSnExp}}, \code{\link{MSnExp}}
#'     and \code{\link{pSet}} for a complete list of inherited methods.
#'
#'     \code{\link{findChromPeaks}} for available peak detection methods
#'     returning a \code{XCMSnExp} object as a result.
#'
#'     \code{\link{groupChromPeaks}} for available peak grouping
#'     methods and \code{\link{featureDefinitions}} for the method to extract
#'     the feature definitions representing the peak grouping results.
#'     \code{\link{adjustRtime}} for retention time adjustment methods.
#'
#'     \code{\link{chromatogram}} to extract MS data as
#'     \code{\link{Chromatogram}} objects.
#'
#'     \code{\link{as}} (\code{as(x, "data.frame")}) in the \code{MSnbase}
#'     package for the method to extract MS data as \code{data.frame}s.
#'
#'     \code{\link{featureSummary}} to calculate basic feature summaries.
#'
#'     \code{\link{featureChromatograms}} to extract chromatograms for each
#'     feature.
#'
#'     \code{\link{chromPeakSpectra}} to extract MS2 spectra with the m/z of
#'     the precursor ion within the m/z range of a peak and a retention time
#'     within its retention time range.
#'
#'     \code{\link{featureSpectra}} to extract MS2 spectra associated with
#'     identified features.
#'
#' @rdname XCMSnExp-class
#'
#' @examples
#'
#' ## Load a test data set with detected peaks
#' library(MSnbase)
#' data(faahko_sub)
#' ## Update the path to the files for the local system
#' dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")
#'
#' ## Disable parallel processing for this example
#' register(SerialParam())
#'
#' ## The results from the peak detection are now stored in the XCMSnExp
#' ## object
#' faahko_sub
#'
#' ## The detected peaks can be accessed with the chromPeaks method.
#' head(chromPeaks(faahko_sub))
#'
#' ## The settings of the chromatographic peak detection can be accessed with
#' ## the processHistory method
#' processHistory(faahko_sub)
#'
#' ## Also the parameter class for the peak detection can be accessed
#' processParam(processHistory(faahko_sub)[[1]])
#'
#' ## The XCMSnExp inherits all methods from the pSet and OnDiskMSnExp classes
#' ## defined in Bioconductor's MSnbase package. To access the (raw) retention
#' ## time for each spectrum we can use the rtime method. Setting bySample = TRUE
#' ## would cause the retention times to be grouped by sample
#' head(rtime(faahko_sub))
#'
#' ## Similarly it is possible to extract the mz values or the intensity values
#' ## using the mz and intensity method, respectively, also with the option to
#' ## return the results grouped by sample instead of the default, which is
#' ## grouped by spectrum. Finally, to extract all of the data we can use the
#' ## spectra method which returns Spectrum objects containing all raw data.
#' ## Note that all these methods read the information from the original input
#' ## files and subsequently apply eventual data processing steps to them.
#' mzs <- mz(faahko_sub, bySample = TRUE)
#' length(mzs)
#' lengths(mzs)
#'
#' ## The full data could also be read using the spectra data, which returns
#' ## a list of Spectrum object containing the mz, intensity and rt values.
#' ## spctr <- spectra(faahko_sub)
#' ## To get all spectra of the first file we can split them by file
#' ## head(split(spctr, fromFile(faahko_sub))[[1]])
#'
#' ############
#' ## Filtering
#' ##
#' ## XCMSnExp objects can be filtered by file, retention time, mz values or
#' ## MS level. For some of these filter preprocessing results (mostly
#' ## retention time correction and peak grouping results) will be dropped.
#' ## Below we filter the XCMSnExp object by file to extract the results for
#' ## only the second file.
#' xod_2 <- filterFile(faahko_sub, file = 2)
#' xod_2
#'
#' ## Now the objects contains only the idenfified peaks for the second file
#' head(chromPeaks(xod_2))
#'
#' ##########
#' ## Coercing to an xcmsSet object
#' ##
#' ## We can also coerce the XCMSnExp object into an xcmsSet object:
#' xs <- as(faahko_sub, "xcmsSet")
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
                 if (!any(ls(object@msFeatureData) == "chromPeakData"))
                     return(paste0("Missing 'chromPeakData'. Please update",
                                   " the object with 'updateObject'"))
             }
             ## 3) Check that the length of the adjustedRtime matches!
             if (any(ls(object@msFeatureData) == "adjustedRtime")) {
                 rt <- rtime(object, bySample = TRUE, adjusted = FALSE)
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

.CHROMPEAKS_REQ_NAMES <- c("rt", "rtmin", "rtmax", "into", "maxo", "sn")
.CHROMPEAKDATA_REQ_NAMES <- c("ms_level", "is_filled")
setClass("XChromatogram",
         slots = c(chromPeaks = "matrix",
                   chromPeakData = "DataFrame"),
         prototype = prototype(
             chromPeaks = matrix(nrow = 0, ncol = length(.CHROMPEAKS_REQ_NAMES),
                                 dimnames = list(character(),
                                                 .CHROMPEAKS_REQ_NAMES)),
             chromPeakData = DataFrame(ms_level = integer(),
                                       is_filled = logical())
         ),
         contains = "Chromatogram",
         validity = .validXChromatogram)

setClass("XChromatograms",
         slots = c(.processHistory = "list",
                   featureDefinitions = "DataFrame"),
         prototype = prototype(.processHistory = list(),
                               featureDefinitions = DataFrame()),
         contains = "MChromatograms",
         validity = .validXChromatograms)

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

setClass("CleanPeaksParam",
         slots = c(maxPeakwidth = "numeric"),
         contains = "Param",
         prototype = prototype(
             maxPeakwidth = 10),
         validity = function(object) {
             msg <- character()
             if (length(object@maxPeakwidth) > 1 || object@maxPeakwidth < 0)
                 msg <- c(msg, paste0("'maxPeakwidth' has to be a positive ",
                                      "number of length 1"))
             if (length(msg))
                 msg
             else TRUE
         })

setClass("MergeNeighboringPeaksParam",
         slots = c(expandRt = "numeric",
                   expandMz = "numeric",
                   ppm = "numeric",
                   minProp = "numeric"),
         contains = "Param",
         prototype = prototype(
             expandRt = 2.0,
             expandMz = 0.0,
             ppm = 10.0,
             minProp = 0.75),
         validity = function(object) {
             msg <- character()
             if (length(object@expandRt) > 1 || !is.finite(object@expandRt))
                 msg <- c(msg, paste0("'expandRt' has to be a (defined) ",
                                      "numeric of length 1"))
             if (length(object@expandMz) > 1 || !is.finite(object@expandMz))
                 msg <- c(msg, paste0("'expandMz' has to be a (defined) ",
                                      "numeric of length 1"))
             if (length(object@ppm) > 1 || !is.finite(object@ppm) ||
                 object@ppm < 0)
                 msg <- c(msg, paste0("'ppm' has to be a positive numeric ",
                                      "of length 1"))
             if (length(object@minProp) > 1 || !is.finite(object@minProp) ||
                 object@minProp < 0)
                 msg <- c(msg, paste0("'minProp' has to be a positive ",
                                      "number of length 1"))
             if (length(msg))
                 msg
             else TRUE
         })

setClass("FilterIntensityParam",
         slots = c(threshold = "numeric",
                   nValues = "integer",
                   value = "character"),
         contains = "Param",
         prototype = prototype(
             threshold = 0,
             nValues = 1L,
             value = "maxo"),
         validity = function(object) {
             msg <- character()
             if (length(object@threshold) > 1 || object@threshold < 0)
                 msg <- c(msg, paste0("'threshold' has to be a positive ",
                                      "number of length 1"))
             if (length(object@nValues) > 1 || object@nValues < 1)
                 msg <- c(msg, paste0("'nValues' has to be a positive ",
                                      "number of length 1"))
             if (length(object@value) > 1)
                 msg <- c(msg, paste0("'value' has to be a character ",
                                      "of length 1"))
             if (length(msg))
                 msg
             else TRUE
         })
