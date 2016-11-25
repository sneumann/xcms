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
                     msg <- validMsg(msg, "Slot '.processHistory' should only contain 'ProcessHistory' objects!")
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
.PROCSTEPS <- c(
    .PROCSTEP.UNKNOWN,
    .PROCSTEP.FEATURE.DETECTION
)

############################################################
## ProcessHistory
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
             error = NULL,
             new("Versioned", versions = c(ProcessHistory = "0.0.2"))
         ),
         validity = function(object) {
             msg <- validMsg(NULL, NULL)
             ## check type:
             if (!any(object@type == .PROCSTEPS))
                 msg <- validMsg(msg, paste0("Got invalid type '", object@type,
                                             "'! Allowd are: ",
                                             paste0("\"", .PROCSTEPS, "\"",
                                                    collapse = ", ")))
             if (is.null(msg)) TRUE
             else msg
         }
         )

## BasicParam class
## CentWaveParam
setClass("Param",
         representation = representation("VIRTUAL"),
         contains = c("Versioned"))

## General detectFeatures method.
##' @title Feature detection methods.
##'
##' @description The \code{detectFeature} methods are part of the modernized
##' \code{xcms} user interface.
##'
##' The implemented feature detection methods are:
##' \describe{
##' \item{centWave}{: feature detection using the \emph{centWave} method.
##' See \code{\link{centWave}} for more details.}
##'
##' \item{matchedFilter}{: peak detection in chromatographic space. See
##' \code{\link{matchedFilter}} for more details.}
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
##' ## Create a CentWaveParam object
##' cwp <- CentWaveParam(ppm = 20)
##' ## Change snthresh parameter
##' snthresh(cwp) <- 5
##' cwp
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
             ## Check the values.
             ## ppm positive numeric of length 1.
             if (length(object@ppm) != 1 | object@ppm < 0)
                 msg <- validMsg(msg, paste0("'ppm' has to be positive numeric",
                                             " of length 1."))
             ## peakwidth
             if (length(object@peakwidth) != 2 | any(object@peakwidth < 0))
                 msg <- validMsg(msg, paste0("'peakwidth' has to be a numeric",
                                             " of length 2 with only positive",
                                             " values."))
             ## snthresh positive numeric of length 1.
             if (length(object@snthresh) != 1 | object@snthresh < 0)
                 msg <- validMsg(msg, paste0("'snthresh' has to be a positive",
                                             " numeric of length 1."))
             ## prefilter: numeric of length 2.
             if (length(object@prefilter) != 2)
                 msg <- validMsg(msg, paste0("'prefilter' has to be a numeric",
                                             " of length 2."))
             ## mzCenterFun: check if method exists.
             allowed_vals <- c("wMean", "mean", "apex", "wMeanApex3",
                               "meanApex3")
             if (!(object@mzCenterFun) %in% allowed_vals)
                 msg <- validMsg(msg, paste0("'mzCenterFun' has to be one of ",
                                             paste0("'", allowed_vals, "'",
                                             collapse = ", "), "."))
             ## integrate: 1 or 2.
             if (!(object@integrate %in% c(1L, 2L)))
                 msg <- validMsg(msg, paste0("'integrate' has to be either 1",
                                             " or 2."))
             ## mzdiff: length 1.
             if (length(object@mzdiff) != 1)
                 msg <- validMsg(msg, paste0("'mzdiff' has to be a numeric of",
                                             " length 1."))
             ## noise: length 1.
             if (length(object@noise) != 1)
                 msg <- validMsg(msg, paste0("'noise' has to be a numeric of",
                                             " length 1."))
             ## fitgauss: length 1.
             if (length(object@fitgauss) != 1)
                 msg <- validMsg(msg, paste0("'fitgauss' has to be a numeric of",
                                             " length 1."))
             ## verboseColumns: length 1.
             if (length(object@verboseColumns) != 1)
                 msg <- validMsg(msg, paste0("'verboseColumns' has to be a ",
                                             "numeric of length 1."))
             ## firstBaselineCheck: length 1.
             if (length(object@firstBaselineCheck) != 1)
                 msg <- validMsg(msg, paste0("'firstBaselineCheck' has to be a",
                                             " numeric of length 1."))
             ## roiList: check
             if (length(object@roiList) > 0) {
                 doHaveExpectedEls <- function(z) {
                     need <- c("scmax", "scmin", "mzmin", "mzmax", "length",
                               "intensity")
                     ## Each element should be of length (nrow) 1 and should
                     ## have scmax, scmin, mzmin, mzmax, length, intensity.
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
             ## roiScales: same length then roiList.
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
##' mfp <- MatchedFilterParam(binSize = 0.2)
##' ## Change snthresh parameter
##' snthresh(mfp) <- 15
##' mfp
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
             ## Check the values.
             ## binSize positive numeric of length 1.
             if (length(object@binSize) != 1 | object@binSize < 0)
                 msg <- validMsg(msg, paste0("'binSize' has to be positive",
                                             " numeric of length 1."))
             ## impute
             if (!any(c("none", "lin", "linbase") == object@impute))
                 msg <- validMsg(msg,
                                 paste0("Only values 'none', 'lin' and ",
                                        "'linbase' are allowed for'impute'"))
             ## baseValue
             if (length(object@baseValue) > 1)
                 msg <- validMsg(msg, paste0("'baseValue' has to be a",
                                             " numeric of length 1."))
             ## distance
             if (length(object@distance) > 1)
                 msg <- validMsg(msg, paste0("'distance' has to be a numeric",
                                             " of length 1."))
             ## fwhm
             if (length(object@fwhm) != 1)
                 msg <- validMsg(msg, paste0("'fwhm' has to be a numeric",
                                             " of length 1."))
             ## sigma
             if (length(object@sigma) != 1)
                 msg <- validMsg(msg, paste0("'sigma' has to be a numeric",
                                             " of length 1."))
             ## max
             if (length(object@max) != 1)
                 msg <- validMsg(msg, paste0("'max' has to be a numeric",
                                             " of length 1."))
             ## snthresh
             if (length(object@snthresh) != 1)
                 msg <- validMsg(msg, paste0("'snthresh' has to be a numeric",
                                             " of length 1."))
             ## steps
             if (length(object@steps) != 1)
                 msg <- validMsg(msg, paste0("'steps' has to be a numeric",
                                             " of length 1."))
             ## mzdiff
             if (length(object@mzdiff) != 1)
                 msg <- validMsg(msg, paste0("'mzdiff' has to be a numeric",
                                             " of length 1."))
             ## index
             if (length(object@index) != 1)
                 msg <- validMsg(msg, paste0("'index' has to be a logical",
                                             " of length 1."))
             if (is.null(msg)) {
                 return(TRUE)
             } else {
                 return(msg)
             }
         })
