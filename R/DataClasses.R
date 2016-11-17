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

##' @title Feature detection using the centWave method
##'
##' @description The \code{CentWaveParam} class allows to specify all settings for
##' a feature detection using the centWave method.
##'
##' @slot ppm,peakwidth,snthresh,prefilter,mzCenterFun,integrate,mzdiff,fitgauss,noise,verboseColumns,roiList,firstBaselinCheck,roiScales See corresponding parameter above.
## @slot peakwidth See corresponding parameter below.
## @slot snthresh See corresponding parameter below.
## @slot prefilter See corresponding parameter below.
## @slot mzCenterFun See corresponding parameter below.
## @slot integrate See corresponding parameter below.
## @slot mzdiff See corresponding parameter below.
## @slot fitgauss See corresponding parameter below.
## @slot noise See corresponding parameter below.
## @slot verboseColumns See corresponding parameter below.
## @slot roiList See corresponding parameter below.
## @slot firstBaselineCheck See corresponding parameter below.
## @slot roiScales See corresponding parameter below.
##' @rdname featureDetection-centWave
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
                 msg <- valisMsg(msg, paste0("'mzdiff' has to be a numeric of",
                                             " length 1."))
             ## noise: length 1.
             if (length(object@noise) != 1)
                 msg <- valisMsg(msg, paste0("'noise' has to be a numeric of",
                                             " length 1."))
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
