#' @title Next Generation `xcms` Result Object
#'
#' @aliases XcmsExperiment-class
#'
#' @description
#'
#' The `XcmsExperiment` is a data container for `xcms` preprocessing results
#' (i.e. results from chromatographic peak detection, alignment and
#' correspondence analysis).
#'
#' It provides the same functionality than the [XCMSnExp] object, but uses the
#' more advanced and modern MS infrastructure provided by the `MsExperiment`
#' and `Spectra` Bioconductor packages. With this comes a higher flexibility on
#' how and where to store the data.
#'
#' Documentation of the various functions for `XcmsExperiment` objects are
#' grouped by topic and provided in the sections below.
#'
#' @section Subsetting and filtering:
#'
#' @section Functionality related to chromatographic peaks:
#'
#' - `chromPeaks`: returns a `numeric` matrix with the identified
#'   chromatographic peaks. Each row represents a chromatographic peak
#'   identified in one sample (file). The number of columns depends on the
#'   peak detection algorithm (see [findChromPeaks()]) but most methods return
#'   the following columns: `"mz"` (intensity-weighted mean of the m/z values
#'   of all mass peaks included in the chromatographic peak), `"mzmin"` (
#'   smallest m/z value of any mass peak in the chromatographic peak), `"mzmax"`
#'   (largest m/z value of any mass peak in the chromatographic peak), `"rt"`
#'   (retention time of the peak apex), `"rtmin"` (retention time of the first
#'   scan/mass peak of the chromatographic peak), `"rtmax"` (retention time of
#'   the last scan/mass peak of the chromatographic peak), `"into"` (integrated
#'   intensity of the chromatographic peak), `"maxo"` (maximal intensity of any
#'   mass peak of the chromatographic peak), `"sample"` (index of the sample
#'   in `object` in which the peak was identified). Parameters `rt`, `mz`,
#'   `ppm`, `msLevel` and `type` allow to extract subsets of identified
#'   chromatographic peaks from the `object`. See parameter description below
#'   for details.
#'
#' - `chromPeakData`: returns a `DataFrame` with potential additional
#'   *annotations* for the identified chromatographic peaks. Each row in this
#'   `DataFrame` corresponds to a row (same index and row name) in the
#'   `chromPeaks` matrix. The default *annotations* are `"ms_level"` (the MS
#'   level in which the peak was identified) and `"is_filled"` (whether the
#'   chromatographic peak was *detected* (by `findChromPeaks`) or *filled-in*
#'   (by `fillChromPeaks`).
#'
#' - `dropChromPeaks`: removes (all) chromatographic peak detection results
#'   from `object`. This will also remove any correspondence results (i.e.
#'   features) and eventually present adjusted retention times from the object.
#'   Alignment results (adjusted retention times) can be retained if parameter
#'   `keepAdjustedRtime` is set to `TRUE`.
#'
#' - `findChromPeaks`: perform chromatographic peak detection. See
#'   [findChromPeaks()] for details.
#'
#' - `hasChromPeaks`: whether the object contains peak detection results.
#'   Parameter `msLevel` allows to check whether peak detection results are
#'   available for the specified MS level(s).
#'
#' @param keepAdjustedRtime `logical(1)`: whether adjusted retention times (if
#'     present) should be retained.
#'
#' @param msLevel `integer` defining the MS level (or multiple MS level if the
#'     function supports it).
#'
#' @param object An `XcmsExperiment` object.
#'
#' @name XcmsExperiment
#'
#' @author Johannes Rainer
#'
#' @exportClass XcmsExperiment
#'
#' @md
NULL

.empty_chrom_peaks <- function(sample = TRUE) {
    cols <- .REQ_PEAKS_COLS
    if (!sample)
        cols <- cols[cols != "sample"]
    matrix(numeric(), ncol = length(cols), nrow = 0,
           dimnames = list(character(), cols))
}

setClass("XcmsExperiment",
         contains = "MsExperiment",
         slots = c(chromPeaks = "matrix",
                   chromPeakData = "data.frame",
                   processHistory = "list"),
         prototype = prototype(chromPeaks = .empty_chrom_peaks(),
                               chromPeakData = data.frame(
                                   ms_level = integer(),
                                   is_filled = logical()),
                               processHistory = list()))

setValidity("XcmsExperiment", function(object) {
    msg <- c(.mse_valid_chrom_peaks(object@chromPeaks),
             .mse_valid_chrom_peak_data(object@chromPeakData),
             .mse_same_rownames(object@chromPeaks, object@chromPeakData))
    if (is.null(msg)) TRUE
    else msg
})

setMethod("show", "XcmsExperiment", function(object) {
    callNextMethod()
    cat(" xcms results:\n")
    if (hasChromPeaks(object))
        cat(" - chromatographic peaks:", nrow(object@chromPeaks),
            "in MS level(s):",
            paste(unique(object@chromPeakData$ms_level), collapse = ", "), "\n")
})

################################################################################
## Filtering and subsetting
################################################################################



################################################################################
## chromatographic peaks
################################################################################

#' @rdname findChromPeaks
setMethod("findChromPeaks", signature(object = "MsExperiment", param = "Param"),
          function(object, param, msLevel = 1L, chunkSize = 2L, ...,
                   BPPARAM = bpparam()) {
              if (length(msLevel) > 1)
                  stop("Currently only peak detection in a single MS level is ",
                       "supported", call. = FALSE)
              if (chunkSize < 0) {
                  res <- .mse_find_chrom_peaks(
                      object, msLevel = msLevel, param = param,
                      BPPARAM = BPPARAM)
              } else {
                  res <- .mse_find_chrom_peaks_chunks(
                      object, msLevel = msLevel, param = param,
                      chunkSize = chunkSize, BPPARAM = BPPARAM)
              }
              ## Assign/define peak IDs.
              pkd <- data.frame(ms_level = rep(msLevel, nrow(res)),
                                is_filled = rep(FALSE, nrow(res)))
              ph <- XProcessHistory(param = param,
                                    type. = .PROCSTEP.PEAK.DETECTION,
                                    fileIndex. = seq_along(object),
                                    msLevel = msLevel)
              if (!is(object, "XcmsExperiment"))
                  object <- as(object, "XcmsExperiment")
              object@processHistory <- c(object@processHistory, list(ph))
              object <- .mse_add_chrom_peaks(object, res, pkd)
              validObject(object)
              object
          })

#' @rdname findChromPeaks
setMethod("findChromPeaks", signature(object = "XcmsExperiment",
                                      param = "Param"),
          function(object, param, msLevel = 1L, chunkSize = 2L, add = FALSE,
                   ..., BPPARAM = bpparam()) {
              ## if (hasFeatures(object)) {
              ##     message("Remove feature definitions")
              ##     object <- dropFeatureDefinitions(
              ##         object, keepAdjustedRtime = hasAdjustedRtime(object))
              ## }
              if (hasChromPeaks(object) && !add) {
                  message("Remove previously identified chromatographic peaks")
                  object <- dropChromPeaks(
                      object, keepAdjustedRtime = hasAdjustedRtime(object))
              }
              callNextMethod()
          })

#' @rdname XcmsExperiment
setMethod("hasChromPeaks", "XcmsExperiment",
          function(object, msLevel = integer()) {
              if (length(msLevel))
                  any(object@chromPeakData$ms_level %in% msLevel)
              else as.logical(nrow(object@chromPeaks))
})

#' @rdname XcmsExperiment
setMethod(
    "dropChromPeaks", "XcmsExperiment",
    function(object, keepAdjustedRtime = FALSE) {
        if (hasChromPeaks(object)) {
            object@processHistory <- dropProcessHistoriesList(
                object@processHistory,
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT))
            object@chromPeaks <- .empty_chrom_peaks()
            object@chromPeakData <- data.frame(ms_level = integer(),
                                               is_filled = logical())
            ## if (hasAdjustedRtime(object) && !keepAdjustedRtime) {
            ##     object <- dropAdjustedRtime(object)
            ## }
            ## if (hasFeatureDefinitions(object)) {
            ##     object <- dropFeatureDefinitions(object)
            ## }
        }
        object
    })

#' @rdname XcmsExperiment
setMethod("chromPeaks", "XCMSnExp", function(object, rt = numeric(),
                                             mz = numeric(), ppm = 0,
                                             msLevel = integer(),
                                             type = c("any", "within",
                                                      "apex_within"),
                                             isFilledColumn = FALSE) {
    type <- match.arg(type)
    pks <- object@chromPeaks
    if (isFilledColumn)
        pks <- cbind(
            pks, is_filled = as.numeric(object@chromPeakData$is_filled))
    if (length(msLevel))
        pks <- pks[which(object@chromPeakData$ms_level %in% msLevel), ,
                   drop = FALSE]
    ## Select peaks within rt range.
    if (nrow(pks) && length(rt)) {
        rt <- range(as.numeric(rt))
        if (type == "any")
            keep <- which(pks[, "rtmin"] <= rt[2] & pks[, "rtmax"] >= rt[1])
        if (type == "within")
            keep <- which(pks[, "rtmin"] >= rt[1] & pks[, "rtmax"] <= rt[2])
        if (type == "apex_within")
            keep <- which(pks[, "rt"] >= rt[1] & pks[, "rt"] <= rt[2])
        pks <- pks[keep, , drop = FALSE]
    }
    ## Select peaks within mz range, considering also ppm
    if (nrow(pks) && length(mz)) {
        mz <- range(as.numeric(mz))
        ## Increase mz by ppm.
        if (is.finite(mz[1]))
            mz[1] <- mz[1] - mz[1] * ppm / 1e6
        if (is.finite(mz[2]))
            mz[2] <- mz[2] + mz[2] * ppm / 1e6
        if (type == "any")
            keep <- which(pks[, "mzmin"] <= mz[2] & pks[, "mzmax"] >= mz[1])
        if (type == "within")
            keep <- which(pks[, "mzmin"] >= mz[1] & pks[, "mzmax"] <= mz[2])
        if (type == "apex_within")
            keep <- which(pks[, "mz"] >= mz[1] & pks[, "mz"] <= mz[2])
        pks <- pks[keep, , drop = FALSE]
    }
    pks
})

#' @rdname XcmsExperiment
setMethod("chromPeakData", "XcmsExperiment", function(object) {
    object@chromPeakData
})
