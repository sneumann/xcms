#' @title Next Generation `xcms` Result Object
#'
#' @aliases XcmsExperiment-class show,XcmsExperiment-method
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
#' - `[`: subset an `XcmsExperiment` by **sample** (parameter `i`). Subsetting
#'   will by default drop correspondence results (as subsetting by samples will
#'   obviously affect the feature definition) and alignment results (adjusted
#'   retention times) while identified chromatographic peaks (for the selected
#'   samples) will be retained. Which preprocessing results should be
#'   kept or dropped can also be configured with optional parameters
#'   `keepChromPeaks` (by default `TRUE`), `keepAdjustedRtime` (by default
#'   `FALSE`) and `keepFeatures` (by default `FALSE`).
#'
#' - `filterFile`: filter an `XcmsExperiment` (or `MsExperiment`) by *file*
#'   (sample). The index of the samples to which the data should be subsetted
#'   can be specified with parameter `file`. The sole purpose of this function
#'   is to provide backward compatibility with the `MSnbase` package. Wherever
#'   possible, the `[` function should be used instead for any sample-based
#'   subsetting.
#'   Note also that in contrast to `[`, `filterFile` does not support subsetting
#'   in arbitrary order.
#'
#' - `filterRt`: filter an `XcmsExperiment` keeping only data within the
#'   specified retention time range (parameter `rt`). This function will keep
#'   all preprocessing results present within the retention time range: all
#'   identified chromatographic peaks with the retention time of the apex
#'   position within the retention time range `rt` are retained.
#'   Parameter `msLevel.` is currently ignored, i.e. filtering will always
#'   performed on **all** MS levels of the object.
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
#' @section Functionality related to alignment:
#'
#' - `adjustRtime`: performs retention time adjustment (alignment) of the data.
#'   See [adjustRtime()] for details.
#'
#' - `hasAdjustedRtime`: whether alignment was performed on the object (i.e.,
#'   the object contains alignment results).
#'
#' @section Functionality for backward compatibility:
#'
#' These functions for `MsExperiment` and `XcmsExperiment` ensure compatibility
#' with the *older* [XCMSnExp()] xcms result object.
#'
#' - `fileNames`: returns the original data file names for the spectra data.
#'   Ideally, the `dataOrigin` or `dataStorage` spectra variables from the
#'   object's `spectra` should be used instead.
#'
#' - `fromFile`: returns the file (sample) index for each spectrum within
#'   `object`. Generally, subsetting by sample using the `[` is the preferred
#'    way to get spectra from a specific sample.
#'
#' - `rtime`: extract retention times of the **spectra** from the
#'   `MsExperiment` or `XcmsExperiment` object. It is thus a shortcut for
#'   `rtime(spectra(object))` which would be the preferred way to extract
#'   retention times from an `MsExperiment`..
#'
#' @section Differences compared to the [XCMSnExp()] object:
#'
#' - Subsetting by `[` supports arbitrary ordering.
#'
#' @param drop For `[`: ignored.
#'
#' @param file For `filterFile`: `integer` with the indices of the samples
#'     (files) to which the data should be subsetted.
#'
#' @param i For `[`: `integer` or `logical` defining the samples/files to
#'     subset.
#'
#' @param isFilledColumn For `chromPeaks`: `logical(1)` whether a column
#'     `"is_filled"` should be included in the returned `matrix` with the
#'     information whether a peak was detected or *only* filled-in. Note that
#'     this information is also provided in the `chromPeakData` data frame.
#'
#' @param j For `[`: not supported.
#'
#' @param keepAdjustedRtime `logical(1)`: whether adjusted retention times (if
#'     present) should be retained.
#'
#' @param msLevel `integer` defining the MS level (or multiple MS level if the
#'     function supports it).
#'
#' @param msLevel. For `filterRt`: ignored. `filterRt` will always filter
#'     by retention times on all MS levels regardless of this parameter.
#'
#' @param mz For `chromPeaks`: `numeric(2)` optionally defining the m/z range
#'     for which chromatographic peaks should be returned. The full m/z range
#'     is used by default.
#'
#' @param ppm For `chromPeaks`: optional `numeric(1)` specifying the ppm by
#'     which the m/z range (defined by `mz` should be extended. For a value of
#'     `ppm = 10`, all peaks within `mz[1] - ppm / 1e6` and
#'     `mz[2] + ppm / 1e6` are returned.
#'
#' @param object An `XcmsExperiment` object.
#'
#' @param rt For `chromPeaks`: `numeric(2)` defining the retention time range
#'     for which chromatographic peaks should be returned. The full range is
#'     used by default.
#'
#' @param type For `chromPeaks` and only if either `mz` and `rt` are defined
#'     too: `character(1)`: defining which peaks should be returned. For
#'     `type = "any"`: returns all chromatographic peaks also only partially
#'     overlapping any of the provided ranges. For `type = "within"`: returns
#'     only peaks completely within the region defined by `mz` and/or `rt`.
#'     For `type = "apex_within"`: returns peaks for which the m/z and retention
#'     time of the peak's apex is within the region defined by `mz` and/or `rt`.
#'
#' @param x An `XcmsExperiment` object.
#'
#' @param ... Additional optional parameters.
#'
#' @name XcmsExperiment
#'
#' @author Johannes Rainer
#'
#' @exportClass XcmsExperiment
#'
#' @md
#'
#' @examples
#'
#' ## Creating a MsExperiment object representing the data from an LC-MS
#' ## experiment.
#' library(MsExperiment)
#' mse <- MsExperiment()
#'
#' ## Defining the raw data files
#' fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko16.CDF', package = "faahKO"),
#'          system.file('cdf/KO/ko18.CDF', package = "faahKO"))
#' fls <- normalizePath(fls)
#'
#' ## Defining a data frame with the sample characterization
#' df <- data.frame(mzML_file = basename(fls),
#'                 dataOrigin = fls,
#'                 sample = c("ko15", "ko16", "ko18"))
#' ## Adding the sample definition to the experiment
#' sampleData(mse) <- DataFrame(df)
#'
#' ## Load the MS data as a Spectra object and add it to the experiment
#' spectra(mse) <- Spectra::Spectra(fls)
#'
#' ## Linking samples to MS spectra. This step is required to define
#' ## which spectra (and eventually files) belong to which sample
#' mse <- linkSampleData(mse, with = "sampleData.dataOrigin = spectra.dataOrigin")
#'
#' ## Perform peak detection on the data using the centWave algorith. Note
#' ## that the parameters are chosen to reduce the run time of the example.
#' p <- CentWaveParam(noise = 10000, snthresh = 40, prefilter = c(3, 10000))
#' xmse <- findChromPeaks(mse, param = p)
#' xmse
#'
#' ## Have a quick look at the identified chromatographic peaks
#' head(chromPeaks(xmse))
#'
#' ## Extract chromatographic peaks identified between 3000 and 3300 seconds
#' chromPeaks(xmse, rt = c(3000, 3300), type = "within")
#'
#' ## Subsetting the data to the results (and data) for the second sample
#' a <- xmse[2]
#' nrow(chromPeaks(xmse))
#' nrow(chromPeaks(a))
#'
#' ## Filtering the result by retention time: keeping all spectra and
#' ## chromatographic peaks within 3000 and 3500 seconds.
#' xmse_sub <- filterRt(xmse, rt = c(3000, 3500))
#' xmse_sub
#' nrow(chromPeaks(xmse_sub))
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
    if (hasAdjustedRtime(object))
        cat(" - adjusted retention times: mean absolute difference",
            format(mean(abs(rtime(spectra(object)) -
                           spectra(object)$rtime_adjusted)),
                  digits = 3), "seconds\n")
})

################################################################################
## Filtering and subsetting
################################################################################

#' @rdname XcmsExperiment
setMethod("[", "XcmsExperiment", function(x, i, j, ...) {
    if (!missing(j))
        stop("subsetting by j not supported")
    .subset_xcms_experiment(x, i = i, ...)
})

#' @rdname XcmsExperiment
setMethod(
    "filterRt", "XcmsExperiment",
    function(object, rt, msLevel.) {
        if (missing(rt))
            return(object)
        rt <- range(rt)
        if (!missing(msLevel.))
            warning("Parameter 'msLevel.' currently ignored.", call. = FALSE)
        msLevel. <- uniqueMsLevels(spectra(object))
        ## Subset chrom peaks
        if (hasChromPeaks(object)) {
            crt <- object@chromPeaks[, "rt"]
            keep <- which(crt >= rt[1] & crt <= rt[2])
            object@chromPeaks <- object@chromPeaks[keep, , drop = FALSE]
        }
        ## for features: keep all features for the chrom peaks that we kept.
        callNextMethod(object = object, rt = rt, msLevel. = msLevel.)
    })

################################################################################
## chromatographic peaks
################################################################################

#' @rdname findChromPeaks
setMethod(
    "findChromPeaks",
    signature(object = "MsExperiment", param = "Param"),
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
setMethod(
    "findChromPeaks",
    signature(object = "XcmsExperiment", param = "Param"),
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
setMethod("chromPeaks", "XcmsExperiment", function(object, rt = numeric(),
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
    DataFrame(object@chromPeakData)
})

################################################################################
## alignment
################################################################################

#' @rdname adjustRtime
setMethod(
    "adjustRtime", signature(object = "MsExperiment", param = "ObiwarpParam"),
    function(object, param, chunkSize = 2L, BPPARAM = bpparam()) {
        msLevel <- 1L
        res <- .mse_obiwarp_chunks(
            object, param = param, chunkSize = chunkSize, BPPARAM = BPPARAM)
        ## Saving adjusted rtimes into $rtime_adjusted
        rt_adj <- rep(NA_real_, length(spectra(object)))
        rt_adj[object@sampleDataLinks[["spectra"]][, 2L]] <-
            unlist(res, use.names = FALSE)
        object@spectra$rtime_adjusted <- rt_adj
        if (!is(object, "XcmsExperiment"))
            object <- as(object, "XcmsExperiment")

        ph <- XProcessHistory(param = param,
                              type. = .PROCSTEP.RTIME.CORRECTION,
                              fileIndex. = seq_along(object),
                              msLevel = msLevel)
        object@processHistory <- c(object@processHistory, list(ph))
        validObject(object)
        object
})

## For XcmsExperiment: need also to adjust chrom peaks.
## applyAdjustedRtime,XcmsExperiment
## dropAdjustedRtime,XcmsExperiment

#' @rdname XcmsExperiment
setMethod("hasAdjustedRtime", "MsExperiment", function(object) {
    any(spectraVariables(spectra(object)) == "rtime_adjusted")
})

################################################################################
## utility and unsorted methods
################################################################################
