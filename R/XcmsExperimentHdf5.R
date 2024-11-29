#' @include hidden_aliases.R

#' @title xcms result object for very large data sets
#'
#' @aliases XcmsExperimentHdf5-class
#'
#' @name XcmsExperimentHdf5
#'
#' @description
#'
#' The *xcms* result objects [XcmsExperiment()] and [XCMSnExp()] keep all
#' preprocessing results in memory and can thus (depending on the size of the
#' data set) require a large amount of memory. The `XcmsExperimentHdf5` class,
#' by using an on-disk data storage mechanism, has a much lower memory
#' footprint allowing thus also to represent and analyze very large data sets
#' on regular computer systems. With some exceptions, including additional
#' parameters, the functionality and usability of this object is identical to
#' the one of `XcmsExperiment` objects. This help page lists only functions
#' that have additional or different parameters than the *default* ones for
#' [XcmsExperiment()] objects.
#'
#' @section Using the HDF5 file-based on-disk data storage:
#'
#' Calling [findChromPeaks()] on an `MsExperiment` using the parameter
#' `hdf5File` will return an instance of the `XcmsExperimentHdf5` class and
#' hence use the on-disk data storage mode described on this page. All data
#' is stored in a file with the name specifyied with parameter `hdf5File`.
#'
#' @details
#'
#' The `XcmsExperimentHdf5` object stores all preprocessing results (except
#' adjusted retention times, which are stored as an additional spectra variable
#' in the object's [Spectra()] object), in a file in HDF5 format.
#'
#' `XcmsExperimentHdf5` uses a different naming scheme for chromatographic
#' peaks: for efficiency reasons, chromatographic peak data is organized by
#' sample and MS level. The chrom peak IDs are hence in the format
#' *CP<MS level>S<sample id><chrom peak index>* with <MS level> being the MS
#' level in which the chromatographic peaks were detected and <sample id>
#' the ID of the sample (usually related to the index in the original
#' `MsExperiment` object) and the <chrom peak index> the index of the
#' chromatographic peak in the chrom peak matrix **of that sample** and
#' MS level.
#'
#' HDF5 does not support parallel processing, thus preprocessing results need
#' to be loaded sequentially.
#'
#' All functionality for `XcmsExperimentHdf5` objects is optimized to reduce
#' memory demand at the cost of eventually lower performance.
#'
#' @section Functionality related to chromatographic peaks:
#'
#' - `chromPeaks()` gains parameter `bySample = FALSE` that, if set to `TRUE`
#'   returns a `list` of `chromPeaks` matrices, one for each sample. Due to
#'   the way data is organized in `XcmsExperimentHdf5` objects this is more
#'   efficient than `bySample = FALSE`. Thus, in cases where chrom peak data
#'   is subsequently evaluated or processed by sample, it is suggested to
#'   use `bySample = TRUE`.
#'
#' - `chromPeakData()` gains a new parameter `peaks = character()` which allows
#'   to specify from which chromatographic peaks data should be returned.
#'   For these chromatographic peaks the ID (row name in `chromPeaks()`)
#'   should be provided with the `peaks` parameter. This can reduce the memory
#'   requirement for cases in which only data of some selected chromatographic
#'   peaks needs to be extracted. Also, `chromPeakData()` supports the
#'   `bySample` parameter described for `chromPeaks()` above.
#'
#' @section Retention time alignment:
#'
#' - `adjustRtimePeakGroups()` and `adjustRtime()` with `PeakGroupsParam`:
#'   parameter `extraPeaks` of `PeakGroupsParam` is **ignored**. Anchor peaks
#'   are thus only defined using the `minFraction` and the optional `subset`
#'   parameter.
#'
#' @section Correspondence analysis results:
#'
#' - `featureDefinitions()`: similarly to `featureDefinitions()` for
#'   [XcmsExperiment] objects, this method returns a `data.frame` with the
#'   characteristics for the defined LC-MS features. The function for
#'   `XcmsExperimentHdf5` does however **not** return the `"peakidx"` column
#'   with the indices of the chromatographic peaks per feature. Also, the
#'   columns are returned in alphabetic order.
#'
#' - `featureValues()`: parameter `value = "index"` (i.e. returning the index
#'   of the chromatographic peaks per feature) is not supported.
#'
#' @author Johannes Rainerr, Philippine Louail
NULL

setClass("XcmsExperimentHdf5",
         contains = "XcmsExperiment",
         slots = c(hdf5_file = "character",
                   hdf5_mod_count = "integer",
                   sample_id = "character",
                   chrom_peaks_ms_level = "integer", # use to keep track of the
                                        # MS levels for which chrom peaks are
                                        # available
                   gap_peaks_ms_level = "integer", # gap-filled chrom peaks
                   features_ms_level = "integer"),
         prototype = prototype(
             hdf5_file = character(),
             hdf5_mod_count = 0L,
             sample_id = character(),
             chrom_peaks_ms_level = integer(),
             gap_peaks_ms_level = integer(),
             features_ms_level = integer()
         ))

setValidity("XcmsExperimentHdf5", function(object) {
    .h5_require_rhdf5()
    if (length(object@hdf5_file)) {
        if(!file.exists(object@hdf5_file))
            return(paste0("Data storage file \"", object@hdf5_file,
                          "\" does not exist!"))
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
    }
    if (length(object@sample_id) != nrow(sampleData(object)))
        return(paste0("Corrupt data: number of samples does not match ",
                      "length of sample IDs."))
    TRUE
})

#' @rdname hidden_aliases
setMethod("show", "XcmsExperimentHdf5", function(object) {
    validObject(object)
    getMethod("show", "MsExperiment")(object)
    cat(" xcms results:\n")
    if (length(object@hdf5_file)) {
        if (hasChromPeaks(object)) {
            cat("  - chromatographic peaks in MS level(s):",
                paste0(object@chrom_peaks_ms_level, collapse = ", "), "\n")
        }
    if (hasAdjustedRtime(object))
        cat("  - adjusted retention times: mean absolute difference",
            format(mean(abs(rtime(spectra(object)) -
                           spectra(object)$rtime_adjusted)),
                   digits = 3), "seconds\n")
    }
    if (hasFeatures(object))
        cat("  - correspondence results in MS level(s):",
            paste(object@features_ms_level, collapse = ", "), "\n")
})

################################################################################
##
##        SUBSET AND FILTER
##
################################################################################

#' @rdname hidden_aliases
setMethod("[", "XcmsExperimentHdf5", function(x, i, j, ...) {
    if (!missing(j))
        stop("subsetting by j not supported")
    .h5_subset_xcms_experiment(x, i = i, ...)
})

################################################################################
##
##        CHROM PEAKS FUNCTIONALITY
##
################################################################################

#' @rdname hidden_aliases
setMethod(
    "findChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "Param"),
    function(object, param, msLevel = 1L, chunkSize = 2L, add = FALSE,
             ..., BPPARAM = bpparam()) {
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
        if (length(msLevel) > 1)
            stop("Currently only peak detection in a single MS level is ",
                 "supported", call. = FALSE)
        if (!msLevel %in% uniqueMsLevels(object)) {
            warning("No spectra of MS level ", msLevel, " present in 'object'",
                    call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Remove feature definitions")
            object <- dropFeatureDefinitions(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        if (hasChromPeaks(object) && !add) {
            message("Remove previously identified chromatographic peaks")
            object <- dropChromPeaks(
                object, keepAdjustedRtime = hasAdjustedRtime(object))
        }
        res <- .mse_spectrapply_chunks(
            object, FUN = .h5_find_chrom_peaks_chunk, msLevel = msLevel,
            param = param, h5_file = object@hdf5_file,
            sample_id = object@sample_id, chunkSize = chunkSize,
            BPPARAM = BPPARAM, add = add & hasChromPeaks(object, msLevel))
        ph <- XProcessHistory(param = param,
                              type. = .PROCSTEP.PEAK.DETECTION,
                              fileIndex. = seq_along(object),
                              msLevel = msLevel)
        object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
        object@chrom_peaks_ms_level <- unique(
            c(object@chrom_peaks_ms_level, msLevel))
        object <- addProcessHistory(object, ph)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "dropChromPeaks", "XcmsExperimentHdf5",
    function(object, keepAdjustedRtime = FALSE) {
        if (hasChromPeaks(object)) {
            pt <- vapply(object@processHistory, processType, character(1))
            object@processHistory <- dropProcessHistoriesList(
                object@processHistory,
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT))
            object@chrom_peaks_ms_level <- integer()
            if (hasAdjustedRtime(object) && !keepAdjustedRtime) {
                stop("Needs to be implemented")
                ## remove if alignment performed AFTER chrom peaks
                nom <- length(pt) + 1L
                idx_cp <- .match_last(.PROCSTEP.PEAK.DETECTION, pt,
                                      nomatch = nom)
                idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, pt,
                                      nomatch = nom)
                if (idx_al > idx_cp)
                    object <- dropAdjustedRtime(object)
            }
            if (hasFeatures(object)) {
                object <- dropFeatureDefinitions(
                    object, keepAdjustedRtime = keepAdjustedRtime)
            }
        }
        object
    })

#' @rdname hidden_aliases
setMethod("hasChromPeaks", "XcmsExperimentHdf5",
          function(object, msLevel = integer()) {
              .h5_require_rhdf5()
              if (!length(object)) return(FALSE)
              if (length(object@chrom_peaks_ms_level)) {
                  if (length(msLevel))
                      all(msLevel %in% object@chrom_peaks_ms_level)
                  else TRUE
              } else FALSE
          })

#' @rdname hidden_aliases
setReplaceMethod("chromPeaks", "XcmsExperimentHdf5", function(object, value) {
    stop("Not implemented for ", class(object)[1L])
})

#' @rdname hidden_aliases
setMethod(
    "chromPeaks", "XcmsExperimentHdf5",
    function(object, rt = numeric(), mz = numeric(), ppm = 0,
             msLevel = integer(),  type = c("any", "within", "apex_within"),
             isFilledColumn = FALSE, columns = character(), bySample = FALSE) {
        if (isFilledColumn)
            stop("Parameter 'isFilledColumn = TRUE' is not supported")
        type <- match.arg(type)
        .h5_require_rhdf5()
        if (!length(object))
            return(object@chromPeaks)
        .h5_check_mod_count(object@hdf5_file, object@hdf5_mod_count)
        if (!hasChromPeaks(object, msLevel = msLevel))
            return(object@chromPeaks)
        msl <- object@chrom_peaks_ms_level
        if (length(msLevel))
            msl <- msl[msl %in% msLevel]
        .h5_chrom_peaks(object, msLevel = msl, columns = columns,
                        by_sample = bySample, mz = mz, rt = rt, ppm = ppm,
                        type = type)
    })

#' @rdname hidden_aliases
setReplaceMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, value) {
    stop("Not implemented for ", class(object)[1L])
})

#' @rdname XcmsExperimentHdf5
setMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, msLevel = integer(), peaks = character(),
             return.type = c("DataFrame", "data.frame"), bySample = FALSE) {
        return.type <- match.arg(return.type)
        .h5_require_rhdf5()
        if (!length(object))
            return(as(object@chromPeakData, return.type))
        .h5_check_mod_count(object@hdf5_file, object@hdf5_mod_count)
        if (!length(msLevel))
            msLevel <- object@chrom_peaks_ms_level
        if (!hasChromPeaks(object, msLevel = msLevel))
            return(as(object@chromPeakData, return.type))
        if (return.type == "DataFrame")
            as(.h5_chrom_peak_data(object, msLevel, peaks = peaks,
                                   by_sample = FALSE), "DataFrame")
        else .h5_chrom_peak_data(object, msLevel, peaks = peaks,
                                 by_sample = bySample)
    })

## #' @rdname refineChromPeaks
## setMethod(
##     "refineChromPeaks",
##     signature(object = "XcmsExperiment", param = "CleanPeaksParam"),
##     function(object, param = CleanPeaksParam(), msLevel = 1L) {
##         if (!hasChromPeaks(object, msLevel = msLevel)) {
##             warning("No chromatographic peaks for MS level ",
##                     msLevel, " present", call. = FALSE)
##             return(object)
##         }
##         if (hasFeatures(object)) {
##             message("Removing feature definitions")
##             object <- dropFeatureDefinitions(object)
##         }
##         validObject(param)
##         rtw <- .chromPeaks(object)[, "rtmax"] - .chromPeaks(object)[, "rtmin"]
##         keep_ms <- object@chromPeakData$ms_level %in% msLevel
##         keep_rt <- rtw < param@maxPeakwidth & keep_ms
##         keep <- which(keep_rt | !keep_ms)
##         message("Removed ", nrow(.chromPeaks(object)) - length(keep), " of ",
##                 nrow(.chromPeaks(object)), " chromatographic peaks.")
##         object@chromPeaks <- object@chromPeaks[keep, , drop = FALSE]
##         object@chromPeakData <- object@chromPeakData[keep, , drop = FALSE]
##         xph <- XProcessHistory(param = param, date. = date(),
##                                type. = .PROCSTEP.PEAK.REFINEMENT,
##                                fileIndex = seq_along(object),
##                                msLevel = msLevel)
##         object <- addProcessHistory(object, xph)
##         validObject(object)
##         object
##     })

#' @rdname hidden_aliases
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperimentHdf5",
              param = "MergeNeighboringPeaksParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        .h5_require_rhdf5()
        if (!hasChromPeaks(object, msLevel = msLevel)) {
            warning("No chromatographic peaks for MS level ",
                    msLevel, " present", call. = FALSE)
            return(object)
        }
        if (hasFeatures(object)) {
            message("Removing feature definitions")
            object <- dropFeatureDefinitions(object)
        }
        validObject(param)
        .xmse_apply_chunks(
            object, FUN = .h5_xmse_merge_neighboring_peaks,
            msLevel = msLevel, expandRt = param@expandRt,
            expandMz = param@expandMz, ppm = param@ppm,
            minProp = param@minProp, BPPARAM = BPPARAM,
            chunkSize = chunkSize, SUBSET_FUN = .h5_subset_xcms_experiment,
            keepAdjustedRtime = TRUE, ignoreHistory = TRUE)
        ## Update the @hdf5_mod_count with the one from the file.
        object@hdf5_mod_count <- rhdf5::h5read(object@hdf5_file,
                                               "/header/modcount")[1L]
        xph <- XProcessHistory(param = param, date. = date(),
                               type. = .PROCSTEP.PEAK.REFINEMENT,
                               fileIndex = seq_along(object),
                               msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

## #' @rdname refineChromPeaks
## setMethod(
##     "refineChromPeaks",
##     signature(object = "XcmsExperiment", param = "FilterIntensityParam"),
##     function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
##         if (!hasChromPeaks(object, msLevel = msLevel)) {
##             warning("No chromatographic peaks for MS level ",
##                     msLevel, " present", call. = FALSE)
##             return(object)
##         }
##         if (hasFeatures(object)) {
##             message("Removing feature definitions")
##             object <- dropFeatureDefinitions(object)
##         }
##         npks_orig <- nrow(.chromPeaks(object))
##         validObject(param)
##         if (param@nValues == 1L) {
##             if (!any(colnames(.chromPeaks(object)) == param@value))
##                 stop("Column '", param@value, "' not available.")
##             keep <- .chromPeaks(object)[, param@value] >= param@threshold |
##                 .chromPeakData(object)$ms_level != msLevel
##         } else
##             keep <- unlist(.xmse_apply_chunks(
##                 object, .xmse_filter_peaks_intensities,nValues = param@nValues,
##                 threshold = param@threshold, msLevel = msLevel,
##                 keepAdjustedRtime = TRUE, ignoreHistory = TRUE,
##                 BPPARAM = BPPARAM, chunkSize = chunkSize), use.names = FALSE)
##         object@chromPeaks <- object@chromPeaks[keep, , drop = FALSE]
##         object@chromPeakData <- object@chromPeakData[keep, ]
##         message("Reduced from ", npks_orig, " to ", nrow(.chromPeaks(object)),
##                 " chromatographic peaks.")
##         xph <- XProcessHistory(param = param, date. = date(),
##                                type. = .PROCSTEP.PEAK.REFINEMENT,
##                                fileIndex = seq_along(object),
##                                msLevel = msLevel)
##         object <- addProcessHistory(object, xph)
##         validObject(object)
##         object
##     })

#' @rdname hidden_aliases
setMethod("hasFilledChromPeaks", "XcmsExperimentHdf5", function(object) {
    length(object@gap_peaks_ms_level) > 0
})

#' @rdname hidden_aliases
setMethod(
    "fillChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "ChromPeakAreaParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
        if (length(msLevel) != 1)
            stop("Can only perform peak filling for one MS level at a time.")
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level ", msLevel, " present.")
        .h5_valid_file(object@hdf5_file, object@hdf5_mod_count)
        ## Get integration function and other info.
        ph <- .xmse_process_history(object, .PROCSTEP.PEAK.DETECTION,
                                    msLevel = msLevel)
        fill_fun <- .history2fill_fun(ph)
        mzf <- "wMean"
        if (length(ph) && inherits(ph[[1L]], "XProcessHistory")) {
            prm <- ph[[1L]]@param
            if (any(slotNames(prm) == "mzCenterFun"))
                mzf <- prm@mzCenterFun
        } else
            prm <- MatchedFilterParam()
        mzf <- paste0("mzCenter.", gsub("mzCenter.", "", mzf, fixed = TRUE))
        ## Identify for each feature the samples in which there is a missing
        ## value
        fvals <- is.na(featureValues(object, msLevel = msLevel, method = "sum"))
        fvals <- fvals[which(rowSums(fvals) > 0), , drop = FALSE]
        if (nrow(fvals)) {
            message("Defining MS area to integrate signal from")
            fr <- .h5_features_ms_region(
                object, mzmin = param@mzmin, mzmax = param@mzmax,
                rtmin = param@rtmin, rtmax = param@rtmax,
                features = rownames(fvals), ms_level = msLevel)
            ## define the features to integrate signal for each sample
            frl <- lapply(seq_len(ncol(fvals)), function(i) {
                fr[which(fvals[, i]), , drop = FALSE]
            })
            names(frl) <- seq_along(frl)
            idx <- seq_along(object)
            chunks <- split(idx, ceiling(idx / chunkSize))
            message("Integrating signal from raw data files")
            pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                                   "total (:percent) in ",
                                                   ":elapsed"),
                                   total = length(chunks), clear = FALSE)
            pb$tick(0)
            res <- lapply(chunks, function(z, ...) {
                mc <- .h5_xmse_integrate_chrom_peaks(
                    .h5_subset_xcms_experiment(
                        object, i = z, keepChromPeaks = TRUE,
                        keepAdjustedRtime = TRUE, keepFeatures = TRUE,
                        ignoreHistory = TRUE),
                    pal = frl[z], intFun = fill_fun, update_features = TRUE,
                    mzCenterFun = mzf, param = prm, msLevel = msLevel,
                    BPPARAM = BPPARAM)
                pb$tick()
                mc
            })
            object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
            object@gap_peaks_ms_level <- union(object@gap_peaks_ms_level,
                                               msLevel)
            ph <- XProcessHistory(
                param = param, date. = date(), type. = .PROCSTEP.PEAK.FILLING,
                fileIndex = seq_along(object), msLevel = msLevel)
            object <- addProcessHistory(object, ph)
        }
        object
    })

#' @rdname hidden_aliases
setMethod("dropFilledChromPeaks", "XcmsExperimentHdf5", function(object) {
    if (!.hasFilledPeaks(object))
        return(object)
    mc <- object@hdf5_mod_count
    for (msl in object@gap_peaks_ms_level) {
        for (id in object@sample_id) {
            pks <- .h5_read_data(
                object@hdf5_file, id, "chrom_peaks", msl,
                read_colnames = TRUE, read_rownames = TRUE)[[1L]]
            pkd <- .h5_read_data(
                object@hdf5_file, id, "chrom_peak_data", msl,
                read_colnames = TRUE)[[1L]]
            keep <- which(!pkd$is_filled)
            l <- list(pks[keep, , drop = FALSE])
            names(l) <- id
            .h5_write_data(
                object@hdf5_file, l, name = "chrom_peaks", ms_level = msl,
                replace = TRUE, write_colnames = TRUE, write_rownames = TRUE)
            l <- list(extractROWS(pkd, keep))
            names(l) <- id
            .h5_write_data(
                object@hdf5_file, l, name = "chrom_peak_data", ms_level = msl,
                replace = TRUE, write_rownames = FALSE)
            fmap <- .h5_read_data(
                object@hdf5_file, id, "feature_to_chrom_peaks", msl)[[1L]]
            l <- list(fmap[fmap[, 2L] %in% keep, , drop = FALSE])
            names(l) <- id
            mc <- .h5_write_data(
                object@hdf5_file, l, name = "feature_to_chrom_peaks",
                write_colnames = FALSE, write_rownames = FALSE,
                ms_level = msl)
       }
    }
    object@gap_peaks_ms_level <- integer()
    object@hdf5_mod_count <- mc
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.PEAK.FILLING)
                type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                         .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                         .PROCSTEP.PEAK.REFINEMENT)
    object
})

#' TODO: LLLLLL
#' - filterMsLevel
#' - filterRt
#' - filterMz
#'
#' @noRd
NULL

################################################################################
##
##        RETENTION TIME ALIGNMENT
##
################################################################################

#' @rdname XcmsExperimentHdf5
setMethod(
    "adjustRtimePeakGroups", c("XcmsExperimentHdf5", "PeakGroupsParam"),
    function(object, param = PeakGroupsParam(), msLevel = 1L) {
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No features present. Please run 'groupChromPeaks' first.")
        subs <- param@subset
        if (!length(subs)) subs <- seq_along(object)
        if (!all(subs %in% seq_along(object)))
            stop("Parameter 'subset' is out of bounds.")
        object@sample_id <- object@sample_id[subs] # quick subset hack.
        rts <- .h5_feature_values_ms_level(msLevel, object, method = "maxint",
                                           intensity = "into", value = "rt")
        pres <- apply(rts, 1, function(z) sum(!is.na(z)))
        rts <- rts[pres >= param@minFraction * length(subs), , drop = FALSE]
        colnames(rts) <- basename(fileNames(object))[subs]
        rts[order(rowMedians(rts, na.rm = TRUE)), , drop = FALSE]
    })

#' @rdname hidden_aliases
setMethod(
    "adjustRtime", signature(object = "XcmsExperimentHdf5",
                             param = "PeakGroupsParam"),
    function(object, param, msLevel = 1L, ...) {
        if (hasAdjustedRtime(object))
            stop("Alignment results already present. Please either remove ",
                 "them with 'dropAdjustedRtime()' in order to perform an ",
                 "alternative, new, alignment, or use 'applyAdjustedRtime()'",
                 " prior 'adjustRtime()' to perform a second round of ",
                 "alignment.")
        if (any(msLevel != 1L))
            stop("Alignment is currently only supported for MS level 1")
        if (!nrow(peakGroupsMatrix(param))) {
            if (!hasFeatures(object, msLevel = msLevel))
                stop("No feature definitions present in 'object'. Please ",
                     "perform first a correspondence analysis using ",
                     "'groupChromPeaks()'")
            peakGroupsMatrix(param) <- adjustRtimePeakGroups(object, param)
        }
        fidx <- as.factor(fromFile(object))
        rt_raw <- split(rtime(object), fidx)
        rt_adj <- .adjustRtime_peakGroupsMatrix(
            rt_raw, peakGroupsMatrix(param), smooth = smooth(param),
            span = span(param), family = family(param),
            subset = subset(param), subsetAdjust = subsetAdjust(param))
        pt <- vapply(object@processHistory, processType, character(1))
        idx_pg <- .match_last(.PROCSTEP.PEAK.GROUPING, pt, nomatch = -1L)
        if (idx_pg > 0)
            ph <- object@processHistory[idx_pg]
        else ph <- list()
        object <- dropFeatureDefinitions(object)
        object@spectra$rtime_adjusted <- unlist(rt_adj, use.names = FALSE)
        res <- mapply(
            FUN = .h5_update_rt_chrom_peaks_sample,
            object@sample_id, rt_raw, rt_adj,
            MoreArgs = list(ms_level = object@chrom_peaks_ms_level,
                            hdf5_file = object@hdf5_file))
        object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
        xph <- XProcessHistory(
            param = param, type. = .PROCSTEP.RTIME.CORRECTION,
            fileIndex = seq_along(object), msLevel = msLevel)
        object@processHistory <- c(object@processHistory, ph, list(xph))
        validObject(object)
        object
})

#' @rdname hidden_aliases
setMethod("dropAdjustedRtime", "XcmsExperimentHdf5", function(object) {
    if (!hasAdjustedRtime(object))
        return(object)
    ptype <- vapply(object@processHistory, processType, character(1))
    nom <- length(ptype) + 1L
    idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
    idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
    if (hasChromPeaks(object)) {
        fidx <- as.factor(fromFile(object))
        res <- mapply(
            FUN = .h5_update_rt_chrom_peaks_sample,
            object@sample_id,
            split(rtime(object, adjusted = TRUE), fidx),
            split(rtime(object, adjusted = FALSE), fidx),
            MoreArgs = list(ms_level = object@chrom_peaks_ms_level,
                            hdf5_file = object@hdf5_file))
        object@hdf5_mod_count <- max(unlist(res, use.names = FALSE))
    }
    svs <- unique(c(spectraVariables(object@spectra), "mz", "intensity"))
    object@spectra <- selectSpectraVariables(
        object@spectra, svs[svs != "rtime_adjusted"])
    object@processHistory <- dropProcessHistoriesList(
        object@processHistory, type = .PROCSTEP.RTIME.CORRECTION, num = 1L)
    if (hasFeatures(object) && idx_co > idx_al) {
        warning("Had to remove feature definitions along with the adjusted ",
                "retention times because of the dependency between them.")
        object <- dropFeatureDefinitions(object)
    }
    object
})

################################################################################
##
##        CORRESPONDENCE
##
################################################################################

#' @rdname hidden_aliases
setMethod("hasFeatures", "XcmsExperimentHdf5", function(object,
                                                        msLevel = integer()) {
    .h5_require_rhdf5()
    if (!length(object)) return(FALSE)
    if (length(object@features_ms_level)) {
        if (length(msLevel))
            all(msLevel %in% object@features_ms_level)
        else TRUE
    } else FALSE
})

#' @rdname hidden_aliases
setMethod(
    "dropFeatureDefinitions", "XcmsExperimentHdf5",
    function(object, keepAdjustedRtime = FALSE) {
        if (!hasFeatures(object))
            return(object)
        ptype <- vapply(object@processHistory, processType, character(1))
        nom <- length(ptype) + 1L
        idx_al <- .match_last(.PROCSTEP.RTIME.CORRECTION, ptype, nomatch = nom)
        idx_co <- .match_last(.PROCSTEP.PEAK.GROUPING, ptype, nomatch = nom)
        object@processHistory <- dropProcessHistoriesList(
            object@processHistory, type = .PROCSTEP.PEAK.GROUPING, num = 1L)
        object@features_ms_level <- integer()
        if (.hasFilledPeaks(object)) {
            stop("Removal of gap-filled peaks needs to be implemented")
            ## object <- .filter_chrom_peaks(
            ##     object, which(!.chromPeakData(object)$is_filled))
            ## object@processHistory <- dropProcessHistoriesList(
            ##     object@processHistory, type = .PROCSTEP.PEAK.FILLING)
        }
        if (!keepAdjustedRtime && hasAdjustedRtime(object) && idx_al > idx_co) {
            object <- dropAdjustedRtime(object)
        }
        object
    })

#' @rdname hidden_aliases
setMethod(
    "groupChromPeaks",
    signature(object = "XcmsExperimentHdf5", param = "Param"),
    function(object, param, msLevel = 1L, add = FALSE) {
        msLevel <- unique(msLevel)
        if (length(msLevel) != 1)
            stop("Can only perform the correspondence analysis on one MS",
                 " level at a time. Please repeat for other MS levels ",
                 "with parameter `add = TRUE`.")
        if (!hasChromPeaks(object, msLevel))
            stop("No chromatographic peak for MS level ", msLevel,
                 " present. Please perform first a peak detection ",
                 "using the 'findChromPeaks' method.", call. = FALSE)
        if (hasFeatures(object) && !add)
            object <- dropFeatureDefinitions(object)
        if (hasFeatures(object, msLevel) && add)
            stop("Adding new features to existing correspondence results of ",
                 "the same MS level is currently not supported.", call. = FALSE)
        cps <- .h5_chrom_peaks(object, msLevel = msLevel,
                               columns = c("mz", "rt"), by_sample = FALSE)
        res <- .xmse_group_cpeaks(cps, param = param)
        if (!nrow(res))
            return(object)
        object@features_ms_level <- as.integer(
            unique(c(object@features_ms_level, msLevel)))
        cpk_idx <- res$peakidx
        res$peakidx <- NULL
        attr(res, "row.names") <- .featureIDs(
            nrow(res), prefix = paste0("FT", msLevel), min_len = 6)
        ## Save features to "/features/ms_<msLevel>/feature_definitions"
        .h5_write_data(object@hdf5_file, list(features = res),
                       name = "feature_definitions", ms_level = msLevel,
                       replace = TRUE, write_rownames = TRUE)
        ## Define the mapping of features to chrom peaks PER SAMPLE
        i <- ncol(cps)
        cps <- unname(cps)
        sample_f <- factor(
            unlist(lapply(cpk_idx, function(idx) cps[idx, i])),
            levels = seq_along(object@sample_id))
        fts_to_cps <- split.data.frame(
            cbind(rep(seq_len(nrow(res)), lengths(cpk_idx)),
                  unlist(cpk_idx, use.names = FALSE)), f = sample_f)
        ## Fix indices relative to chrom peak matrices WITHIN sample
        ## Maybe add an index -1 for features with no assigned peaks? Just
        ## to ensure all features are listed at least once?
        npks <- as.integer(
            table(factor(cps[, i], levels = seq_along(object@sample_id))))
        fts_to_cps <- mapply(function(a, b) {
            if (nrow(a))
                a[, 2L] <- a[, 2L] - b
            a
        }, fts_to_cps, as.list(c(0, cumsum(npks)[-length(npks)])))
        names(fts_to_cps) <- object@sample_id
        mc <- .h5_write_data(
            object@hdf5_file, fts_to_cps, name = "feature_to_chrom_peaks",
            ms_level = rep(msLevel, length(fts_to_cps)), replace = TRUE,
            write_rownames = FALSE, write_colnames = FALSE)
        object@hdf5_mod_count <- mc
        xph <- XProcessHistory(param = param, type. = .PROCSTEP.PEAK.GROUPING,
                               fileIndex = seq_along(object), msLevel = msLevel)
        object <- addProcessHistory(object, xph)
        validObject(object)
        object
    })

#' @rdname hidden_aliases
setMethod(
    "featureArea", "XcmsExperimentHdf5",
    function(object, mzmin = min, mzmax = max, rtmin = min,
             rtmax = max, features = character(), msLevel = 1L) {
        if (!hasFeatures(object, msLevel))
            stop("No correspondence results available. Please run ",
                 "'groupChromPeaks' first.", call. = FALSE)
        if (!length(features))
            features <- .h5_feature_definitions_rownames(object, msLevel)[[1L]]
        .h5_features_ms_region(
            object, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
            rtmax = rtmax, features, ms_level = msLevel)
    })

#' @rdname hidden_aliases
setReplaceMethod("featureDefinitions", "XcmsExperimentHdf5",
                 function(object, value) {
                     stop("Not implemented for ", class(object)[1L])
                 })

#' @rdname hidden_aliases
setMethod(
    "featureDefinitions", "XcmsExperimentHdf5",
    function(object, mz = numeric(), rt = numeric(), ppm = 0,
             type = c("any", "within", "apex_within"), msLevel = integer()) {
        if (!hasFeatures(object))
            return(object@featureDefinitions)
        type <- match.arg(type)
        if (length(msLevel))
            msLevel <- intersect(msLevel, object@features_ms_level)
        else msLevel <- object@features_ms_level
        if (!length(msLevel))
            return(object@featureDefinitions)
        fd <- .h5_read_data(object@hdf5_file, rep("features", length(msLevel)),
                            name = "feature_definitions", ms_level = msLevel,
                            read_rownames = TRUE)
        msl <- rep(msLevel, vapply(fd, nrow, 1L))
        fd <- do.call(rbind, fd)
        fd$ms_level <- msl
        .subset_feature_definitions(fd, mz = mz, rt = rt,
                                    ppm = ppm, type = type)
    })

#' @rdname hidden_aliases
setMethod(
    "featureValues", "XcmsExperimentHdf5",
    function(object, method = c("medret", "maxint", "sum"), value = "into",
             intensity = "into", filled = TRUE, missing = NA_real_,
             msLevel = integer()) {
        if (!length(msLevel)) msLevel <- object@features_ms_level
        if (!hasFeatures(object, msLevel = msLevel))
            stop("No feature definitions for MS level(s) ", msLevel," present.")
        method <- match.arg(method)
        if (value == "index")
            stop("'featureValues' for 'XcmsExperimentHdf5' does not support ",
                 "'value = \"index\"'.", call. = FALSE)
        if (is.character(missing) && !(missing %in% c("rowmin_half")))
            stop("if 'missing' is not 'NA' or a numeric it should",
                 " be one of: \"rowmin_half\".")
        msLevel <- intersect(msLevel, object@features_ms_level)
        vals <- do.call(
            rbindFill,
            lapply(msLevel, .h5_feature_values_ms_level, x = object,
                   method = method, value = value, intensity = intensity,
                   filled = filled)
        )
        colnames(vals) <- basename(fileNames(object))
        missing <- missing[1L]
        if (!is.na(missing)) {
            if (is.numeric(missing))
                vals[is.na(vals)] <- missing
            if (missing == "rowmin_half")
                for (i in seq_len(nrow(vals))) {
                    nas <- is.na(vals[i, ])
                    if (any(nas))
                        vals[i, nas] <- min(vals[i, ], na.rm = TRUE) / 2
            }
        }
        vals
    })

################################################################################
##
##        OTHER FUNCTIONALIY
##
################################################################################

#' @rdname hidden_aliases
setMethod(
    "chromatogram", "XcmsExperimentHdf5",
    function(object, rt = matrix(nrow = 0, ncol = 2),
             mz = matrix(nrow = 0, ncol = 2), aggregationFun = "sum",
             msLevel = 1L, chunkSize = 2L, isolationWindowTargetMz = NULL,
             return.type = c("XChromatograms", "MChromatograms"),
             include = character(),
             chromPeaks = c("apex_within", "any", "none"),
             BPPARAM = bpparam()) {
        if (!is.matrix(rt)) rt <- matrix(rt, ncol = 2L)
        if (!is.matrix(mz)) mz <- matrix(mz, ncol = 2L)
        if (length(include)) {
            warning("Parameter 'include' is deprecated, please use ",
                    "'chromPeaks' instead")
            chromPeaks <- include
        }
        if (nrow(mz) && !nrow(rt))
            rt <- cbind(rep(-Inf, nrow(mz)), rep(Inf, nrow(mz)))
        if (nrow(rt) && !nrow(mz))
            mz <- cbind(rep(-Inf, nrow(rt)), rep(Inf, nrow(rt)))
        return.type <- match.arg(return.type)
        chromPeaks <- match.arg(chromPeaks)
        if (!hasChromPeaks(object, msLevel))
            chromPeaks <- "none"
        if (hasAdjustedRtime(object))
            object <- applyAdjustedRtime(object)
        .h5_x_chromatograms(
            object, ms_level = msLevel, chromPeaks = chromPeaks,
            mz = mz, rt = rt, aggregationFun = aggregationFun,
            chunkSize = chunkSize, return.type = return.type,
            isolationWindow = isolationWindowTargetMz,
            BPPARAM = BPPARAM)
    })
