#' @include hidden_aliases.R

#' XcmsExperimentHdf5 uses a different convention for chrom peak IDs: for
#' efficiency reasons, chromatographic peak data is organized by MS level and
#' sample/file. The chrom peak IDs are hence in the format
#' *CP<MS level><sample id><chrom peak index>* with <MS level> being the MS
#' level in which the chromatographic peaks were detected and <sample id>
#' the ID of the sample (usually related to the index in the original
#' `MsExperiment` object) and the <chrom peak index> the index
#' of the chromatographic peak in the chrom peak matrix **of that sample** and
#' MS level.
#'
#' @noRd
NULL

setClass("XcmsExperimentHdf5",
         contains = "XcmsExperiment",
         slots = c(hdf5_file = "character",
                   hdf5_mod_count = "integer",
                   sample_id = "character",
                   has_chrom_peaks = "logical",
                   has_features = "logical"),
         prototype = prototype(
             hdf5_file = character(),
             hdf5_mod_count = 0L,
             sample_id = character(),
             has_chrom_peaks = FALSE,
             has_features = FALSE
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
            msl <- .h5_chrom_peak_ms_levels(object@hdf5_file,
                                            object@sample_id[1L])
            cat("  - chromatographic peaks in MS level(s):",
                paste0(msl, collapse = ", "), "\n")
        }
    if (hasAdjustedRtime(object))
        cat("  - adjusted retention times: mean absolute difference",
            format(mean(abs(rtime(spectra(object)) -
                           spectra(object)$rtime_adjusted)),
                   digits = 3), "seconds\n")
    }
    ## if (hasFeatures(object))
    ##     cat("  - correspondence results:", nrow(object@featureDefinitions),
    ##         "features in MS level(s):",
    ##         paste(unique(object@featureDefinitions$ms_level), collapse = ", "),
    ##         "\n")
})

################################################################################
##
##        CHROM PEAKS FUNCTIONALITY
##
################################################################################

#' @rdname hidden_aliases
setMethod("hasChromPeaks", "XcmsExperimentHdf5",
          function(object, msLevel = integer()) {
              .h5_require_rhdf5()
              if (!length(object)) return(FALSE)
              if (object@has_chrom_peaks) {
                  if (length(msLevel)) {
                      msl <- .h5_chrom_peak_ms_levels(object@hdf5_file,
                                                      object@sample_id[1L])
                      all(msLevel %in% msl)
                  } else TRUE
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
             msLevel = integer(), sample = integer(),
             type = c("any", "within", "apex_within"),
             columns = character(), isFilledColumn = FALSE) {
        type <- match.arg(type)
    stop("Not implemented for ", class(object)[1L])
        ## pks <- object@chromPeaks
        ## if (isFilledColumn)
        ##     pks <- cbind(
        ##         pks, is_filled = as.numeric(object@chromPeakData$is_filled))
        ## pks[.index_chrom_peaks(object, rt = rt, mz = mz, ppm = ppm,
        ##                        msLevel = msLevel, type = type), , drop = FALSE]
    })

#' @rdname hidden_aliases
setReplaceMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, value) {
    stop("Not implemented for ", class(object)[1L])
})

#' @rdname hidden_aliases
setMethod(
    "chromPeakData", "XcmsExperimentHdf5",
    function(object, msLevel = integer(), sample = integer(),
             return.type = c("DataFrame", "data.frame")) {
        return.type <- match.arg(return.type)
        stop("Not implemented for ", class(object)[1L])
        ## if (return.type == "DataFrame")
        ##     as(.chromPeakData(object, msLevel = msLevel), "DataFrame")
        ## else .chromPeakData(object, msLevel = msLevel)
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
            stop("Needs to be implemented")
            ## object <- dropFeatureDefinitions(object)
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

#' chromPeaks needs to iterate through all - will need also a chunkSize for
#' that.
