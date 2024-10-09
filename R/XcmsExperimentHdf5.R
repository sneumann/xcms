#' @include hidden_aliases.R

setClass("XcmsExperimentHdf5",
         contains = "XcmsExperiment",
         slots = c(hdf5_file = "character",
                   hdf5_mod_count = "integer",
                   sample_id = "integer",
                   has_chrom_peaks = "logical",
                   has_features = "logical"),
         prototype = prototype(
             hdf5_file = character(),
             hdf5_mod_count = 0L,
             sample_id = integer(),
             has_chrom_peaks = FALSE,
             has_features = FALSE
         ))

setValidity("XcmsExperimentHdf5", function(object) {
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

################################################################################
##
##        CHROM PEAKS FUNCTIONALITY
##
################################################################################

#' @rdname hidden_aliases
setMethod("hasChromPeaks", "XcmsExperimentHdf5",
          function(object, msLevel = integer()) {
              if (!length(object)) return(FALSE)
              h5 <- rhdf5::H5Fopen(object@hdf5_file)
              on.exit(rhdf5::H5Fclose(h5))
              msl <- .h5_ms_levels(h5)
              if (!length(msl)) return(FALSE)
              has_cp <- vapply(
                  paste0("/ms_", msl, "/", object@sample_id[1L], "/"),
                  function(x) {
                      any(.h5_dataset_names(x, h5) == "chrom_peaks")
                  }, NA)
              names(has_cp) <- msl
              if (length(msLevel))
                  all(has_cp[as.character(msLevel)] %in% TRUE)
              else any(has_cp)
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

#' @rdname refineChromPeaks
setMethod(
    "refineChromPeaks",
    signature(object = "XcmsExperimentHdf5",
              param = "MergeNeighboringPeaksParam"),
    function(object, param, msLevel = 1L, chunkSize = 2L, BPPARAM = bpparam()) {
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
        ## In chunks of data:
        ## run peak refinement, writing the results back to HDF5 and report
        ## the last chrom peak ID. Use these to update the rownames in all
        ## tables of the same MS level.
        res <- .xmse_apply_chunks(
            ## LLLLL can we fit/adapt xmse_apply_chunks?
            object, .h5_xmse_merge_neighboring_peaks, msLevel = msLevel,
            expandRt = param@expandRt, expandMz = param@expandMz,
            ppm = param@ppm, minProp = param@minProp, BPPARAM = BPPARAM,
            keepAdjustedRtime = TRUE, ignoreHistory = TRUE,
            keepSampleIndex = FALSE, chunkSize = chunkSize)
        ## Update the rownames of all data sets of that MS level.
        ## res should be the highest number per subset. use the max of that to
        ## define the names.

        ## Update the @hdf5_mod_count with the one from the file.


        pks <- do.call(rbind, lapply(res, `[[`, 1L))
        pkd <- do.call(rbind.data.frame, c(lapply(res, `[[`, 2L),
                                           make.row.names = FALSE))
        npks <- unlist(lapply(res, `[[`, 3L), use.names = FALSE)
        pks[, "sample"] <- rep(seq_along(npks), npks)
        nas <- which(is.na(rownames(pks))) # merged peaks
        if (!any(colnames(pkd) == "merged"))
            pkd$merged <- FALSE
        pkd$merged[nas] <- TRUE
        ## Fix rownames AAAAAAHHHHHH!
        maxi <- max(as.integer(sub("CP", "", rownames(object@chromPeaks))))
        rownames(pks)[nas] <- .featureIDs(length(nas), "CP", from = maxi + 1L)
        rownames(pkd) <- rownames(pks)
        ## Merge with existing peaks from **other** MS levels
        keep <- object@chromPeakData$ms_level != msLevel
        if (any(keep)) {
            object@chromPeaks <- rbind(object@chromPeaks[keep, ], pks)
            object@chromPeakData <- rbindFill(object@chromPeakData[keep, ], pkd)
        } else {
            object@chromPeaks <- pks
            object@chromPeakData <- pkd
        }
        message("Reduced from ", npks_orig, " to ", nrow(.chromPeaks(object)),
                " chromatographic peaks.")
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


## Subsetting should be fairly easy, because we just need to subset the
## sample_id vector.
## .h5_subset_xcms_experiment <- function(x, i = integer(),
##                                        keepChromPeaks = TRUE,
##                                        keepAdjustedRtime = FALSE,
##                                        keepFeatures = FALSE,
##                                        ignoreHistory = FALSE,
##                                        keepSampleIndex = FALSE,
##                                        ...) {
##     i <- i2index(i, length(x))
##     if (any(i < 0)) {
##         if (all(i < 0))
##             i <- seq_along(x)[i]
##         else stop("Mixing positive and negative indices is not supported.")
##     }
##     stop("Needs to be implemented")
## }

#' chromPeaks needs to iterate through all - will need also a chunkSize for
#' that.
