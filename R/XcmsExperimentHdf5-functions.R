#' Properties of the HDF5 file used for on-disk storage of xcms results:
#' - all preprocessing results are stored within the same file (in HDF5 format).
#' - the file contains an additional field/data set */header/modcount* that is
#'   used to keep track of every data write operation to the file using an
#'   incremental number. Comparing the value of this *modcount* with the one
#'   stored within the `XcmsExperimentHdf5` object in R can be used to validate
#'   the object/data.
#' - storage of chrom peak detection results are organized by sample and
#'   MS level:
#'   /<sample id>/ms_<ms_level>/chrom_peaks (float array)
#'   /<sample id>/ms_<ms_level>/chrom_peaks_rownames (character array)
#'   /<sample id>/ms_<ms_level>/chrom_peaks_colnames (character array)
#'   /<sample id>/ms_<ms_level>/chrom_peak_data (list of arrays).
#' - the feature definitions `data.frame` is stored as a data set with its
#'   rownames as additional (character) array:
#'   /features/feature_definitions (list of arrays)
#'   /features/feature_definitions_rownames (character array)
#' - the information which chrom peaks of a sample are assigned to which
#'   feature is saved along with the chrom peaks data as a two column integer
#'   array, the first row with the indices of the features, the second with
#'   the index of the chrom peak(s) assigned to the respective feature.
#'   /<sample id>/ms_<ms_level>/features_to_chrom_peaks (integer array)
#'
#' Getting feature values requires looping through the samples and extracting
#' the chrom peak values for the respective features.
#'
#' @noRd
NULL


#' Convert a `XcmsExperiment` to an `XcmsExperimentHdf5` object: export all
#' data to a HDF5 file `h5_file` and return the `XcmsExperimentHdf5`.
#'
#' @noRd
.xcms_experiment_to_hdf5 <- function(x, h5_file = character()) {
    if (!length(h5_file))
        stop("Parameter 'h5_file' is mandatory.")
    .h5_initialize_file(h5_file)
    has_chrom_peaks <- hasChromPeaks(x)
    has_features <- hasFeatures(x)
    x <- as(x, "XcmsExperimentHdf5")
    x@sample_id <- .featureIDs(length(x), "S")
    x@hdf5_file <- h5_file
    mod_count <- 0L
    if (has_chrom_peaks) {
        message("Note: reformatting row names for the chromPeaks matrix.")
        ## Memory-efficient export: save the data for one sample at a time. If
        ## that is too slow we could split the data and export all in one go.
        is_sample <- colnames(x@chromPeaks) == "sample"
        msl <- unique(x@chromPeakData$ms_level)

        for (i in seq_along(x@sample_id)) {
            idx <- unname(which(x@chromPeaks[, is_sample] == i))
            pks <- x@chromPeaks[idx, !is_sample, drop = FALSE]
            pkd <- x@chromPeakData[idx, , drop = FALSE]
            f <- factor(pkd$ms_level, levels = msl)
            ## Update chrom peak IDs to the new format
            pks <- split.data.frame(pks, f)
            for (j in length(msl))
                rownames(pks[[j]]) <- .featureIDs(
                    nrow(pks[[j]]), paste0("CP", msl[j], x@sample_id[i]),
                    min_len = 6)
            pkd <- split.data.frame(
                pkd[, colnames(pkd) != "ms_level", drop = FALSE], f)
            names(pks) <- x@sample_id[i]
            names(pkd) <- x@sample_id[i]
            mod_count <- .h5_write_data(
                h5_file, pks, name = "chrom_peaks", ms_level = msl,
                replace = FALSE, write_colnames = TRUE, write_rownames = TRUE)
            mod_count <- .h5_write_data(
                h5_file, pkd, name = "chrom_peak_data", ms_level = msl,
                replace = FALSE, write_rownames = FALSE)
        }
        slot(x, "chromPeaks", check = FALSE) <-
            x@chromPeaks[integer(), , drop = FALSE]
        slot(x, "chromPeakData", check = FALSE) <-
            x@chromPeakData[integer(), , drop = FALSE]
        slot(x, "chrom_peaks_ms_level", check = FALSE) <- msl
    }
    if (has_features) {
        stop("Can not yet save feature definitions to HDF5")
        slot(x, "has_features", check = FALSE) <- TRUE
    }
    x@hdf5_mod_count <- mod_count
    x
}

#' Subset an `XcmsExperimentHdf5` object. Similar to `.subset_xcms_experiment()`
#' for `XcmsExperiment`, but optimized for `XcmsExperimentHdf5`.
#'
#' @noRd
.h5_subset_xcms_experiment <- function(x, i = integer(),
                                       keepChromPeaks = TRUE,
                                       keepAdjustedRtime = FALSE,
                                       keepFeatures = FALSE,
                                       ignoreHistory = FALSE,
                                       ...) {
    i <- i2index(i, length(x))
    if (any(i < 0)) {
        if (all(i < 0))
            i <- seq_along(x)[i]
        else stop("Mixing positive and negative indices is not supported.")
    }
    drop <- character()
    if (!keepAdjustedRtime && hasAdjustedRtime(x)) {
        svs <- unique(c(spectraVariables(x@spectra), "mz", "intensity"))
        x@spectra <- selectSpectraVariables(
            x@spectra, svs[svs != "rtime_adjusted"])
        drop <- c(drop, .PROCSTEP.RTIME.CORRECTION)
    }
    if (!keepFeatures && hasFeatures(x)) {
        stop("Subsetting with features present needs to be implemented")
        drop <- c(drop, .PROCSTEP.PEAK.GROUPING)
    }
    if (!keepChromPeaks && hasChromPeaks(x)) {
        x@chrom_peaks_ms_level <- integer()
        drop <- c(drop, .PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.FILLING,
                  .PROCSTEP.CALIBRATION, .PROCSTEP.PEAK.REFINEMENT)
    }
    if (!ignoreHistory && length(drop))
        x@processHistory <- dropProcessHistoriesList(
            x@processHistory, type = drop)
    x@sample_id <- x@sample_id[i]
    getMethod("[", "MsExperiment")(x, i = i)
}

################################################################################
##
##        CHROM PEAK RELATED THINGS
##
################################################################################

## findChromPeaks ->
## .mse_find_chrom_peaks_chunks ->
## .mse_spectrapply_chunks -> .mse_find_chrom_peaks_chunk (performs peak
## detection with Spectra as input)
## .mse_spectrapply_chunks: needs to get a function that also saves the
## results to hdf5.

#' This is equivalent to .mse_find_chrom_peaks_chunk, but instead of returning
#' the detected peaks saves them to the specified `hdf5_file`
#'
#' @param x `Spectra` of the samples from the present chunk
#'
#' @param h5_file HDF5 file name to which the results should be saved to.
#'
#' @param sample_id `character` with the names/IDs of (all!) samples (i.e.
#'     @sample_id)
#'
#' @param add `logical(1)` whether newly identified chromatographic peaks
#'     should be added to existing chromatographic peaks. If `add = FALSE` (the
#'     default) previous results get replaced. Note: the upstream function
#'     should ensure that there are chrom peaks for that ms level for `TRUE`.
#'
#' @noRd
.h5_find_chrom_peaks_chunk <- function(x, msLevel = 1L, param,
                                       h5_file = character(),
                                       sample_id = character(),
                                       add = FALSE,
                                       ...,
                                       BPPARAM = bpparam()) {
    chunk_sample_index <- unique(x$.SAMPLE_IDX)
    message("sample index: ", paste0(chunk_sample_index, collapse = ", "))
    res <- .mse_find_chrom_peaks_chunk(
        x, msLevel = msLevel, param = param, BPPARAM = BPPARAM)
    names(res) <- sample_id[chunk_sample_index]
    pkdl <- vector("list", length(res)) # chromPeakData list
    names(pkdl) <- names(res)
    for (i in seq_along(res)) {
        sid <- sample_id[chunk_sample_index[i]]
        max_index <- 0L
        rnames <- character()
        nr <- nrow(res[[i]])
        pkd <- data.frame(ms_level = rep(msLevel, nr),
                          is_filled = rep(FALSE, nr))
        if (add) {
            ## Need to load previous results and append to that.
            pks <- .h5_read_data(h5_file, id = sid, name = "chrom_peaks",
                                 ms_level = msLevel, read_colnames = TRUE,
                                 read_rownames = TRUE)[[1L]]
            rnames <- rownames(pks)
            max_index <- max(
                as.integer(sub(paste0("CP", msLevel, sid), "", rnames)))
            res[[i]] <- rbindFill(pks, res[[i]])
            pkd <- rbindFill(.h5_read_data(
                h5_file, id = sid, name = "chrom_peak_data",
                ms_level = msLevel, read_rownames = FALSE)[[1L]], pkd)
        }
        pkdl[[i]] <- pkd
        rownames(res[[i]]) <- c(
            rnames, .featureIDs(nr, paste0("CP", msLevel, sid),
                                from = max_index + 1L, min_len = 6))
    }
    .h5_write_data(h5_file, res, name = "chrom_peaks",
                   ms_level = rep(msLevel, length(res)), replace = TRUE,
                   write_colnames = TRUE, write_rownames = TRUE)
    .h5_write_data(h5_file, pkdl, name = "chrom_peak_data",
                   ms_level = rep(msLevel, length(res)), replace = TRUE,
                   write_rownames = FALSE)
}

#' Similar to `.xmse_merge_neighboring_peaks()` in XcmsExperiment-functions.R,
#' but this does not return the chromatograpic peaks but stores them into
#' the HDF5 file instead. The function needs also to update/define rownames
#' for the newly added chromatographic peaks.
#'
#' @param x `XcmsExperimentHdf5` object with potentially multiple samples.
#'
#' @noRd
.h5_xmse_merge_neighboring_peaks <- function(x, msLevel = 1L, expandRt = 2,
                                             expandMz = 0, ppm = 10,
                                             minProp = 0.75,
                                             BPPARAM = bpparam()) {
    keep <- msLevel(spectra(x)) == msLevel
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    ## Get the list of chromPeak data for x.
    pksl <- .h5_read_data(
        x@hdf5_file, id = x@sample_id, name = "chrom_peaks",
        ms_level = rep(msLevel, length(x@sample_id)),
        read_colnames = TRUE, read_rownames = TRUE)
    ## Get the max index of a chrom peak per sample
    max_index <- integer(length(pksl))
    for (i in seq_along(pksl))
        max_index[i] <- max(
            c(0L, as.integer(sub(paste0("CP", msLevel, x@sample_id[i]), "",
                                 rownames(pksl[[i]])))))
    ## Get the list of chromPeakData for x.
    pkdl <- .h5_read_data(
        x@hdf5_file, id = x@sample_id, name = "chrom_peak_data",
        ms_level = rep(msLevel, length(x@sample_id)), read_rownames = TRUE)
    ## Do refinement (in parallel)
    res <- bpmapply(
        .merge_neighboring_peaks2,
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel),
                        f = factor()), f), pksl, pkdl, split(rt, f),
        MoreArgs = list(expandRt = expandRt, expandMz = expandMz,
                        ppm = ppm, minProp = minProp),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    ## Replace data in hdf5 for samples with changed data.
    for (i in seq_along(res)) {
        l <- list(res[[i]]$chromPeaks)
        nas <- is.na(rownames(l[[1L]]))
        rownames(l[[1L]])[nas] <- .featureIDs(
            sum(nas), paste0("CP", msLevel, x@sample_id[i]), max_index[i] + 1L)
        names(l) <- x@sample_id[i]
        .h5_write_data(h5_file = x@hdf5_file, data_list = l,
                       name = "chrom_peaks", ms_level = msLevel,
                       replace = TRUE, write_colnames = FALSE,
                       write_rownames = TRUE)
        pkd <- res[[i]]$chromPeakData
        if (!any(colnames(pkd) == "merged"))
            pkd$merged <- FALSE
        pkd$merged[grep("^CP", rownames(pkd), invert = TRUE)] <- TRUE
        l <- list(pkd)
        names(l) <- x@sample_id[i]
        .h5_write_data(h5_file = x@hdf5_file, data_list = l,
                       name = "chrom_peak_data", ms_level = msLevel,
                       replace = TRUE, write_rownames = FALSE)
    }
}

#' Internal function to extract the `chromPeaks` `matrix` of `x`. Mandatory
#' variables are `x` and `msLevel`.
#'
#' @param x `XcmsExperimentHdf5` for which the `chromPeaks()` information
#'     should be returned. The function returns data for all samples in the
#'     object.
#'
#' @param msLevel `integer(1)` to restrict the extraction to selected MS
#'     level(s). MS level(s) **have** to be provided.
#'
#' @param columns optional `character` allowing to define a subset of columns
#'     from which the data should be returned.
#'
#' @param by_sample `logical(1)` whether a `list` of `chromPeak` matrices split
#'     per sample should be returned or the *conventional* matrix with an
#'     additional column `"sample"`.
#'
#' @return
#'
#' For `by_sample = TRUE`: a `list` of chrom peak matrices, one element for
#' each sample/MS level. This is useful only for internal functions that
#' process the data per sample to avoid unnecessary merging and splitting.
#'
#' For `by_sample = FALSE`: a `numeric` `matrix` with the chrom peak matrix. A
#' columns `"sample"` is added to indicate the sample from which the data is.
#'
#' @noRd
.h5_chrom_peaks <- function(x, msLevel = integer(), columns = character(),
                            rt = numeric(), mz = numeric(), ppm = 0,
                            type = "any", by_sample = TRUE) {
    if (length(columns)) {
        ## Get column names, convert column names to indices.
        cn <- .h5_chrom_peaks_colnames(x, msLevel = msLevel)
        idx_columns <- match(columns, cn)
        if (anyNA(idx_columns))
            stop("Column(s) ", paste0("\"", columns[is.na(idx_columns)], "\"",
                                      collapse = ", "), " not found",
                 call. = FALSE)
    } else idx_columns <- NULL
    ids <- rep(x@sample_id, length(msLevel))
    msl <- rep(msLevel, each = length(x@sample_id))
    res <- .h5_read_data(x@hdf5_file, id = ids, name = "chrom_peaks",
                         ms_level = msl, read_colnames = TRUE,
                         read_rownames = TRUE, j = idx_columns,
                         rt = rt, mz = mz, ppm = ppm, type = type)
    if (by_sample) {
        names(res) <- ids
        res
    } else {
        l <- vapply(res, nrow, 1L)
        cbind(do.call(rbind, res), sample = rep(match(ids, x@sample_id), l))
    }
}

#' Extract the `chromPeakData` data.frame. Using `peaks` allows to reduce memory
#' demand because only data from the specified chrom peaks is returned. This
#' assumes that `chromPeaks()` was called before to get the IDs of the peaks.
#'
#' @param x `XcmsExperimentHdf5`
#'
#' @param columns optional `character()` to define the columns to extract.
#'
#' @param peaks optional `character()` to define selected chromatographic peaks
#'     for which the data should be returned. If not specified data for all
#'     chrom peaks is returned.
#'
#' @param by_sample `logical(1)` whether results should be `rbind` or returned
#'     as a `list` of `data.frame`.
#'
#' @noRd
.h5_chrom_peak_data <- function(x, msLevel = integer(), columns = character(),
                                peaks = character(), by_sample = TRUE) {
    ids <- rep(x@sample_id, length(msLevel))
    msl <- rep(msLevel, each = length(x@sample_id))
    ## Eventually pass chrom peak ids along to read only specicic data...
    res <- .h5_read_data(x@hdf5_file, id = ids, name = "chrom_peak_data",
                         ms_level = msl, read_rownames = TRUE, peaks = peaks)
    if (by_sample) {
        names(res) <- ids
        res <- mapply(FUN = function(a, b) {
            a$ms_level <- b
            a
        }, res, msl, SIMPLIFY = FALSE)
        res
    } else {
        l <- vapply(res, nrow, 1L)
        cbind(do.call(rbind, res), ms_level = rep(msl, l))
    }
}

.h5_chrom_peaks_colnames <- function(x, msLevel = 1L) {
    rhdf5::h5read(x@hdf5_file,
                  name = paste0("/", x@sample_id[1L], "/ms_",
                                msLevel[1L], "/chrom_peaks_colnames"),
                  drop = TRUE)
}

.h5_chrom_peak_data_colnames <- function(x, msLevel = 1L) {
    h5 <- rhdf5::H5Fopen(x@hdf5_file)
    on.exit(rhdf5::H5Fclose(h5))
    c(.h5_dataset_names(
        paste0("/", x@sample_id[1L], "/ms_", msLevel, "/chrom_peak_data"), h5),
      "ms_level")
}

#' Replace the retention times of chrom peaks with new values, depending
#' on the provided rts. This function is used during retention time alignment
#'
#' @param id `character(1)` with the ID of the sample
#'
#' @param rt_old `numeric` with the original retention times
#'
#' @param rt_new `numeric` with the new retention times
#'
#' @param ms_level `integer` defining for which MS levels the retention times
#'     should be adjusted. Ideally for all!
#'
#' @param hdf5_file `character(1)` with the name of the HDF5 file.
#'
#' @return hdf5_count
#'
#' @noRd
.h5_update_rt_chrom_peaks_sample <- function(id, rt_old, rt_new, ms_level,
                                             hdf5_file) {
    ## loop over MS levels
    cnt <- 0L
    for (msl in ms_level) {
        ## read chrom peaks
        cp <- .h5_read_data(hdf5_file, id = id, name = "chrom_peaks",
                            ms_level = msl, read_colnames = TRUE,
                            read_rownames = FALSE)[[1L]]
        ## adjust chrom peak rt - use .applyRtAdjToChromPeaks for that.
        cp <- .applyRtAdjToChromPeaks(
            cbind(cp, sample = rep(1, nrow(cp))), rtraw = list(rt_old),
            rtadj = list(rt_new))
        l <- list(cp[, colnames(cp) != "sample", drop = FALSE])
        names(l) <- id
        ## replace chrom peaks
        cnt <- .h5_write_data(hdf5_file, data_list = l, "chrom_peaks",
                              ms_level = msl, replace = FALSE,
                              write_colnames = FALSE, write_rownames = FALSE)
    }
    cnt
}


################################################################################
##
##        FEATURES THINGS
##
################################################################################

## ## WE MIGHT ACTUALLY NOT NEED THIS!
## #' Get chrom peaks for features from one sample. Allows to define/subset by
## #' features (`i`) and select column(s) from the chrom peaks matrix to return.
## #'
## #' @param sample_id `character(1)` with the ID of the sample from which to
## #'     return the data.
## #'
## #' @param hdf5_file `character(1)` with the HDF5 file name
## #'
## #' @param ms_level `integer(1)` with the MS level of the features/chrom peaks
## #'
## #' @param i optional `integer` to select the features for which to return the
## #'     data.
## #'
## #' @param j optional `integer` defining the index of the column(s) to return.
## #'
## #' @return `matrix` with the chrom peak data, first column being the feature
## #'     index.
## #'
## #' @importFrom S4Vectors findMatches to
## #'
## #' @noRd
## .h5_feature_chrom_peaks_sample <- function(sample_id, hdf5_file, ms_level,
##                                            i = integer(), j = NULL) {
##     fidx <- .h5_read_data(hdf5_file, sample_id, name = "feature_to_chrom_peaks",
##                           ms_level = ms_level)[[1L]]
##     if (length(i)) {
##         hits <- findMatches(i, fidx[, 1L])
##         fidx <- fidx[to(hits), , drop = FALSE]
##     }
##     vals <- .h5_read_data(hdf5_file, sample_id, name = "chrom_peaks",
##                           ms_level = ms_level, i = fidx[, 2L], j = j)[[1L]]
##     cbind(fidx[, 1L], vals)
## }

#' Extracts feature values for one sample summing intensities for features
#' with multiple peaks assigned.
#'
#' @param hdf5_file `character(1)` with the HDF5 file name.
#'
#' @param sample_id `character(1)` with the sample ID.
#'
#' @param ms_level `integer(1)` with the MS level.
#'
#' @param n_features `integer(1)` with the total number of features for that
#'     MS level.
#'
#' @param method `character(1)` defining the method to be used to tackle
#'     features with multiple peaks.
#'
#' @param col_idx `integer` with the index of the peak columns that should be
#'     loaded and processed by the functions. The first index **must** be the
#'     index of the `value` column (i.e. the column that should be reported).
#'     For `method = "maxint"`, the second column should be the column defined
#'     with parameter `intensity`, i.e. the column with the intensity values
#'     to select the *larger* peak. For `method = "rtmed"` it should be the
#'     index of the column `"rt"`.
#'
#' @param filled `logical(1)` whether gap-filled values should be reported or
#'     removed.
#'
#' @param rtmed `numeric` with the `"rtmed"` column of the feature definitions.
#'     Only used (but required) for `method = "medret"`.
#'
#' @noRd
.h5_feature_values_sample <- function(sample_id, hdf5_file, ms_level,
                                      n_features, method,
                                      col_idx = integer(),
                                      filled = TRUE, rtmed, ...) {
    res <- rep(NA_real_, n_features)
    sid <- paste0("/", sample_id, "/ms_", ms_level)
    vals <- .h5_read_data(hdf5_file, sample_id, name = "chrom_peaks",
                          ms_level = ms_level, j = col_idx)[[1L]]
    fidx <- .h5_read_data(hdf5_file, sample_id, name = "feature_to_chrom_peaks",
                          ms_level = ms_level)[[1L]]
    ## remove gap-filled values
    if (!filled) {
        is_filled <- rhdf5::h5read(hdf5_file,
                                   paste0(sid, "/chrom_peak_data/is_filled"),
                                   drop = TRUE)
        vals[is_filled, 1L] <- NA_real_
    }
    ## set/assign single and multiple values.
    res[fidx[, 1L]] <- vals[fidx[, 2L], 1L]
    if (method == "medret") {
        ## calculate difference between feature and peak rt
        vals[fidx[, 2L], 2L] <- vals[fidx[, 2L], 2L] - rtmed[fidx[, 1L]]
    }
    ## handle duplicates
    f <- factor(fidx[, 1L], levels = seq_len(n_features))
    pk_idx <- split(fidx[, 2L], f)
    idx_multi <- which(lengths(pk_idx) > 1L)
    if (length(idx_multi)) {
        FUN <- switch(
            method,
            sum = function(z) sum(vals[z, 1L]),
            maxint = function(z) vals[z, 1L][which.max(vals[z, 2L])],
            medret = function(z) vals[z, 1L][which.min(abs(vals[z, 2L]))])
        res[idx_multi] <- vapply(pk_idx[idx_multi], FUN, 1.1)
    }
    res
}

#' Get feature values for a specific MS level.
#'
#' @noRd
.h5_feature_values_ms_level <- function(ms_level, x, method, value, intensity,
                                        filled = TRUE) {
    cn <- .h5_chrom_peaks_colnames(x, ms_level)
  col <- switch(method,
                  sum = value,
                  medret = c(value, "rt"),
                  maxint = c(value, intensity))
    if (!all(col %in% cn))
        stop("Not all requested columns available. Please make sure 'value' ",
             "and 'intensity' (if defined) are available columns in the ",
             "chrom peak matrix.", call. = FALSE)
    col_idx <- match(col, cn)
    rtmed <- rhdf5::h5read(x@hdf5_file,
                           paste0("/features/ms_", ms_level,
                                  "/feature_definitions/rtmed"), drop = TRUE)
    rn <- rhdf5::h5read(x@hdf5_file,
                        paste0("/features/ms_", ms_level,
                               "/feature_definitions_rownames"), drop = TRUE)
    res <- do.call(
        cbind, lapply(x@sample_id, .h5_feature_values_sample,
                      hdf5_file = x@hdf5_file, ms_level = ms_level,
                      n_features = length(rtmed), method = method,
                      col_idx = col_idx, filled = filled, rtmed = rtmed))
    rownames(res) <- rn
    res
}

################################################################################
##
##        ALIGNMENT RELATED FUNCTIONALITY
##
################################################################################


################################################################################
##
##        EIC/CHROMATOGRAMS FUNCTIONALITY
##
################################################################################
#' Read chromatograms for a set of samples (chunk) and adds chromatographic
#' peaks.
#'
#' @param x `XcmsExperimentHdf5` for one subset/chunk of data from which the
#'     data should be extracted
#'
#' @param index `integer` with the index of the current subset `x` in the *full*
#'     data set.
#'
#' @param ms_level `integer(1)` with the MS level.
#' @noRd
.h5_x_chromatogram <- function(x, index = seq_along(x), ms_level = 1L,
                               mz, rt, ppm = 0, chromPeaks = "any",
                               BPPARAM = bpparam()) {
    ## Get the chromatograms in parallel.
    chr <- as(chromatogram(as(x, "MsExperiment"), mz = mz,
                        rt = rt, BPPARAM = BPPARAM), "XChromatograms")
    js <- seq_len(nrow(chr))
    message("Processing chromatographic peaks")
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = ncol(chr) + 1L, clear = FALSE)
    for (i in seq_along(x)) {
        cp <- .h5_read_data(x@hdf5_file, x@sample_id[i], "chrom_peaks",
                            ms_level = ms_level, read_colnames = TRUE,
                            read_rownames = TRUE)[[1L]]
        for (j in js) {
            idx <- which(.is_chrom_peaks_within_mz_rt(
                cp, rt[j, ], mz[j, ], ppm, chromPeaks), useNames = FALSE)
            a <- cbind(cp[idx, , drop = FALSE], sample = rep(i, length(idx)))
            b <- .h5_read_data(
                x@hdf5_file, x@sample_id[i], "chrom_peak_data",
                ms_level = ms_level, read_colnames = TRUE, i = idx,
                read_rownames = FALSE)[[1L]]
            b$ms_level <- rep(ms_level, length(idx))
            rownames(b) <- rownames(a)
            tmp <- chr@.Data[j, i][[1L]]
            slot(tmp, "chromPeaks", check = FALSE) <- a
            slot(tmp, "chromPeakData", check = FALSE) <- as(b, "DataFrame")
            chr@.Data[j, i][[1L]] <- tmp
        }
        pb$tick()
    }
    pb$tick()
    if (hasFeatures(x, ms_level)) {
        stop("Not yet implemented")
        ## Somehow add features.
    }
    chr@.processHistory <- x@processHistory
    chr
}


################################################################################
##
##        HDF5 FUNCTIONALITY
##
################################################################################

.h5_have_rhdf5 <- function() {
    return(requireNamespace("rhdf5", quietly = TRUE))
}

.h5_require_rhdf5 <- function() {
    if (!.h5_have_rhdf5())
        stop("Package 'rhdf5' is required for this functionality. Please ",
             "install using 'BiocManager::install(\"rhdf5\")' and try again.")
}

##  --------  READING  --------

#' Reads a single `matrix` from the HDF5 file. Depending on
#' `read_colnames` and `read_rownames` also the rownames and colnames are read.
#' Not reading them has performance advantages.
#'
#' This function supports reading only subsets of the data from the HDF5 file,
#' which, surprisingly, has a negative impact on performance. Maybe chunking
#' might improve that behaviour. Thus, for now, the `.h5_read_matrix2()`
#' function should be used instead.
#'
#' @param name `character(1)` with the name of the data set to read.
#'
#' @param h5 HDF5 file handle
#'
#' @param index `list` with `integer` indices of the rows and columns to read
#'     only a subset of the data. Is passed directly to `rhdf5::h5read()`.
#'
#' @param read_colnames `logical(1)` whether column names should also be read
#'     and set.
#'
#' @param read_rownames `logical(1)` whether row names should be read and set.
#'
#' @return numeric `matrix`
#'
#' @noRd
.h5_read_matrix <- function(name, h5, index = list(NULL, NULL),
                            read_colnames = FALSE,
                            read_rownames = FALSE,
                            rownames = paste0(name, "_rownames")) {
    d <- rhdf5::h5read(h5, name = name, index = index)
    if (read_rownames)
        rownames(d) <- rhdf5::h5read(h5, name = rownames, drop = TRUE,
                                     index = index[1L])
    if (read_colnames)
        colnames(d) <- rhdf5::h5read(h5, name = paste0(name, "_colnames"),
                                     drop = TRUE, index = index[2L])
    d
}

#' Same functionality as `.h5_read_matrix`, but this one does the subsetting
#' in R, i.e. reads first the **full** matrix into R and does the subsetting
#' in R. For data matrices up to 100,000 rows and 3 columns this function
#' is faster.
#'
#' @noRd
.h5_read_matrix2 <- function(name, h5, index = list(NULL, NULL),
                             read_colnames = FALSE,
                             read_rownames = FALSE,
                             rownames = paste0(name, "_rownames")) {
    d <- rhdf5::h5read(h5, name = name)
    if (!is.null(index[[1L]]))
        d <- d[index[[1L]], , drop = FALSE]
    if (!is.null(index[[2L]]))
        d <- d[, index[[2L]], drop = FALSE]
    if (read_rownames)
        rownames(d) <- rhdf5::h5read(h5, name = rownames, drop = TRUE,
                                     index = index[1L])
    if (read_colnames)
        colnames(d) <- rhdf5::h5read(h5, name = paste0(name, "_colnames"),
                                     drop = TRUE, index = index[2L])
    d
}

.h5_read_chrom_peaks_matrix <- function(name, h5, index = list(NULL, NULL),
                                        read_colnames = FALSE,
                                        read_rownames = FALSE,
                                        rownames = paste0(name, "_rownames"),
                                        rt = numeric(), mz = numeric(),
                                        ppm = 0, type = "any") {
    read_colnames <- read_colnames || length(rt) > 0 || length(mz) > 0
    d <- .h5_read_matrix2(name, h5, index, read_colnames, read_rownames,
                          rownames)
    if (length(rt) | length(mz))
        d[.is_chrom_peaks_within_mz_rt(d, rt = rt, mz = mz,
                                       ppm = ppm, type = type), , drop = FALSE]
    else d
}

#' Read a single `data.frame` from the HDF5 file. With
#' `read_rownames = TRUE` also the row names are read and set, which requires
#' an additional reading step. Note that for a `data.frame` each column
#' is stored as a separate array/DATASET type for the `data.frame` GROUP. Thus,
#' a single column can be simply read by specifying the `name` accordingly,
#' e.g. <data set name>/<column name>.
#'
#' @param name `character(1)` the name of the data set.
#'
#' @param h5 HDF5 file handle
#'
#' @param index `list` with integer indices passed to rhdf5::h5read to read
#'     only a subset of the data. Since `rhdf5::h5read()` seems to not support
#'     parameter `index` for data sets the subsetting is done in R - but only
#'     for `index[[1L]]`, i.e. rows. `index[[2L]]` is currently IGNORED.
#'
#' @param read_rownames `logical(1)` whether rownames should be read and set.
#'
#' @param rownames `character(1)` defining the name of the HDF5 array
#'     containing the rownames.
#'
#' @noRd
.h5_read_data_frame <- function(name, h5, index = list(NULL, NULL),
                                read_rownames = FALSE,
                                rownames = paste0(name, "_rownames"), ...) {
    d <- rhdf5::h5read(h5, name = name)
    if (is.list(d))
        d <- lapply(d, as.vector)
    d <- as.data.frame(d)
    if (read_rownames)
        rownames(d) <- rhdf5::h5read(h5, rownames, drop = TRUE)
    if (is.null(index[[1L]]))
        d
    else d[index[[1L]], , drop = FALSE]
}

.h5_read_chrom_peak_data <- function(name, h5, index = list(NULL, NULL),
                                     read_rownames = FALSE, peaks = character(),
                                     ...) {
    cd <- .h5_read_data_frame(
        name, h5, read_rownames = read_rownames || length(peaks) > 0,
        index = index, rownames = sub("_data", "s_rownames", name))
    if (length(peaks)) {
        idx <- match(peaks, rownames(cd))
        cd <- cd[idx[!is.na(idx)], , drop = FALSE]
        if (!read_rownames)
            rownames(cd) <- NULL
    }
    cd
}

#' Reads the names of the data sets of a group. This can for example be used
#' to get the column names of a `data.frame` that was saved as a GROUP of
#' DATASETs.
#'
#' @noRd
.h5_dataset_names <- function(name, h5, recursive = FALSE) {
    g <- rhdf5::H5Gopen(h5, name)
    on.exit(rhdf5::H5Gclose(g))
    rhdf5::h5ls(g, recursive = recursive, datasetinfo = FALSE)$name
}


.h5_ms_levels <- function(h5, sample_id) {
    nms <- .h5_dataset_names(paste0("/", sample_id), h5)
    as.integer(unique(sub("ms_", "", grep("^ms", nms, value = TRUE))))
}

.h5_chrom_peak_ms_levels <- function(h5_file, sample_id) {
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(rhdf5::H5Fclose(h5))
    msl <- .h5_ms_levels(h5, sample_id)
    has_cp <- vapply(
        paste0("/", sample_id, "/ms_", msl, "/"),
        function(x) {
            any(.h5_dataset_names(x, h5) == "chrom_peaks")
        }, NA)
    msl[has_cp]
}

#' Read selected datasets from a HDF5 file.
#'
#' @note
#'
#' Setting row and column names requires additional read steps and is hence
#' considerably slower than *just* importing the data. Also, if possible,
#' consider importing only the required column(s) using parameter `column`.
#'
#' @param h5_file `character(1)` with the HDF5 file name
#'
#' @param id `character` with the ID(s) of the data sets to read.
#'
#' @param name `character(1)` specifying which data should be read.
#'
#' @param ms_level `integer` with the MS level of each sample/data set that
#'     should be read. Has to have the same length than `id`.
#'
#' @param read_colnames `logical(1)` whether column names should be read and
#'     set for each `matrix`.
#'
#' @param read_rownames `logical(1)` whether row names should be read and
#'     set for each `matrix` or `data.frame`
#'
#' @param rownames `logical(1)` defining the name of the HDF5 array containing
#'     the row names.
#'
#' @param i optional index with the rows to read.
#'
#' @param j For `name = "chrom_peak_data"`: `character(1)` allowing to
#'     select a **single** column to read. For `name = "chrom_peaks"`: `integer`
#'     with the indices of the column(s) that should be imported.
#'
#' @param ... additional parameters passed to `FUN`
#'
#' @return `list()` with the read datasets. Will be a `list` of `numeric`
#'     matrices for `name = "chrom_peaks"` or a `list` with `data.frame`s for
#'     `name = "chrom_peak_data"`.
#'
#' @noRd
.h5_read_data <- function(h5_file = character(),
                          id = character(),
                          name = c("chrom_peaks", "chrom_peak_data",
                                   "feature_definitions",
                                   "feature_to_chrom_peaks"),
                          ms_level = integer(),
                          read_colnames = FALSE,
                          read_rownames = FALSE,
                          i = NULL, j = NULL, ...) {
    if (!length(id)) return(list())
    stopifnot(length(ms_level) == length(id))
    name <- match.arg(name)
    FUN <- switch(name,
                  chrom_peak_data = .h5_read_chrom_peak_data,
                  feature_definitions = .h5_read_data_frame,
                  chrom_peaks = .h5_read_chrom_peaks_matrix,
                  .h5_read_matrix2)
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    d <- paste0("/", id, "/ms_", ms_level, "/", name)
    index <- list(i, j)
    if (is.character(j) && length(j) == 1L) {
        d <- paste0(d, "/", j)
        index <- list(i, NULL)
    }
    lapply(d, FUN = FUN, read_colnames = read_colnames,
           read_rownames = read_rownames, index = index, h5 = h5, ...)
}

##  --------  VALIDITY  --------

#' Compares the "mod_count" attribute from an h5 file with the expected
#' one. Throws an error if it is different.
#'
#' The `mod_count` gets incremented by any operation that writes data to the
#' HDF5 file. If the `mod_count` of the xcms result object is different than
#' the one in the HDF5 file an error is thrown.
#'
#' @param h5 file handle for the HDF5 file
#'
#' @param mod_count `integer(1)` the expected value for the modification
#'     counter.
#'
#' @return `TRUE` or throws an error.
#'
#' @noRd
.h5_check_mod_count <- function(h5, mod_count = 0L) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L]
    if (mc != mod_count)
        stop("The HDF5 file was changed by a different process. This xcms ",
             "result object/variable is no longer valid.")
    TRUE
}

#' Checks the HDF5 file for validity.
#'
#' @param h5_file `character(1)` with the file name
#'
#' @param mod_count `integer(1)` with the expected value for the modification
#'     counter
#'
#' @return `TRUE` if the file is valid or throws an error if the file is
#'     not valid.
#'
#' @noRd
.h5_valid_file <- function(h5_file = character(), mod_count = 0L) {
    if (!length(h5_file))
        stop("'hdf5_file' missing with no default")
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    if (!rhdf5::H5Lexists(h5, "/header/package"))
        stop("File \"", h5_file, "\" is not in correct format.")
    res <- rhdf5::h5read(h5, "/header/package")
    if (res != "package:xcms")
        stop("File \"", h5_file, "\" is not in correct format.")
    .h5_check_mod_count(h5, mod_count = mod_count)
    TRUE
}

##  --------  WRITING  --------

.h5_compression_level <- function() 0L

.h5_filter <- function() "NONE"

#' Initializes the HDF5 file
#'
#' @noRd
.h5_initialize_file <- function(x, mod_count = 0L) {
    if (file.exists(x))
        stop("File \"", x, "\" already exists. Please choose a different name ",
             "or remove that file first.")
    h5 <- rhdf5::H5Fcreate(x)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    comp_level <- .h5_compression_level()
    flt <- .h5_filter()
    rhdf5::h5createGroup(h5, "header")
    rhdf5::h5write("package:xcms", h5, "/header/package",
                   level = comp_level)
    rhdf5::h5write(mod_count, h5, "/header/modcount",
                   level = comp_level)
}

#' Every writing operation should increate the "mod count", i.e. the
#' count of data modifications. This function increases the mod count by +1
#' and returns this value.
#'
#' @noRd
.h5_increment_mod_count <- function(h5) {
    mc <- rhdf5::h5read(h5, "/header/modcount")[1L] + 1L
    rhdf5::h5write(mc, h5, "/header/modcount",
                   level = .h5_compression_level())
    mc
}

#' Bare writing function of a `matrix`
#'
#' @param x `matrix` with e.g. the chrom peak results **of a single sample**
#'
#' @param h5 HDF5 file handle.
#'
#' @param name `character(1)` with the name for the data set (e.g.
#'     `"/S1/ms_1/chrom_peaks"`.
#'
#' @param level `integer(1)` with the compression level.
#'
#' @param write_colnames `logical(1)` whether to write the column names of `x`
#'     as an additional data set `paste0(name, "_colnames")`.
#'
#' @param write_rownames `logical(1)` whether to write the rownames of `x` as
#'     an additional data set `paste0(name, "_rownames")`.
#'
#' @param replace `logical(1)` whether an eventually existing data set should
#'     be replaced or updated. Data sets can only be updated if their
#'     dimensions are identical.
#'
#' @noRd
.h5_write_matrix <- function(x, h5, name, level, write_colnames = TRUE,
                             write_rownames = TRUE, replace = TRUE) {
    if (replace && rhdf5::H5Lexists(h5, name))
        rhdf5::h5delete(h5, name)
    rhdf5::h5write(x, h5, name = name, level = level,
                   write.attributes = FALSE,
                   createnewfile = FALSE)
    if (write_rownames) {
        dn <- paste0(name, "_rownames")
        if (replace && rhdf5::H5Lexists(h5, dn))
            rhdf5::h5delete(h5, dn)
        rhdf5::h5write(rownames(x), h5, name = dn,
                       level = level, createnewfile = FALSE)
    }
    if (write_colnames) {
        dn <- paste0(name, "_colnames")
        if (replace && rhdf5::H5Lexists(h5, dn))
            rhdf5::h5delete(h5, dn)
        rhdf5::h5write(colnames(x), h5, name = dn,
                       level = level, createnewfile = FALSE)
    }
}

#' Bare writing function of a data.frame
#'
#' @param x `data.frame`
#'
#' @param h5 HDF5 file handle.
#'
#' @param name `character(1)` defining the name for the data set.
#'
#' @param level `integer(1)` with the compression level.
#'
#' @param replace `logical(1)` whether an eventually existing data set should
#'     be replaced or updated. Data sets can only be updated if their
#'     dimensions are identical.
#'
#' @noRd
.h5_write_data_frame <- function(x, h5, name, level, replace = TRUE,
                                 write_rownames = FALSE, ...) {
    if (replace && rhdf5::H5Lexists(h5, name))
        rhdf5::h5delete(h5, name)
    rhdf5::h5writeDataset(x, h5, name, level = level,
                          DataFrameAsCompound = FALSE)
    if (write_rownames) {
        dn <- paste0(name, "_rownames")
        if (replace && rhdf5::H5Lexists(h5, dn))
            rhdf5::h5delete(h5, dn)
        rhdf5::h5write(rownames(x), h5, name = dn,
                       level = level, createnewfile = FALSE)
    }
}

#' Writes data (matrix, data.frame) to a HDF5 file organized by sample and
#' MS level:
#'
#' /<sample id>/ms_<MS level>/<data set>
#'
#' Example:
#'
#' /S001/ms_1/chrom_peaks
#' /S001/ms_2/chrom_peaks
#' /S002/ms_1/chrom_peaks
#' /S002/ms_2/chrom_peaks
#' ...
#'
#' @param h5_file the HDF5 file.
#'
#' @param data_list `list` of data sets (`matrix` or `data.frame`) that should
#'     be written. The names of the `list` are used as sample identifier.
#'
#' @param name `character(1)` with the name of the data set. This might be
#'     changed in future.
#'
#' @param ms_level `integer` with **same length** than `data_list` defining
#'     the MS level of the data sets.
#'
#' @param replace `logical(1)` whether eventual existing data sets should be
#'     replaced or updated. Data must be replaced if the dimensions are
#'     different.
#'
#' @noRd
.h5_write_data <- function(h5_file = character(),
                           data_list = list(),
                           name = c("chrom_peaks", "chrom_peak_data",
                                    "feature_definitions",
                                    "feature_to_chrom_peaks"),
                           ms_level = integer(),
                           replace = TRUE,
                           write_colnames = TRUE,
                           write_rownames = TRUE) {
    if (!length(data_list)) return(TRUE)
    stopifnot(length(ms_level) == length(data_list))
    stopifnot(!is.null(names(data_list)))
    name <- match.arg(name)
    h5 <- rhdf5::H5Fopen(h5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    FUN <- .h5_write_matrix
    if (name %in% c("chrom_peak_data", "feature_definitions"))
        FUN <- .h5_write_data_frame
    comp_level <- .h5_compression_level()
    flt <- .h5_filter()
    for (i in seq_along(data_list)) {
        group_sample <- paste0("/", names(data_list)[i])
        if (!rhdf5::H5Lexists(h5, group_sample))
            rhdf5::h5createGroup(h5, group_sample)
        group_ms <- paste0(group_sample, "/ms_", ms_level[i])
        if (!rhdf5::H5Lexists(h5, group_ms))
            rhdf5::h5createGroup(h5, group_ms)
        group_data <- paste0(group_ms, "/", name)
        FUN(data_list[[i]], h5, group_data, comp_level,
            write_colnames = write_colnames, write_rownames = write_rownames,
            replace = replace)
    }
    .h5_increment_mod_count(h5)
}