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

#' @rdname XcmsExperimentHdf5
toXcmsExperimentHdf5 <- function(object, hdf5File = tempfile()) {
    if (missing(object) || !is(object, "XcmsExperiment"))
        stop("'object' is expected to be an instance of class 'XcmsExperiment'")
    .xcms_experiment_to_hdf5(object, h5_file = hdf5File)
}

#' @rdname XcmsExperimentHdf5
toXcmsExperiment <- function(object, ...) {
    if (missing(object) || !is(object, "XcmsExperimentHdf5"))
        stop("'object' is expected to be an instance of ",
             "class 'XcmsExperimentHdf5'")
    .h5_to_xcms_experiment(object)
}

#' Coerce from XcmsExperimentHdf5 to XcmsExperiment
#'
#' @noRd
.h5_to_xcms_experiment <- function(from) {
    res <- as(as(from, "MsExperiment"), "XcmsExperiment")
    if (hasChromPeaks(from)) {
        res@chromPeaks <- chromPeaks(from)
        res@chromPeakData <- chromPeakData(from, return.type = "data.frame")
    }
    if (hasFeatures(from)) {
        map <- do.call(rbind, lapply(from@features_ms_level, function(m) {
            tmp <- .h5_read_data(from@hdf5_file, from@sample_id,
                                 "feature_to_chrom_peak",
                                 rep(m, length(from@sample_id)))
            fid <- .h5_feature_definitions_rownames(from, m)[[1L]]
            cid <- .h5_chrom_peaks_rownames(from, m)
            do.call(rbind, mapply(tmp, cid, FUN = function(a, b) {
                cbind(fid[a[, 1L]], b[a[, 2L]])
            }))
        }))
        fd <- featureDefinitions(from)
        pkidx <- match(map[, 2L], rownames(res@chromPeaks))
        fd$peakidx <- split(pkidx, factor(map[, 1L], levels = rownames(fd)))
        res@featureDefinitions <- fd
    }
    res@processHistory <- from@processHistory
    res
}

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
    gap_peaks_ms_level <- integer()
    x <- as(x, "XcmsExperimentHdf5")
    x@sample_id <- .featureIDs(length(x), "S")
    x@hdf5_file <- h5_file
    mod_count <- 0L
    if (has_features)
        pidx <- split(x@featureDefinitions$peakidx,
                      x@featureDefinitions$ms_level)
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
            for (j in seq_along(msl)) {
                rownames(pks[[j]]) <- .featureIDs(
                    nrow(pks[[j]]), paste0("CP", msl[j], x@sample_id[i]),
                    min_len = 6)
            }
            pkd <- split.data.frame(
                pkd[, colnames(pkd) != "ms_level", drop = FALSE], f)
            names(pks) <- rep(x@sample_id[i], length(pks))
            names(pkd) <- rep(x@sample_id[i], length(pkd))
            gaps <- vapply(pkd, function(z) any(z$is_filled), NA)
            gap_peaks_ms_level <- union(gap_peaks_ms_level, msl[gaps])
            mod_count <- .h5_write_data(
                h5_file, pks, name = "chrom_peaks", ms_level = msl,
                replace = FALSE, write_colnames = TRUE, write_rownames = TRUE)
            mod_count <- .h5_write_data(
                h5_file, pkd, name = "chrom_peak_data", ms_level = msl,
                replace = FALSE, write_rownames = FALSE)
            if (has_features) {
                ## Create chrom peak to feature mapping per MS level
                idx <- split(idx, f)
                for (j in names(pidx)) {
                    map <- lapply(pidx[[j]], function(z) {
                        mtch <- match(z, idx[[j]])
                        mtch[!is.na(mtch)]
                    })
                    l <- list(cbind(rep(seq_along(map), lengths(map)),
                                    unlist(map, use.names = FALSE)))
                    names(l) <- x@sample_id[i]
                    .h5_write_data(
                        h5_file, l, name = "feature_to_chrom_peaks",
                        write_colnames = FALSE, write_rownames = FALSE,
                        ms_level = j)
                }
            }
        }
        slot(x, "chromPeaks", check = FALSE) <-
            x@chromPeaks[integer(), , drop = FALSE]
        slot(x, "chromPeakData", check = FALSE) <-
            x@chromPeakData[integer(), , drop = FALSE]
        slot(x, "chrom_peaks_ms_level", check = FALSE) <- msl
        slot(x, "gap_peaks_ms_level", check = FALSE) <- gap_peaks_ms_level
    }
    if (has_features) {
        message("Note: reformatting feature IDs.")
        ## Export the feature definitions, split by MS level
        msl <- unique(x@featureDefinitions$ms_level)
        for (m in msl) {
            fd <- extractROWS(x@featureDefinitions,
                              which(x@featureDefinitions$ms_level == m))
            fd$peakidx <- NULL
            fd$ms_level <- NULL
            attr(fd, "row.names") <- .featureIDs(
                nrow(fd), prefix = paste0("FT", m), min_len = 6)
            mod_count <- .h5_write_data(
                h5_file, list(features = fd), "feature_definitions",
                m, replace = TRUE, write_rownames = TRUE)
        }
        x@features_ms_level <- msl
        slot(x, "featureDefinitions", check = FALSE) <-
            x@featureDefinitions[integer(), , drop = FALSE]
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
    if (!keepChromPeaks && hasChromPeaks(x)) {
        x@chrom_peaks_ms_level <- integer()
        x@gap_peaks_ms_level <- integer()
        drop <- c(drop, .PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.FILLING,
                  .PROCSTEP.CALIBRATION, .PROCSTEP.PEAK.REFINEMENT)
        keepFeatures <- FALSE
    }
    if (!keepFeatures && hasFeatures(x)) {
        x@features_ms_level <- integer()
        drop <- c(drop, .PROCSTEP.PEAK.GROUPING)
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
                        f = factor(), return.type = "list"), f),
        pksl, pkdl, split(rt, f),
        MoreArgs = list(expandRt = expandRt, expandMz = expandMz,
                        ppm = ppm, minProp = minProp),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    ## Replace data in hdf5 for samples with changed data.
    for (i in seq_along(res)) {
        l <- list(res[[i]]$chromPeaks)
        nas <- is.na(rownames(l[[1L]]))
        rownames(l[[1L]])[nas] <- .featureIDs(
            sum(nas), paste0("CP", msLevel, x@sample_id[i]),
            max_index[i] + 1L, min_len = 6)
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

#' Perform peak integration on provided m/z - RT regions. Can be called for
#' gap filling or manual definition of chrom peaks.
#'
#' The newly defined chrom peaks will be appended to eventually existing chrom
#' peak array in the HDF5 file. Also, content to the chrom peak data of the
#' respective peaks will be added to the HDF5 file.
#'
#' For gap-filling, the rownames of the individual matrices in `pal` need to
#' represent the feature IDs and `update_features` needs to be set to `TRUE`.
#'
#' @param x `XcmsExperimentHdf5` with data from potentially multiple files.
#'
#' @param pal `list` of peak area matrices defining (for each individual sample
#'     in `x`) the MS area from which the signal should be integrated.
#'     `lengh(pal)` must be equal to `length(x)`. Names should represent the
#'     file/sample **index**!
#'
#' @param msLevel `integer(1)` with the MS level on which to integrate the data.
#'
#' @param intFun function to be used for the integration.
#'
#' @param mzCenterFun function to calculate the m/z value
#'
#' @param update_features `logical(1)` whether feature to chrom peak mappings
#'     should be updated too. If `TRUE` (i.e. for gap filling with
#'     `fillChromPeaks()`) the rownames of the individual peak area definitions
#'     in `pal` need to correspond to feature IDs.
#'
#' @param is_filled `logical(1)` with the value to be used in the
#'     `chromPeakData`'s `"is_filled"` column. Should be `TRUE` for gap-filling
#'     and `FALSE` for *manual chrom peaks*.
#'
#' @return `integer(1)` being the hdf5 counter that keeps track of the number
#'     of writing operations to the HDF5 file.
#'
#' @noRd
.h5_xmse_integrate_chrom_peaks <-
    function(x, pal, msLevel = 1L, intFun = .chrom_peak_intensity_centWave,
             mzCenterFun = "mzCenter.wMean", param = MatchedFilterParam(),
             BPPARAM = bpparam(), update_features = FALSE, is_filled = TRUE,
             ...) {
        keep <- which(msLevel(spectra(x)) == msLevel)
        f <- as.factor(fromFile(x)[keep])
        if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
        else rt <- rtime(spectra(x))[keep]
        cn <- c(.h5_chrom_peaks_colnames(x, msLevel = msLevel), "sample")
        res <- bpmapply(
            split(peaksData(filterMsLevel(spectra(x), msLevel), f = factor(),
                            return.type = "list"), f),
            split(rt, f),
            pal,
            as.integer(names(pal)),
            FUN = intFun,
            MoreArgs = list(mzCenterFun = mzCenterFun, cn = cn, param = param),
            SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
        if (update_features) {
            fids <- .h5_feature_definitions_rownames(x, msLevel)[[1L]]
            feature_idx <- lapply(res, function(z) match(rownames(z), fids))
        }
        ## Update HDF5 file content per sample
        mc <- x@hdf5_mod_count
        for (i in seq_along(res)) {
            ## chrom peaks
            pks <- .h5_read_data(
                x@hdf5_file, id = x@sample_id[i], ms_level = msLevel,
                name = "chrom_peaks", read_colnames=TRUE,
                read_rownames = TRUE)[[1L]]
            prefix <- paste0("CP", msLevel, x@sample_id[i])
            max_index <- max(as.integer(sub(prefix, "", rownames(pks))))
            rownames(res[[i]]) <- .featureIDs(
                nrow(res[[i]]), prefix, max_index + 1L, min_len = 6)
            l <- list(rbindFill(pks, res[[i]][, colnames(res[[i]]) !="sample",
                                              drop = FALSE]))
            names(l) <- x@sample_id[i]
            rm(pks)
            mc <- .h5_write_data(
                x@hdf5_file, l, name = "chrom_peaks", ms_level = msLevel,
                replace = TRUE, write_colnames = TRUE, write_rownames = TRUE)
            ## chrom peak data
            pkd <- .h5_read_data(
                x@hdf5_file, id = x@sample_id[i], ms_level = msLevel,
                name = "chrom_peak_data", read_colnames = TRUE)[[1L]]
            l <- list(rbindFill(
                pkd, data.frame(is_filled = rep(is_filled, nrow(res[[i]])),
                                merged = rep(FALSE, nrow(res[[i]])))))
            names(l) <- x@sample_id[i]
            mc <- .h5_write_data(
                x@hdf5_file, l, name = "chrom_peak_data", ms_level = msLevel,
                replace = TRUE, write_rownames = FALSE)
            ## feature to chrom peak mapping
            if (update_features) {
                fmap <- rbind(
                    .h5_read_data(
                        x@hdf5_file, x@sample_id[i],
                        "feature_to_chrom_peaks", msLevel)[[1L]],
                    cbind(feature_idx[[i]],
                          seq(nrow(pkd) + 1, length.out = nrow(res[[i]])))
                )
                l <- list(fmap[order(fmap[, 1L]), , drop = FALSE])
                names(l) <- x@sample_id[i]
                mc <- .h5_write_data(
                    x@hdf5_file, l, name = "feature_to_chrom_peaks",
                    write_colnames = FALSE, write_rownames = FALSE,
                    ms_level = msLevel)
            }
        }
        mc
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
    if (by_sample) {
        sample_idx <- integer()
    } else {
        ## pass the sample index to the import function to add the "sample" col
        sample_idx <- match(ids, x@sample_id)
        names(sample_idx) <- paste0("/", ids, "/ms_", msl, "/chrom_peaks")
    }
    res <- .h5_read_data(x@hdf5_file, id = ids, name = "chrom_peaks",
                         ms_level = msl, read_colnames = TRUE,
                         read_rownames = TRUE, j = idx_columns,
                         rt = rt, mz = mz, ppm = ppm, type = type,
                         sample_index = sample_idx)
    if (by_sample) {
        names(res) <- ids
        res
    } else {
        do.call(base::rbind, res)
    }
}

#' Extract the `chromPeakData` data.frame. Using `peaks` allows to reduce memory
#' demand because only data from the specified chrom peaks is returned. This
#' assumes that `chromPeaks()` was called before to get the IDs of the peaks.
#' We're using the data.table::rbindlist to combine the data.frames because its
#' much faster - but unfortunately drops also the rownames.
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
#' @importFrom data.table rbindlist
#'
#' @noRd
.h5_chrom_peak_data <- function(x, msLevel = integer(), columns = character(),
                                peaks = character(), by_sample = TRUE) {
    ids <- rep(x@sample_id, length(msLevel))
    msl <- rep(msLevel, each = length(x@sample_id))
    names(msl) <- paste0("/", ids, "/ms_", msl, "/chrom_peak_data")
    res <- .h5_read_data(x@hdf5_file, id = ids, name = "chrom_peak_data",
                         ms_level = msl, read_rownames = TRUE, peaks = peaks,
                         ms_levels = msl)
    if (by_sample) {
        names(res) <- ids
        res
    } else {
        rn <- unlist(lapply(res, rownames), use.names= FALSE, recursive = FALSE)
        res <- base::as.data.frame(data.table::rbindlist(res))
        attr(res, "row.names") <- rn
        res
    }
}

.h5_chrom_peaks_colnames <- function(x, msLevel = 1L) {
    rhdf5::h5read(x@hdf5_file,
                  name = paste0("/", x@sample_id[1L], "/ms_",
                                msLevel[1L], "/chrom_peaks_colnames"),
                  drop = TRUE)
}

.h5_chrom_peaks_rownames <- function(x, msLevel = x@chrom_peaks_ms_level) {
    ids <- rep(x@sample_id, length(msLevel))
    msl <- rep(msLevel, each = length(x@sample_id))
    h5 <- rhdf5::H5Fopen(x@hdf5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    d <- paste0("/", ids, "/ms_", msl, "/chrom_peaks_rownames")
    lapply(d, FUN = rhdf5::h5read, file = h5, drop = TRUE)
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

#' Filter the chrom peak HDF5 entries based on a user provided filter function.
#' If present, this function updates also the chrom peak to feature mapping and
#' feature definitions entries.
#'
#' Functions that can be used for `FUN`:
#' - `.which_in_range()`: filter based on rt or m/z range.
#'
#' @param x `XcmsExperimentHdf5`.
#'
#' @param FUN function to filter each chrom peak matrix. Additional parameters
#'     to this functions are provided through `...`. This function is expected
#'     to take either the chrom peaks `matrix` or chrom peak data `data.frame`
#'     and returns an `integer` with the index of the rows to keep.
#'
#' @param chrom_peak_data `logical(1)` whether FUN should be applied to the
#'     chrom peak `matrix` (the default) or the chrom peak data `data.frame`
#'     (`chrom_peak_data = TRUE`).
#'
#' @return the function returns the *mod counter* representing eventual updates
#'     to the file content.
#'
#' @noRd
.h5_filter_chrom_peaks <- function(x, msLevel, FUN = NULL,
                                   chrom_peak_data = FALSE, ...) {
    mc <- x@hdf5_mod_count
    keep_features <- vector("list", length(x))
    names(keep_features) <- x@sample_id
    h5f <- x@hdf5_file
    data_changed <- FALSE # keep track if there was any change in data
    for (msl in msLevel) {
        handle_features <- any(x@features_ms_level %in% msl)
        for (id in x@sample_id) {
            pks <- .h5_read_data(h5f, id, "chrom_peaks", msl,
                                 read_colnames =TRUE, read_rownames =TRUE)[[1L]]
            pkd <- .h5_read_data(h5f, id, "chrom_peak_data", msl,
                                 read_colnames = TRUE)[[1L]]
            if (chrom_peak_data) idx <- FUN(pkd, ...)
            else idx <- FUN(pks, ...)
            if (handle_features) {
                fmap <- .h5_read_data(
                    h5f, id, "feature_to_chrom_peaks", msl)[[1L]]
                fmap <- fmap[fmap[, 2L] %in% idx, , drop = FALSE]
                keep_features[[id]] <- cbind(fmap[, 1L], match(fmap[, 2L], idx))
            }
            if (.is_equal(idx, 1:nrow(pks))) next # skip export
            ## Subset and export
            data_changed <- TRUE
            l <- list(pks[idx, , drop = FALSE])
            names(l) <- id
            .h5_write_data(h5f, l, "chrom_peaks", msl, replace = TRUE,
                           write_colnames = TRUE, write_rownames = TRUE)
            l <- list(extractROWS(pkd, idx))
            names(l) <- id
            mc <- .h5_write_data(h5f, l, "chrom_peak_data", msl,
                                 replace = TRUE, write_rownames = FALSE)
        }
        if (handle_features && data_changed) {
            fd <- .h5_read_data(h5f, "features", "feature_definitions",
                                read_rownames = TRUE, ms_level = msl)[[1L]]
            keep <- sort(unique(unlist(
                lapply(keep_features, `[`, , j = 1L), use.names = FALSE)))
            fd <- extractROWS(fd, keep)
            .h5_write_data(h5f, list(features = fd), "feature_definitions",
                           msl, replace = TRUE, write_rownames = TRUE)
            keep_features <- lapply(keep_features,
                                    function(z, i) {
                                        z[, 1L] <- match(z[, 1L], keep)
                                        z
                                    }, i = keep)
            mc <- .h5_write_data(
                h5f, keep_features, "feature_to_chrom_peaks",
                rep(msl, length(x)), write_colnames = FALSE,
                write_rownames = FALSE)
        }
    }
    mc
}

#' filter feature definitions. The function replaces the feature definitions
#' entry in the HDF5 file and iterates over all samples to update the feature
#' to chrom peak assignments.
#'
#' @param x `XcmsExperimentHdf5`
#'
#' @param msLevel `integer` with the MS level(s) in which to subset the
#'     features
#'
#' @param feature_id `character` with the IDs (rownames) of the features that
#'     should be retained. All other features are removed.
#'
#' @noRd
.h5_filter_feature_definitions <- function(x, msLevel, feature_id) {
    mc <- x@hdf5_mod_count
    h5f <- x@hdf5_file
    for (msl in msLevel) {
        fd <- .h5_read_data(h5f, "features", "feature_definitions",
                            read_rownames = TRUE, ms_level = msl)[[1L]]
        keep <- which(rownames(fd) %in% feature_id)
        fd <- extractROWS(fd, keep)
        mc <- .h5_write_data(h5f, list(features = fd), "feature_definitions",
                             msl, replace = TRUE, write_rownames = TRUE)
        for (id in x@sample_id) {
            fmap <- .h5_read_data(
                h5f, id, "feature_to_chrom_peaks", msl)[[1L]]
            fmap <- fmap[fmap[, 1L] %in% keep, , drop = FALSE]
            fmap[, 1L] <- match(fmap[, 1L], keep)
            fmap <- list(fmap)
            names(fmap) <- id
            mc <- .h5_write_data(h5f, fmap, "feature_to_chrom_peaks",
                                 msl, write_colnames = FALSE,
                                 write_rownames = FALSE)
        }
    }
    mc
}

#' Get spectra for chromatographic peaks of a specified sample. This function
#' is used by `chromPeakSpectra()` and `featureSpectra()`
#'
#' @param h5_file `character(1)` with the HDF5 file name
#'
#' @param id `character(1)` with the ID of the sample from which to process
#'     the data.
#'
#' @param s `Spectra` of the **current** sample (`id`) only!
#'
#' @param method `character(1)` defining which spectra to extract.
#'
#' @param msLevel `integer(1)` with the MS level of the spectra to extract
#'
#' @param chromPeaksMsLevel `integer(1)` with the MS level of the chrom peaks
#'     for which `Spectra` should be retrieved.
#'
#' @param expandRd, expandMz, ppm: expand ranges.
#'
#' @param skipFilled `logical(1)` whether to skip gap-filled peaks.
#'
#' @param peaks `integer()` with the index of the chrom peaks to extract.
#'     Will process all chrom peaks if not provided. If `character()` is
#'     provided it will be matched against the row names of the chrom peak
#'     matrix. If provided as `integer` indices, these are expected to be
#'     between 1 and `nrow()` of the chrom peaks for that particular sample. If
#'     provided as `character` it can also be IDs of chrom peaks from another
#'     sample.
#'
#' @param chromPeakColumns `character` with optional columns to add to the
#'     `Spectra`.
#'
#' @return always a `Spectra` object - even if empty because either the
#'     specified peaks are not available for that sample or because no
#'     matching `Spectra` could be found.
#'
#' @noRd
.h5_chrom_peak_spectra_sample <- function(h5_file, id, s, method, msLevel,
                                          chromPeaksMsLevel = 1L,
                                          expandRt = 0, expandMz = 0, ppm = 0,
                                          skipFilled = FALSE,
                                          peaks = integer(),
                                          chromPeakColumns = c("rt", "mz")) {
    p <- .h5_read_data(h5_file, id = id, name = "chrom_peaks",
                       ms_level = chromPeaksMsLevel, read_colnames = TRUE,
                       read_rownames = TRUE)[[1L]]
    if (length(peaks)) {
        if (is.character(peaks))
            peaks <- which(rownames(p) %in% peaks)
        p <- p[peaks, , drop = FALSE]
    }
    if (skipFilled && nrow(p)) {
        pkd <- .h5_read_data(h5_file, id = id, name = "chrom_peak_data",
                             ms_level = chromPeaksMsLevel,
                             read_colnames = TRUE)[[1L]]
        if (length(peaks))
            pkd <- pkd[peaks, , drop = FALSE]
        p <- p[!pkd$is_filled, , drop = FALSE]
    }
    if (nrow(p)) {
        if (ppm != 0)
            expandMz <- expandMz + p[, "mz"] * ppm / 1e6
        if (expandMz[1L] != 0) {
            p[, "mzmin"] <- p[, "mzmin"] - expandMz
            p[, "mzmax"] <- p[, "mzmax"] + expandMz
        }
        if (expandRt != 0) {
            p[, "rtmin"] <- p[, "rtmin"] - expandRt
            p[, "rtmax"] <- p[, "rtmax"] + expandRt
        }
        s <- filterMsLevel(s, msLevel)
        idx <- switch(
            method,
            all = .spectra_index_list(s, p, msLevel),
            closest_rt = .spectra_index_list_closest_rt(s, p, msLevel),
            closest_mz = .spectra_index_list_closest_mz(s, p, msLevel),
            largest_tic = .spectra_index_list_largest_tic(s,p,msLevel),
            largest_bpi = .spectra_index_list_largest_bpi(s,p,msLevel))
        ids <- rep(rownames(p), lengths(idx))
        s <- s[unlist(idx)]
        pk_data <- as.data.frame(p[ids, chromPeakColumns, drop = FALSE])
        pk_data$id <- ids
        colnames(pk_data) <- paste0("chrom_peak_", colnames(pk_data))
        s <- .add_spectra_data(s, pk_data)
        s
    } else
        Spectra()
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
    rn <- .h5_feature_definitions_rownames(x, ms_level)[[1L]]
    res <- do.call(
        cbind, lapply(x@sample_id, .h5_feature_values_sample,
                      hdf5_file = x@hdf5_file, ms_level = ms_level,
                      n_features = length(rtmed), method = method,
                      col_idx = col_idx, filled = filled, rtmed = rtmed))
    rownames(res) <- rn
    res
}

#' Returns the rtmin, rtmax, mzmin and mzmax for each feature depending on the
#' associated chrom peaks. For XcmsExperimentHdf5 we have to:
#' - loop over each sample
#' - load the mapping of features to chrom peaks
#' - load the chrom peak data
#' - extract the required information per feature
#' - after the loop: aggregate the information.
#'
#' We could also use `featureValues()`, but would need to call that 4 times,
#' i.e. load the data 4 times.
#'
#' For the final calculation of the region boundaries using the functions
#' defined with parameters `mzmin`, `mzmax`, `rtmin` and `rtmax` we consider
#' here, if multiple chrom peaks are assigned in a sample to a feature,
#' the min `"rtmin"`, `"mzmin"`, max `"rtmax"`, `"mzmax"` for each feature
#' in each sample (in contrast to the `.features_ms_region()` function that
#' considered all values for all chrom peaks of a feature).
#'
#' @param features `character` with the feature IDs.
#'
#' @noRd
.h5_features_ms_region <- function(x, mzmin, mzmax, rtmin, rtmax, features,
                                   ms_level = 1L) {
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = length(x@sample_id) + 1L, clear = FALSE)
    pb$tick(0)
    fids <- .h5_feature_definitions_rownames(x, ms_level)[[1L]]
    feature_idx <- unique(match(features, fids))
    if (anyNA(feature_idx))
        stop("Some of the provided feature IDs were not found in the",
             " data set.", call. = FALSE)
    feature_idx <- sort(feature_idx)
    cn <- .h5_chrom_peaks_colnames(x, ms_level)
    cn_idx <- match(c("mzmin", "mzmax", "rtmin", "rtmax"), cn)
    template <- rep(NA_real_, length(feature_idx))
    rtmin_mat <- rtmax_mat <- mzmin_mat <- mzmax_mat <-
        vector("list", length(x@sample_id))
    for (i in seq_along(x@sample_id)) {
        sample_id <- x@sample_id[i]
        sid <- paste0("/", sample_id, "/ms_", ms_level)
        fidx <- .h5_read_data(x@hdf5_file, sample_id,
                              name = "feature_to_chrom_peaks",
                              ms_level = ms_level)[[1L]]
        fidx <- fidx[fidx[, 1L] %in% feature_idx, , drop = FALSE]
        if (nrow(fidx)) {
            vals <- .h5_read_data(
                x@hdf5_file, sample_id, name = "chrom_peaks",
                ms_level = ms_level, i = fidx[, 2L], j = cn_idx,
                read_colnames = FALSE, read_rownames = FALSE)[[1L]]
            idx <- as.factor(match(fidx[, 1L], feature_idx))
            mzmin_mat[[i]] <- .h5_features_ms_region_values(
                template, vals[, 1L], idx, dups = base::min)
            mzmax_mat[[i]] <- .h5_features_ms_region_values(
                template, vals[, 2L], idx, dups = base::max)
            rtmin_mat[[i]] <- .h5_features_ms_region_values(
                template, vals[, 3L], idx, dups = base::min)
            rtmax_mat[[i]] <- .h5_features_ms_region_values(
                template, vals[, 4L], idx, dups = base::max)
        } else
            mzmin_mat[[i]] <- mzmax_mat[[i]] <- rtmin_mat[[i]] <-
                rtmax_mat[[i]] <- template
        pb$tick()
    }
    mzmin_mat <- do.call(base::cbind, mzmin_mat)
    mzmax_mat <- do.call(base::cbind, mzmax_mat)
    rtmin_mat <- do.call(base::cbind, rtmin_mat)
    rtmax_mat <- do.call(base::cbind, rtmax_mat)
    res <- cbind(mzmin = apply(mzmin_mat, 1L, mzmin, na.rm = TRUE),
                 mzmax = apply(mzmax_mat, 1L, mzmax, na.rm = TRUE),
                 rtmin = apply(rtmin_mat, 1L, rtmin, na.rm = TRUE),
                 rtmax = apply(rtmax_mat, 1L, rtmax, na.rm = TRUE))
    rownames(res) <- fids[feature_idx]
    res[features, , drop = FALSE]
}

.h5_features_ms_region_values <- function(x, vals, map, dups = base::min) {
    vl <- base::split(vals, map)
    l <- lengths(vl) > 1L
    x[as.integer(names(vl[!l]))] <- unlist(vl[!l], FALSE, FALSE)
    x[as.integer(names(vl[l]))] <- vapply(vl[l], dups, x[1L], USE.NAMES = FALSE)
    x
}

.h5_feature_definitions_rownames <- function(x, msLevel = x@features_ms_level) {
    ids <- paste0("/features/ms_", msLevel, "/feature_definitions_rownames")
    h5 <- rhdf5::H5Fopen(x@hdf5_file)
    on.exit(invisible(rhdf5::H5Fclose(h5)))
    lapply(ids, FUN = rhdf5::h5read, file = h5, drop = TRUE)
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
#' peaks. Chromatographic peaks are added/processed by sample reading the
#' respective information from the hdf5 file. Features are added by processing
#' the feature and chrom peak IDs of the selected chrom peaks.
#'
#'
#' @param x `XcmsExperimentHdf5` for one subset/chunk of data from which the
#'     data should be extracted
#'
#' @param index `integer` with the index of the current subset `x` in the *full*
#'     data set.
#'
#' @param ms_level `integer(1)` with the MS level.
#' @noRd
.h5_x_chromatograms <- function(x, index = seq_along(x), aggregationFun = "sum",
                               ms_level = 1L, mz, rt, isolationWindow = NULL,
                               chromPeaks = "any", chunkSize = 2L,
                               return.type = "XChromatograms",
                               BPPARAM = bpparam()) {
    message("Extracting chromatographic data")
    chr <- as(.mse_chromatogram(
        as(x, "MsExperiment"), rt = rt, mz = mz, msLevel = ms_level,
        aggregationFun = aggregationFun, isolationWindow = isolationWindow,
        chunkSize = chunkSize, BPPARAM = BPPARAM), return.type)
    if (return.type == "MChromatograms" || chromPeaks == "none")
        return(chr)
    message("Processing chromatographic peaks")
    js <- seq_len(nrow(chr))
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = length(x), clear = FALSE)
    has_features <- hasFeatures(x, ms_level)
    if (has_features) {
        fd <- .h5_read_data(x@hdf5_file, "features", "feature_definitions",
                            read_rownames = TRUE, ms_level = ms_level)[[1L]]
        f2p <- vector("list", 1000) # initialize with an educated guess;
        cnt <- 1L
    }
    mat <- chr@.Data
    slot(chr, ".Data", check = FALSE) <- matrix(ncol = ncol(chr),
                                                nrow = nrow(chr))
    for (i in seq_along(x)) { # iterate over samples
        cp <- .h5_read_data(x@hdf5_file, x@sample_id[i], "chrom_peaks",
                            ms_level = ms_level, read_colnames = TRUE,
                            read_rownames = TRUE)[[1L]]
        cd <- .h5_read_data(x@hdf5_file, x@sample_id[i], "chrom_peak_data",
                            ms_level = ms_level, read_colnames = TRUE,
                            read_rownames = FALSE)[[1L]]
        if (has_features)
            fidx <- .h5_read_data(x@hdf5_file, x@sample_id[i],
                                  "feature_to_chrom_peaks",
                                  ms_level = ms_level)[[1L]]
        for (j in js) { # iterate over ranges/EICs
            idx <- which(.is_chrom_peaks_within_mz_rt(
                cp, rt[j, ], mz[j, ], type = chromPeaks), useNames = FALSE)
            li <- length(idx)
            a <- cbind(cp[idx, , drop = FALSE], sample = rep(i, li))
            b <- extractROWS(cd, idx)
            b$ms_level <- rep(ms_level, li)
            attr(b, "row.names") <- rownames(a)
            tmp <- mat[j, i][[1L]]
            slot(tmp, "chromPeaks", check = FALSE) <- a
            slot(tmp, "chromPeakData", check = FALSE) <- as(b, "DataFrame")
            mat[j, i][[1L]] <- tmp # this does not seem to cause copying
            ## Add mapping of features and chrom peaks for that sample/EIC
            if (li && has_features) {
                is_feature <- fidx[, 2L] %in% idx
                if (any(is_feature)) {
                    f2p[[cnt]] <- cbind(rownames(fd)[fidx[is_feature, 1L]],
                                        rownames(cp)[fidx[is_feature, 2L]],
                                        rep(as.character(j), sum(is_feature)))
                    cnt <- cnt + 1L
                }
            }
        }
        pb$tick()
    }
    slot(chr, ".Data", check = FALSE) <- mat
    if (has_features) {
        message("Processing features")
        ## Define the featureDefinitions and assign chrom peaks to each,
        ## matching the same `"row"` (!).
        f2p <- do.call(base::rbind, f2p) # mapping of features and chrom peaks
        cp <- chromPeaks(chr)            # chrom peaks to calculate index
        cp_row_id <- paste0(rownames(cp), ".", cp[, "row"])
        ft_cp_row_id <- paste0(f2p[, 2], ".", f2p[, 3])
        ft_row_id <- paste0(f2p[, 1], ".", f2p[, 3])
        peakidx <- lapply(split(ft_cp_row_id, ft_row_id),
                          base::match, table = cp_row_id)
        id <- strsplit(names(peakidx), ".", fixed = TRUE)
        fts <- fd[vapply(id, `[`, i = 1L, FUN.VALUE = NA_character_), ]
        fts$peakidx <- unname(peakidx)
        fts$row <- vapply(id, function(z) as.integer(z[2L]), NA_integer_)
        fts$ms_level <- ms_level
        slot(chr, "featureDefinitions", check = FALSE) <-
            DataFrame(extractROWS(fts, order(fts$row)))
    }
    slot(chr, ".processHistory", check = FALSE) <- x@processHistory
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
                                        ppm = 0, type = "any",
                                        sample_index = integer()) {
    read_colnames <- read_colnames || length(rt) > 0 || length(mz) > 0
    d <- .h5_read_matrix2(
        name, h5, index, read_colnames, read_rownames, rownames)
    if (length(rt) | length(mz))
        d <- d[.is_chrom_peaks_within_mz_rt(d, rt = rt, mz = mz,
                                            ppm = ppm, type = type), ,
               drop = FALSE]
    ## If sample_index is provided add a column "sample" with the index.
    if (length(sample_index))
        d <- cbind(d, sample = rep(sample_index[name], nrow(d)))
    d
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
        attr(d, "row.names") <- rhdf5::h5read(h5, rownames, drop = TRUE)
    if (is.null(index[[1L]]))
        d
    else d[index[[1L]], , drop = FALSE]
}

#' Function to read the chromPeakData `data.frame` for one sample.
#'
#' @param peaks optional `character` with the chrom peak IDs for which the
#'     data should be extracted.
#'
#' @param ms_levels optional **named** `integer` with the MS levels of all
#'     imported data sets. At least one of the `names(ms_levels)` should match
#'     param `name`. If provided, a column `$ms_level` is added to the result.
#'
#' @noRd
.h5_read_chrom_peak_data <- function(name, h5, index = list(NULL, NULL),
                                     read_rownames = FALSE, peaks = character(),
                                     ms_levels = integer(), ...) {
    cd <- .h5_read_data_frame(
        name, h5, read_rownames = read_rownames || length(peaks) > 0,
        index = index, rownames = sub("_data", "s_rownames", name))
    if (length(peaks)) {
        idx <- match(peaks, rownames(cd))
        cd <- cd[idx[!is.na(idx)], , drop = FALSE]
        if (!read_rownames)
            rownames(cd) <- NULL
    }
    if (length(ms_levels))
        cd$ms_level <- rep(unname(ms_levels[name]), nrow(cd))
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
    mc <- .h5_mod_count(h5)
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

.h5_mod_count <- function(h5) {
    rhdf5::h5read(h5, "/header/modcount")[1L]
}

#' Every writing operation should increate the "mod count", i.e. the
#' count of data modifications. This function increases the mod count by +1
#' and returns this value.
#'
#' @noRd
.h5_increment_mod_count <- function(h5) {
    mc <- .h5_mod_count(h5) + 1L
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
        rn <- rownames(x)
        if (is.null(rn)) rn <- character()
        dn <- paste0(name, "_rownames")
        if (replace && rhdf5::H5Lexists(h5, dn))
            rhdf5::h5delete(h5, dn)
        rhdf5::h5write(rn, h5, name = dn, level = level, createnewfile = FALSE)
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
        rn <- rownames(x)
        if (is.null(rn)) rn <- character()
        dn <- paste0(name, "_rownames")
        if (replace && rhdf5::H5Lexists(h5, dn))
            rhdf5::h5delete(h5, dn)
        rhdf5::h5write(rn, h5, name = dn, level = level, createnewfile = FALSE)
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
