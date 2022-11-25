#' Adding chrom peaks to an XcmsExperiment. If chrom peaks are present they
#' will be added. Rownames will be fixed/added too.
#'
#' @noRd
.mse_add_chrom_peaks <- function(object, pks, pkd) {
    ## Keeping also old chromPeaks.
    nr <- nrow(object@chromPeaks)
    if (nr) {
        ## Keep the original chrom peak names.
        idx <- max(as.integer(sub("CP", "", rownames(object@chromPeaks))))
        rownames(pks) <- rownames(pkd) <- .featureIDs(nrow(pks), prefix = "CP",
                                                      from = idx + 1L)
        object@chromPeaks <- suppressWarnings(
            rbindFill(object@chromPeaks, pks))
        object@chromPeakData <- suppressWarnings(
            rbindFill(object@chromPeakData, pkd))
    } else {
        rownames(pks) <- rownames(pkd) <- .featureIDs(nrow(pks), prefix = "CP")
        object@chromPeaks <- pks
        object@chromPeakData <- pkd
    }
    object
}

.mse_valid_chrom_peaks <- function(x) {
    if (!all(.REQ_PEAKS_COLS %in% colnames(x)))
        "One or more required columns missing in 'chromPeaks'"
    else NULL
}

.mse_valid_chrom_peak_data <- function(x) {
    if (!all(c("ms_level", "is_filled") %in% colnames(x)))
        "One or more required columns missing in 'chromPeakData'"
    else NULL
}

.mse_valid_feature_def <- function(x) {
    msg <- NULL
    if (!is.data.frame(x))
        msg <- "featureDefinitions has to be a data.frame"
    if (!all(.REQ_PEAKG_COLS %in% colnames(x)))
        msg <- "featureDefinitions lacks one or more required columns"
    msg
}

.mse_same_rownames <- function(a, b) {
    if (!all(rownames(a) == rownames(b)))
        "Row names of 'chromPeaks' and 'chromPeakData' don't match"
    else NULL
}

#' @title Subsetting an XcmsExperiment
#'
#' @description
#'
#' Subsetting the xcms results by sample. In contrast to the XCMSnExp we allow
#' here a little more flexibility:
#' - keepChromPeaks: do not automatically drop the chrom peaks.
#' - keepAdjustedRtime: eventually keep adjusted retention times.
#' - keepFeatures: keep features.
#'
#' Also, we allow subsetting in arbitrary order, re-assigning chrom peaks etc.
#'
#' @param ignoreHistory `logical(1)` whether process history should not be
#'     updated. Setting to `TRUE` would be faster.
#'
#' @note Unit tests are in test_XcmsExperiment.R!
#'
#' @author Johannes Rainer
#'
#' @noRd
.subset_xcms_experiment <- function(x, i = integer(),
                                    keepChromPeaks = TRUE,
                                    keepAdjustedRtime = FALSE,
                                    keepFeatures = FALSE,
                                    ignoreHistory = FALSE,
                                    keepSampleIndex = FALSE, ...) {
    ## if (!length(i))
    ##     return(x)
    i <- MsCoreUtils::i2index(i, length(x))
    ## This is a special case that would make life (=performance) miserable
    if (length(i) != length(unique(i)))
        stop("Duplicated indices are not (yet) supported for ",
             "'[,XcmsExperiment'", call. = FALSE)
    if (!keepAdjustedRtime && hasAdjustedRtime(x)) {
        svs <- unique(c(spectraVariables(object@spectra), "mz", "intensity"))
        object@spectra <- selectSpectraVariables(
            object@spectra, svs[svs != "rtime_adjusted"])
    }
    if (keepFeatures && hasFeatures(x))
        stop("Not implemented yet", call. = FALSE)
    if (hasChromPeaks(x)) {
        if (keepChromPeaks) {
            keep <- x@chromPeaks[, "sample"] %in% i
            x@chromPeaks <- x@chromPeaks[keep, , drop = FALSE]
            x@chromPeakData <- x@chromPeakData[keep, , drop = FALSE]
            if (!keepSampleIndex)
                x@chromPeaks[, "sample"] <- match(x@chromPeaks[, "sample"], i)
        } else {
            x@chromPeaks <- .empty_chrom_peaks()
            x@chromPeakData <- data.frame(ms_level = integer(),
                                          is_filled = logical())
            if (!ignoreHistory)
                x@processHistory <- dropProcessHistoriesList(
                    x@processHistory,
                    type = c(.PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.GROUPING,
                             .PROCSTEP.PEAK.FILLING, .PROCSTEP.CALIBRATION,
                             .PROCSTEP.PEAK.REFINEMENT))
        }
    }
    getMethod("[", "MsExperiment")(x, i = i)
}

#' @param x `chromPeaks` `matrix`.
#'
#' @param param The parameter object with the settings for the peak grouping.
#'
#' @author Johannes Rainer
#'
#' @noRd
.xmse_group_cpeaks <- function(cp, param, index = seq_len(nrow(cp))) {
    pclass <- class(param)[1L]
    if (!nrow(cp))
        return(.empty_feature_definitions())
    switch(
        pclass,
        PeakDensityParam = {
            do_groupChromPeaks_density(
                cp, sampleGroups = sampleGroups(param), bw = bw(param),
                minFraction = minFraction(param),
                minSamples = minSamples(param), binSize = binSize(param),
                maxFeatures = maxFeatures(param), index = index)
        },
        MzClustParam = {
            tmp <- do_groupPeaks_mzClust(
                cp, sampleGroups = sampleGroups(param), ppm = ppm(param),
                absMz = absMz(param), minFraction = minFraction(param),
                minSamples = minSamples(param))
            res <- as.data.frame(tmp$featureDefinitions)
            res$peakidx <- tmp$peakIndex
            res
        },
        NearestPeaksParam = {
            tmp <- do_groupChromPeaks_nearest(
                cp, sampleGroups = sampleGroups(param),
                mzVsRtBalance = mzVsRtBalance(param),
                absMz = absMz(param), absRt = absRt(param), kNN = kNN(param))
            res <- as.data.frame(tmp$featureDefinitions)
            res$peakidx <- tmp$peakIndex
            res
        },
        stop("No correspondence analysis method for '", pclass,
             "' available", call. = FALSE))
}

#' @param x peak matrix (columns `"mz"` and `"intensity"`).
#'
#' @param mzr m/z range within which to aggregate intensities
#'
#' @param INTFUN function to aggregate the intensities.
#'
#' @param ... additional parameters for `FUN`.
#'
#' @note
#'
#' `matrixStats` `sum2`, `colSums2` etc could be nice alternatives that don't
#' subset the matrix because they iterate on the indices. Somehow `sum` seems
#' however to be faster than `sum2` - even on subsetting etc.
#'
#' @author Johannes Rainer
#'
#' @noRd
.aggregate_intensities <- function(x, mzr = numeric(), INTFUN = sumi, ...) {
    ## x could also be just a numeric vector - to support Chromatogram...
    if (length(mzr))
        vals <- x[between(x[, "mz"], mzr), "intensity"]
    else vals <- x[, "intensity"]
    INTFUN(vals, ...)
}

#' Filter a peaks matrix based on m/z value range. Returns a `matrix` with 0
#' rows if no value is within the mz range.
#'
#' @noRd
.pmat_filter_mz <- function(x, mzr = numeric()) {
    if (!length(mzr))
        return(x)
    x[between(x[, "mz"], mzr), , drop = FALSE]
}

#' @title Sum MS Intensity Values
#'
#' @description
#'
#' `sumi` sums MS intensity values ignoring missing values but returning
#' `NA_real_` in case all intensity values are `NA`. This is in contrast to
#' `sum(x, na.rm = TRUE)` that, while ignoring also missing values, returns
#' `0` if all values in `x` are `NA`.
#'
#' @param x `numeric` with the intensity values to sum.
#'
#' @return `numeric(1)` with the sum of `x` or `NA_real_` if **all** values in
#'     `x` are `NA`.
#'
#' @author Johannes Rainer
#'
#' @noRd
sumi <- function(x) {
    ## Note: any(!is.na(x)) would be faster, but ONLY if there are NAs,
    ## otherwise this version is slightly faster.
    if (all(is.na(x)))
        return(NA_real_)
    sum(x, na.rm = TRUE)
}

#' Same as functions-Chromatogram.R::.chrom_merge_neighboring_peaks, but this
#' takes a list of peak matrices as input and does not require creation of a
#' Chromatogram object.
#'
#' @param x `list` with `numeric` matrices with columns `mz` and `intensity`
#'     representing the MS data. This has to be from a single sample and from
#'     a single MS level.
#'
#' @param rt `numeric` with the retention times for `x`. HAS TO BE ADJUSTED RT
#'
#' @note
#'
#' The m/z range of the new merged peak considers also `expandMz` and `ppm`!
#' The m/z range of this version represents the m/z range of the intensity
#' actually considered for the integration. This is in contrast to
#' functions-Chromatogram.R/.chrom_merge_neighboring_peaks that simply uses
#' the m/z range of all candidates that should be evaluated for merging.
#'
#' @noRd
.merge_neighboring_peak_candidates <- function(x, rt, pks, pkd, minProp = 0.75,
                                                expandMz = 0, ppm = 10,
                                                diffRt = 0) {
    if (!length(pks) || nrow(pks) < 2)
        return(list(chromPeaks = pks, chromPeakData = pkd))
    idx <- order(pks[, "rtmin"])
    pks <- pks[idx, , drop = FALSE]
    pkd <- pkd[idx, ]
    rownames(pkd) <- NULL
    pks_new <- pks
    pks_new[ , ] <- NA_real_
    rownames(pks_new) <- rep(NA_character_, nrow(pks))
    pks_new[1L, ] <- pks[1L, ]
    rownames(pks_new)[1L] <- rownames(pks)[1L]
    current_peak <- 1L # point always to the current *new* (merged) peak.
    drop_cols <- !(colnames(pks_new) %in% c("mz", "mzmin", "mzmax", "rt",
                                            "rtmin", "rtmax", "into",
                                            "maxo", "sn", "sample"))
    ## Similar to the original code we define the "expanded m/r range" for the
    ## full set of peaks based on the m/z min and max of **all** candidate peaks
    ## Adjusting the m/z range for each tested candidate peak individually would
    ## eventually result in wrong intensity estimation. The current approach
    ## is more greedy, but, using *reasonable* settings it is supposed ot be
    ## correct.
    ## The reported m/z range for merged candidates represents the full m/z
    ## range of all intensities considered in the calculation of the "into".
    ppme <- ppm * 1e-6
    mzr <- range(pks[, c("mzmin", "mzmax")])
    mzr_e <- c(mzr[1L] - expandMz - ppme * mzr[1L],
               mzr[2L] + expandMz + ppme * mzr[2L])
    ## pre-selecting the peaksData x *might* speed up things?
    rtidx <- which(rt >= min(pks[, "rtmin"]) & rt <= max(pks[, "rtmax"]))
    rt <- rt[rtidx]
    x <- lapply(x[rtidx], .pmat_filter_mz, mzr = mzr_e)
    for (i in 2:nrow(pks)) {
        if ((pks[i, "rtmin"] - pks_new[current_peak, "rtmax"]) < diffRt) {
            ## skip if second peak contained within first
            if (pks[i, "rtmin"] >= pks_new[current_peak, "rtmin"] &
                pks[i, "rtmax"] <= pks_new[current_peak, "rtmax"])
                next
            rt_mid <- (pks[i, "rtmin"] + pks_new[current_peak, "rtmax"]) / 2
            ## If rt_mid is NOT between the peaks, take the midpoint between
            ## the apexes instead.
            apexes <- range(c(pks[i, "rt"], pks_new[current_peak, "rt"]))
            if (rt_mid < apexes[1L] || rt_mid > apexes[2L])
                rt_mid <- sum(apexes) / 2
            ## Calculate the mean of the 3 data points closest to rt_mid.
            rt_idx <- order(abs(rt - rt_mid))[1:3]
            mid_vals <- vapply(x[rt_idx], .aggregate_intensities, numeric(1))
            if (!all(is.na(mid_vals)) &&
                mean(mid_vals, na.rm = TRUE) >
                min(pks_new[current_peak, "maxo"], pks[i, "maxo"]) * minProp) {
                ## Merge the existing peak with the new one, re-calculating the
                ## intensity and the m/z range from the original data.
                pks_new[current_peak, drop_cols] <- NA_real_
                rtmin <- pks_new[current_peak, "rtmin"]
                rtmax <- max(pks[i, "rtmax"], pks_new[current_peak, "rtmax"])
                idx <- which(rt >= rtmin & rt <= rtmax)
                peak_width <- (rtmax - rtmin) / (idx[length(idx)] - idx[1L])
                pkm <- do.call(rbind, x[idx])
                pks_new[current_peak, c("rtmax", "mzmin", "mzmax", "into")] <-
                    c(rtmax, range(pkm[, "mz"], na.rm = TRUE),
                      sum(pkm[, "intensity"], na.rm = TRUE) * peak_width)
                if (pks[i, "maxo"] > pks_new[current_peak, "maxo"]) {
                    pks_new[current_peak, c("mz", "rt", "maxo", "sn")] <-
                        pks[i, c("mz", "rt", "maxo", "sn")]
                    pkd[current_peak, ] <- pkd[i, ] # replace peak data with new
                }
                rownames(pks_new)[current_peak] <- NA_character_
            } else {
                current_peak <- current_peak + 1
                pks_new[current_peak, ] <- pks[i, ]
                rownames(pks_new)[current_peak] <- rownames(pks)[i]
                pkd[current_peak, ] <- pkd[i, ]
            }
        } else {
            current_peak <- current_peak + 1
            pks_new[current_peak, ] <- pks[i, ]
            rownames(pks_new)[current_peak] <- rownames(pks)[i]
            pkd[current_peak, ] <- pkd[i, ]
        }
    }
    keep <- which(!is.na(pks_new[, "rt"]))
    list(chromPeaks = pks_new[keep, , drop = FALSE],
         chromPeakData = pkd[keep, ])
}

#' similar to functions-XCMSnExp.R/.merge_neigboring_peaks but works on only
#' low level data from one file/sample. All data is expected to be from a
#' single sample and the correct (single) MS level.
#'
#' @param x `list` of peak matrices (as returned by `Spectra::peaksData`)
#'
#' @param pks `matrix` with chromatographic peaks (`chromPeaks` for one MS)
#'
#' @param pkd `data.frame` (`chromPeakData`)
#'
#' @param rt `numeric` with the retention times.
#'
#' @return `list` with the peaks matrix and peaks data containing all peaks as
#'     well as the merged peaks.
#'
#' @noRd
.merge_neighboring_peaks2 <- function(x, pks, pkd, rt, expandRt = 2,
                                      expandMz = 0, ppm = 10, minProp = 0.75) {
    cands <- .define_merge_candidates(pks, expandMz, ppm, expandRt)
    if (!length(cands))
        return(list(chromPeaks = pks, chromPeakData = pkd))
    cands <- cands[[2L]]
    pks_new <- pkd_new <- vector("list", length(cands))
    for (i in seq_along(cands)) {
        res <- .merge_neighboring_peak_candidates(
            x, rt = rt, pks[cands[[i]], , drop = FALSE],
            pkd[cands[[i]], , drop = FALSE], diffRt = 2 * expandRt,
            minProp = minProp, expandMz = expandMz, ppm = ppm)
        pks_new[[i]] <- res$chromPeaks
        pkd_new[[i]] <- res$chromPeakData
    }
    pks_new <- do.call(rbind, pks_new)
    pkd_new <- do.call(rbind, pkd_new)
    ## drop peaks that were candidates, but that were not returned (i.e.
    ## were either merged or dropped)
    keep <- !(rownames(pks) %in% setdiff(unlist(cands, use.names = FALSE),
                                         rownames(pks_new)))
    pks <- pks[keep, , drop = FALSE]
    pkd <- pkd[keep, ]
    ## add merged peaks
    news <- is.na(rownames(pks_new))
    if (any(news)) {
        pks <- rbind(pks, pks_new[news, , drop = FALSE])
        pkd <- rbind(pkd, pkd_new[news, ])
    }
    list(chromPeaks = pks, chromPeakData = pkd, npeaks = nrow(pks))
}

#' @param x `XcmsExperiment` with potentially multiple files, samples.
#'
#' @noRd
.xmse_merge_neighboring_peaks <- function(x, msLevel = 1L, expandRt = 2,
                                          expandMz = 0, ppm = 10,
                                          minProp = 0.75, BPPARAM = bpparam()) {
    keep <- msLevel(spectra(x)) == msLevel
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    f_peaks <- factor(chromPeaks(x, msLevel = msLevel)[, "sample"],
                      levels = levels(f))
    res <- bpmapply(
        .merge_neighboring_peaks2,
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel)), f),
        split.data.frame(chromPeaks(x, msLevel = msLevel), f = f_peaks),
        split.data.frame(chromPeakData(
            x, msLevel = msLevel, return.type = "data.frame"), f = f_peaks),
        split(rt, f),
        MoreArgs = list(expandRt = expandRt, expandMz = expandMz,
                        ppm = ppm, minProp = minProp),
        SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    list(chromPeaks = do.call(rbind, lapply(res, `[[`, 1L)),
         chromPeakData = do.call(rbind, lapply(res, `[[`, 2L)),
         npeaks = vapply(res, `[[`, i = 3L, integer(1)))
}

#' @param x `XcmsExperiment` (multiple files)
#'
#' @return `logical` or length equal `nrow(chromPeaks(x))`.
#'
#' @noRd
.xmse_filter_peaks_intensities <- function(x, nValues = 1L, threshold = 0,
                                           msLevel = 1L, BPPARAM = bpparam()) {
    keep <- msLevel(spectra(x)) == msLevel
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    f_peaks <- factor(chromPeaks(x)[, "sample"], levels = levels(f))
    res <- bpmapply(
        function(x, rt, pks, msLevels, nValues, threshold, msLevel) {
            if (!length(pks)) return(logical())
            vapply(seq_len(nrow(pks)), function(i) {
                if (msLevels[i] != msLevel)
                    FALSE
                else {
                    rtidx <- between(rt , pks[i, c("rtmin", "rtmax")])
                    vals <- vapply(x[rtidx], .aggregate_intensities,
                                   mzr = pks[i, c("mzmin", "mzmax")],
                                   numeric(1))
                    sum(vals >= threshold, na.rm = TRUE) >= nValues
                }
            }, logical(1))
        },
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel)), f),
        split(rt, f),
        split.data.frame(chromPeaks(x), f = f_peaks),
        split(chromPeakData(x)$ms_level, f = f_peaks),
        MoreArgs = list(nValues = nValues, threshold = threshold,
                        msLevel = msLevel),
        USE.NAMES = FALSE, BPPARAM = BPPARAM)
    unlist(res, use.names = FALSE)
}

#' Apply any function `FUN` to chunks of an `XcmsExperiment`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.xmse_apply_chunks <- function(x, FUN, ..., keepChromPeaks = TRUE,
                               keepAdjustedRtime = FALSE, keepFeatures = FALSE,
                               ignoreHistory = FALSE, keepSampleIndex = FALSE,
                               chunkSize = 1L) {
    idx <- seq_along(x)
    chunks <- split(idx, ceiling(idx / chunkSize))
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = length(chunks), clear = FALSE)
    pb$tick(0)
    lapply(chunks, function(z, ...) {
        suppressMessages(
            res <- FUN(.subset_xcms_experiment(
                x, i = z, keepChromPeaks = keepChromPeaks,
                keepAdjustedRtime = keepAdjustedRtime,
                keepFeatures = keepFeatures, ignoreHistory = ignoreHistory,
                keepSampleIndex = keepSampleIndex), ...)
        )
        pb$tick()
        res
    }, ...)
}

#' Perform peak integration on provided m/z - RT regions. Can be called for
#' gap filling or manual definition of chrom peaks. Requires the data (as
#' `XcmsExperiment` or `MsExperiment`, the definition of regions and the
#' integration function. This function simply splits the input object
#' and calls the (provided) integration function.
#'
#' @param x `XcmsExperiment` with data from potentially multiple files.
#'
#' @param pal `list` of peak area matrices defining (for each individual sample
#'     in `x`) the MS area from which the signal should be integrated. Names
#'     should represent the file/sample index!
#'
#' @param msLevel `integer(1)` with the MS level on which to integrate the data.
#'
#' @param intFun function to be used for the integration.
#'
#' @param mzCenterFun function to calculate the m/z value
#'
#' @return `matrix` with the newly integrated peaks (NO chromPeakData!).
#'
#' @noRd
.xmse_integrate_chrom_peaks <- function(x, pal, msLevel = 1L,
                                        intFun = .chrom_peak_intensity_centWave,
                                        mzCenterFun = "mzCenter.wMean",
                                        param = MatchedFilterParam(),
                                        BPPARAM = bpparam()) {
    keep <- which(msLevel(spectra(x)) == msLevel)
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    cn <- colnames(chromPeaks(x))
    res <- bpmapply(split(peaksData(filterMsLevel(spectra(x), msLevel)), f),
                    split(rt, f),
                    pal,
                    as.integer(names(pal)),
                    FUN = intFun,
                    MoreArgs = list(mzCenterFun = mzCenterFun, cn = cn,
                                    param = param),
                    SIMPLIFY = FALSE, USE.NAMES = FALSE, BPPARAM = BPPARAM)
    do.call(rbind, res)
}

#' Integrates MS intensities for a chromatographic peak. This is equivalent
#' to the `.getChromPeakData` function.
#'
#' @param x `list` of peak matrices (from a single MS level and from a single
#'     file/sample).
#'
#' @param rt retention time for each peak matrix.
#'
#' @param peakArea `matrix` defining the chrom peak area.
#'
#' @author Johannes Rainer
#'
#' @noRd
.chrom_peak_intensity_centWave <- function(x, rt, peakArea,
                                           mzCenterFun = "weighted.mean",
                                           sampleIndex = integer(),
                                           cn = character(), ...) {
    res <- matrix(NA_real_, ncol = length(cn), nrow = nrow(peakArea))
    rownames(res) <- rownames(peakArea)
    colnames(res) <- cn
    res[, "sample"] <- sampleIndex
    res[, c("rtmin", "rtmax", "mzmin", "mzmax")] <-
        peakArea[, c("rtmin", "rtmax", "mzmin", "mzmax")]
    for (i in seq_len(nrow(res))) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        keep <- which(between(rt, rtr))
        if (length(keep)) {
            xsub <- lapply(x[keep], .pmat_filter_mz,
                           mzr = peakArea[i, c("mzmin", "mzmax")])
            mat <- do.call(rbind, xsub)
            if (nrow(mat)) {
                ## can have 0, 1 or x values per rt; repeat rt accordingly
                rts <- rep(rt[keep], vapply(xsub, nrow, integer(1L)))
                maxi <- which.max(mat[, 2L])[1L]
                mmz <- do.call(mzCenterFun, list(mat[, 1L], mat[, 2L]))
                if (is.na(mmz)) mmz <- mat[maxi, 1L]
                res[i, c("rt", "mz", "maxo", "into")] <- c(
                    rts[maxi], mmz, mat[maxi, 2L],
                    sum(mat[, 2L], na.rm = TRUE) *
                    ((rtr[2L] - rtr[1L]) / max(1L, (length(keep) - 1L)))
                )
            }
        }
    }
    res[!is.na(res[, "maxo"]), , drop = FALSE]
}

.chrom_peak_intensity_matchedFilter <- function(x, rt, peakArea,
                                                mzCenterFun = "weighted.mean",
                                                sampleIndex = integer(),
                                                cn = character(),
                                                param = MatchedFilterParam(),
                                                ...) {
    ## Need to calculate the profile matrix and then work on that.
    stop("Not yet implemented")
    pmat <- .peaksdata_profmat(x)
}

.chrom_peak_intensity_msw <- function(x, rt, peakArea,
                                      mzCenterFun = "weighted.mean",
                                      sampleIndex = integer(),
                                      cn = character(), ...) {
    stop("Not yet implemented")
}

.xmse_process_history <- function(x, type = character(), msLevel = integer()) {
    if (!length(msLevel) && !length(type))
        return(x@processHistory)
    if (length(msLevel))
        keep_msl <- vapply(
            x@processHistory, function(z) msLevel(z) %in% msLevel, logical(1))
    else keep_msl <- rep(TRUE, length(x@processHistory))
    if (length(type))
        keep_t <- vapply(
            x@processHistory, function(z) processType(z) %in% type, logical(1))
    else keep_t <- rep(TRUE, length(x@processHistory))
    x@processHistory[keep_msl & keep_t]
}

.history2fill_fun <- function(x = list()) {
    if (!length(x))
        cl <- "CentWaveParam"
    else cl <- class(x[[1L]]@param)[1L]
    switch(cl,
           MatchedFilterParam = .chrom_peak_intensity_matchedFilter,
           MSWParam = .chrom_peak_intensity_msw,
           .chrom_peak_intensity_centWave)
}
