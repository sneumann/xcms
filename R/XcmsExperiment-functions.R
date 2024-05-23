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
    i <- i2index(i, length(x))
    if (any(i < 0)) {
        if (all(i < 0))
            i <- seq_along(x)[i]
        else stop("Mixing positive and negative indices is not supported.")
    }
    ## This is a special case that would make life (=performance) miserable
    if (length(i) != length(unique(i)))
        stop("Duplicated indices are not (yet) supported for ",
             "'[,XcmsExperiment'", call. = FALSE)
    drop <- character()
    if (!keepAdjustedRtime && hasAdjustedRtime(x)) {
        svs <- unique(c(spectraVariables(x@spectra), "mz", "intensity"))
        x@spectra <- selectSpectraVariables(
            x@spectra, svs[svs != "rtime_adjusted"])
        drop <- c(drop, .PROCSTEP.RTIME.CORRECTION)
    }
    if (!keepFeatures && hasFeatures(x)) {
        x@featureDefinitions <- .empty_feature_definitions()
        drop <- c(drop, .PROCSTEP.PEAK.GROUPING)
    }
    if (hasChromPeaks(x)) {
        if (keepChromPeaks) {
            x <- .filter_chrom_peaks(x, which(x@chromPeaks[, "sample"] %in% i))
            if (!keepSampleIndex)
                x@chromPeaks[, "sample"] <- match(x@chromPeaks[, "sample"], i)
        } else {
            x@chromPeaks <- .empty_chrom_peaks()
            x@chromPeakData <- data.frame(ms_level = integer(),
                                          is_filled = logical())
            drop <- c(drop, .PROCSTEP.PEAK.DETECTION, .PROCSTEP.PEAK.FILLING,
                      .PROCSTEP.CALIBRATION, .PROCSTEP.PEAK.REFINEMENT)
        }
    }
    if (!ignoreHistory && length(drop))
        x@processHistory <- dropProcessHistoriesList(
            x@processHistory, type = drop)
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
                maxFeatures = maxFeatures(param), ppm = ppm(param),
                index = index)
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
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel),
                        f = factor()), f),
        split.data.frame(chromPeaks(x, msLevel = msLevel), f = f_peaks),
        split.data.frame(.chromPeakData(x, msLevel = msLevel), f = f_peaks),
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
    f_peaks <- factor(.chromPeaks(x)[, "sample"], levels = levels(f))
    res <- bpmapply(
        function(x, rt, pks, msLevels, nValues, threshold, msLevel) {
            if (!length(pks)) return(logical())
            rtc <- c("rtmin", "rtmax")
            mzc <- c("mzmin", "mzmax")
            vapply(seq_len(nrow(pks)), function(i) {
                if (msLevels[i] != msLevel)
                    FALSE
                else {
                    rtidx <- between(rt , pks[i, rtc])
                    vals <- vapply(x[rtidx], .aggregate_intensities,
                                   mzr = pks[i, mzc],
                                   numeric(1))
                    sum(vals >= threshold, na.rm = TRUE) >= nValues
                }
            }, logical(1))
        },
        split(peaksData(filterMsLevel(spectra(x), msLevel = msLevel),
                        f = factor()), f),
        split(rt, f),
        split.data.frame(.chromPeaks(x), f = f_peaks),
        split(.chromPeakData(x)$ms_level, f = f_peaks),
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
                                        BPPARAM = bpparam(), ...) {
    keep <- which(msLevel(spectra(x)) == msLevel)
    f <- as.factor(fromFile(x)[keep])
    if (hasAdjustedRtime(x)) rt <- spectra(x)$rtime_adjusted[keep]
    else rt <- rtime(spectra(x))[keep]
    cn <- colnames(.chromPeaks(x))
    res <- bpmapply(split(peaksData(filterMsLevel(spectra(x), msLevel),
                                    f = factor()), f),
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
            ## length of xsub is the number of spectra, the number of peaks can
            ## however be 0 if no peak was found. Maybe we should/need to
            ## consider adding 0 or NA intensity for those.
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
                if ("beta_cor" %in% cn)
                    res[i, c("beta_cor", "beta_snr")] <- .get_beta_values(
                        mat[, 2L], rts)
            }
        }
    }
    res[!is.na(res[, "maxo"]), , drop = FALSE]
}

#' Difference to the original code is that the weighted mean is also calculated
#' if some of the peak intensities in the profile matrix are 0
#'
#' Note also that the rt of a peak is estimated differently, i.e. not from the
#' gaussian fitted curve, but simply the rt of the largest signal.
#'
#' @noRd
.chrom_peak_intensity_matchedFilter <- function(x, rt, peakArea,
                                                mzCenterFun = "weighted.mean",
                                                sampleIndex = integer(),
                                                cn = character(),
                                                param = MatchedFilterParam(),
                                                ...) {
    res <- matrix(NA_real_, ncol = length(cn), nrow = nrow(peakArea))
    rownames(res) <- rownames(peakArea)
    colnames(res) <- cn
    res[, "sample"] <- sampleIndex
    res[, c("rtmin", "rtmax", "mzmin", "mzmax")] <-
        peakArea[, c("rtmin", "rtmax", "mzmin", "mzmax")]
    basespc <- NULL
    if (length(distance(param)) > 0) {
        mzr <- range(vapply(x, function(z) range(z[, 1L]), numeric(2)))
        mass <- seq(floor(mzr[1L] / binSize(param)) * binSize(param),
                    ceiling(mzr[2L] / binSize(param)) * binSize(param),
                    by = binSize(param))
        bin_size <- (max(mass) - min(mass)) / (length(mass) - 1)
        basespc <- distance(param) * bin_size
    }
    pmat <- .peaksdata_profmat(
        x, method = .impute2method(param), step = binSize(param),
        baselevel = baseValue(param), basespace = basespc, baseValue = 0,
        returnBreaks = TRUE)
    brks <- pmat$breaks
    pmat <- pmat$profMat
    bin_size <- diff(brks[1:2])
    bin_half <- bin_size / 2
    mzs <- brks[-length(brks)] + bin_half ## midpoint for the breaks
    for (i in seq_len(nrow(res))) {
        rtr <- peakArea[i, c("rtmin", "rtmax")]
        mzr <- peakArea[i, c("mzmin", "mzmax")]
        idx_rt <- which(between(rt, rtr))
        idx_mz <- which(between(mzs, c(mzr[1L] - bin_half, mzr[2L] + bin_half)))
        if (length(idx_rt) && length(idx_mz)) {
            imat <- pmat[idx_mz, idx_rt, drop = FALSE]
            if (any(imat > 0)) {
                rti <- colMax(imat, na.rm = TRUE)
                max_idx <- which.max(rti)
                idx_rt_range <- idx_rt[c(1, length(idx_rt))]
                rt_width <- diff(rt[idx_rt_range]) / diff(idx_rt_range)
                res[i, c("mz", "rt", "into", "maxo")] <- c(
                    weighted.mean(mzs[idx_mz], rowSums(imat), na.rm = TRUE),
                    rt[idx_rt[max_idx]],
                    rt_width * sum(rti, na.rm = TRUE),
                    rti[max_idx]
                )
            }
        }
    }
    res[!is.na(res[, "maxo"]), , drop = FALSE]
}

.chrom_peak_intensity_msw <- function(x, rt, peakArea,
                                      mzCenterFun = "weighted.mean",
                                      sampleIndex = integer(),
                                      cn = character(), ...) {
    res <- matrix(NA_real_, ncol = length(cn), nrow = nrow(peakArea))
    colnames(res) <- cn
    rownames(res) <- rownames(peakArea)
    res[, "sample"] <- sampleIndex
    mzc <- c("mzmin", "mzmax")
    res[, mzc] <- peakArea[, mzc]
    res[, c("rt", "rtmin", "rtmax")] <- -1
    mzs <- unlist(lapply(x, function(z) z[, "mz"]), use.names = FALSE)
    ints <- unlist(lapply(x, function(z) z[, "intensity"]), use.names = FALSE)
    for (i in 1:nrow(res)) {
        mz_area <- which(between(mzs, peakArea[i, mzc]))
        mtx <- cbind(time = -1, mz = mzs[mz_area], intensity = ints[mz_area])
        if (length(mz_area) && !all(is.na(mtx[, 3]))) {
            mzDiff <- abs(mtx[, 2] - peakArea[i, "mzmed"])
            mz_idx <- which(mzDiff == min(mzDiff))
            maxi <- mz_idx[which.max(mtx[mz_idx, 3])]
            res[i, c("mz", "maxo", "into")] <-
                c(mtx[maxi, 2:3], sum(mtx[, 3], na.rm = TRUE))
        }
    }
    res
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

#' Function to subset (eventually re-order) chromPeaks and update in addition
#' also the feature definitions with the correct indices.
#'
#' @noRd
.filter_chrom_peaks <- function(x, idx = integer()) {
    cpn <- rownames(x@chromPeaks)
    x@chromPeaks <- x@chromPeaks[idx, , drop = FALSE]
    x@chromPeakData <- x@chromPeakData[idx, ]
    if (hasFeatures(x) && (nrow(x@chromPeaks) != length(cpn)))
        x@featureDefinitions <- .update_feature_definitions(
            x@featureDefinitions, cpn, rownames(.chromPeaks(x)))
    x
}

#' Given a `matrix` with several defined regions determine for each (row)
#' the indices of the spectra in `x` that have a retention time within
#' column `"rtmin"` and `"rtmax"` and, for
#' `msLevel > 1` also a precursor m/z within `region[i, c("mzmin", "mzmax")]`
#'
#' @param x `Spectra`
#'
#' @param region `matrix` with the definition of the *region*.
#'
#' @param msLevel `integer` defining the MS level
#'
#' @return `list` with `integer` indices of the matching spectra for each
#'     region.
#'
#' @noRd
.spectra_index_list <- function(x, region, msLevel) {
    rt <- rtime(x)
    rtc <- c("rtmin", "rtmax")
    mzc <- c("mzmin", "mzmax")
    if (msLevel > 1) {
        pmz <- precursorMz(x)
        lapply(seq_len(nrow(region)), function(z) {
            which(between(rt, region[z, rtc]) &
                  between(pmz, region[z, mzc]))
        })
    } else
        lapply(seq_len(nrow(region)), function(z) {
            which(between(rt, region[z, rtc]))
        })
}

#' Same as `.spectra_index_list_closest_rt` but returns only a single index
#' per region, i.e., the index of the spectra with the retention time that
#' is closest to `region[i, "rt"]`.
#'
#' @noRd
.spectra_index_list_closest_rt <- function(x, region, msLevel) {
    rt <- rtime(x)
    mzc <- c("mzmin", "mzmax")
    rtc <- c("rtmin", "rtmax")
    if (msLevel > 1) {
        pmz <- precursorMz(x)
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]) &
                         between(pmz, region[z, mzc]))
            idx[which.min(abs(rt[idx] - region[z, "rt"]))]
        })
    } else
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]))
            idx[which.min(abs(rt[idx] - region[z, "rt"]))]
        })
}

#' Same as `.spectra_index_list_closest_rt` but returns only a single index
#' per region, i.e., the index of the spectra with the precursor m/z that
#' is closest to `region[i, "mz"]`.
#'
#' @noRd
.spectra_index_list_closest_mz <- function(x, region, msLevel) {
    rt <- rtime(x)
    rtc <- c("rtmin", "rtmax")
    mzc <- c("mzmin", "mzmax")
    pmz <- precursorMz(x)
    lapply(seq_len(nrow(region)), function(z) {
        idx <- which(between(rt, region[z, rtc]) &
                     between(pmz, region[z, mzc]))
        idx[which.min(abs(pmz[idx] - region[z, "mz"]))]
    })
}

#' Get for each region the spectrum with the largest sum of intensity
#'
#' @noRd
.spectra_index_list_largest_tic <- function(x, region, msLevel) {
    rt <- rtime(x)
    tic <- ionCount(x)
    mzc <- c("mzmin", "mzmax")
    rtc <- c("rtmin", "rtmax")
    if (msLevel > 1) {
        pmz <- precursorMz(x)
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]) &
                         between(pmz, region[z, mzc]))
            idx[which.max(tic[idx])]
        })
    } else
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]))
            idx[which.max(tic[idx])]
        })
}

#' Get for each region the spectrum with the largest intensity
#'
#' @noRd
.spectra_index_list_largest_bpi <- function(x, region, msLevel) {
    rt <- rtime(x)
    bpi <- max(intensity(x), na.rm = TRUE)
    mzc <- c("mzmin", "mzmax")
    rtc <- c("rtmin", "rtmax")
    if (msLevel > 1) {
        pmz <- precursorMz(x)
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]) &
                         between(pmz, region[z, mzc]))
            idx[which.max(bpi[idx])]
        })
    } else
        lapply(seq_len(nrow(region)), function(z) {
            idx <- which(between(rt, region[z, rtc]))
            idx[which.max(bpi[idx])]
        })
}

#' returns a `Spectra` with all spectra for the provided region.
#'
#' @param x `XcmsExperiment` with all data (for all samples).
#'
#' @param peaks `integer` with the indices of the chromatographic peaks or
#'     `integer()` for all peaks.
#'
#' @noRd
.mse_spectra_for_peaks <- function(x, method = c("all", "closest_rt",
                                                 "closest_mz", "largest_tic",
                                                 "largest_bpi"),
                                   msLevel = 2L, expandRt = 0, expandMz = 0,
                                   ppm = 0, skipFilled = FALSE,
                                   peaks = integer(), BPPARAM = bpparam()) {
    method <- match.arg(method)
    pks <- .chromPeaks(x)[, c("mz", "mzmin", "mzmax", "rt",
                              "rtmin", "rtmax", "maxo", "sample")]
    if (ppm != 0)
        expandMz <- expandMz + pks[, "mz"] * ppm / 1e6
    if (expandMz[1L] != 0) {
        pks[, "mzmin"] <- pks[, "mzmin"] - expandMz
        pks[, "mzmax"] <- pks[, "mzmax"] + expandMz
    }
    if (expandRt != 0) {
        pks[, "rtmin"] <- pks[, "rtmin"] - expandRt
        pks[, "rtmax"] <- pks[, "rtmax"] + expandRt
    }
    if (length(peaks)) keep <- peaks
    else {
        keep <- rep(TRUE, nrow(pks))
        if (skipFilled && any(.chromPeakData(x)$is_filled))
            keep <- !.chromPeakData(x)$is_filled
    }
    pks <- pks[keep, , drop = FALSE]
    f <- as.factor(as.integer(pks[, "sample"]))
    res <- bpmapply(
        split.data.frame(pks, f),
        split(spectra(x), factor(fromFile(x), levels = levels(f))),
        FUN = function(pk, sp, msLevel, method) {
            sp <- filterMsLevel(sp, msLevel)
            idx <- switch(
                method,
                all = .spectra_index_list(sp, pk, msLevel),
                closest_rt = .spectra_index_list_closest_rt(sp, pk, msLevel),
                closest_mz = .spectra_index_list_closest_mz(sp, pk, msLevel),
                largest_tic = .spectra_index_list_largest_tic(sp, pk, msLevel),
                largest_bpi = .spectra_index_list_largest_bpi(sp, pk, msLevel))
            ids <- rep(rownames(pk), lengths(idx))
            res <- sp[unlist(idx)]
            res$peak_id <- ids
            res
        },
        MoreArgs = list(msLevel = msLevel, method = method),
        BPPARAM = BPPARAM)
    Spectra:::.concatenate_spectra(res)
}

#' @param peaks `matrix` with chrom peaks.
#'
#' @param peakIdx `list` of `integer` indices defining which chromatographic
#'     peaks should be combined to a feature.
#'
#' @noRd
.manual_feature_definitions <- function(peaks, peakIdx) {
    peakIdx <- lapply(peakIdx, as.integer)
    res <- matrix(nrow = length(peakIdx), ncol = 7)
    colnames(res) <- c("mzmed", "mzmin", "mzmax", "rtmed",
                       "rtmin", "rtmax", "npeaks")
    if (!all(unlist(peakIdx) %in% seq_len(nrow(peaks))))
        stop("Some of the provided indices are out of bounds. 'peakIdx' needs",
             " to be a list of integer between 1 and ", nrow(peaks),
             call. = FALSE)
    for (i in seq_along(peakIdx)) {
        cp <- peaks[peakIdx[[i]], , drop = FALSE]
        res[i, ] <- c(median(cp[, "mz"]), min(cp[, "mz"]), max(cp[, "mz"]),
                      median(cp[, "rt"]), min(cp[, "rt"]), max(cp[, "rt"]),
                      nrow(cp))
    }
    res <- as.data.frame(res)
    res$peakidx <- peakIdx
    res
}

#' Function to help extracting the *old* MChromatograms and XChromatograms
#' from an XcmsExperiment.
#'
#' @param object `XcmsExperiment` or `XCMSnExp` object.
#'
#' @noRd
.xmse_extract_chromatograms_old <- function(object, rt, mz, aggregationFun,
                                            msLevel, isolationWindow = NULL,
                                            chunkSize, chromPeaks,
                                            return.type, BPPARAM) {
    chrs <- as(.mse_chromatogram(
        as(object, "MsExperiment"), rt = rt, mz = mz,
        aggregationFun = aggregationFun, msLevel = msLevel,
        isolationWindow = isolationWindow, chunkSize = chunkSize,
        BPPARAM = BPPARAM), return.type)
    if (return.type == "MChromatograms" || chromPeaks == "none")
        return(chrs)
    js <- seq_len(ncol(chrs))
    fd <- fData(chrs)
    rtc <- c("rtmin", "rtmax")
    mzc <- c("mzmin", "mzmax")
    if (inherits(object, "XcmsExperiment"))
        cpd <- .chromPeakData(object)
    else cpd <- as.data.frame(chromPeakData(object))
    message("Processing chromatographic peaks")
    pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                           "total (:percent) in ",
                                           ":elapsed"),
                           total = nrow(chrs) + 1L, clear = FALSE)
    for (i in seq_len(nrow(chrs))) {
        idx <- .index_chrom_peaks(
            object, rt = fd[i, rtc], mz = fd[i, mzc],
            msLevel = chrs[i, 1]@msLevel, type = chromPeaks)
        f_s <- factor(.chromPeaks(object)[idx, "sample"], levels = js)
        pkl <- split.data.frame(.chromPeaks(object)[idx, , drop = FALSE], f_s)
        cpl <- split.data.frame(cpd[idx, , drop = FALSE], f_s)
        pb$tick()
        for (j in js) {
            tmp <- chrs@.Data[i, j][[1L]]
            slot(tmp, "chromPeaks", check = FALSE) <- pkl[[j]]
            slot(tmp, "chromPeakData", check = FALSE) <- as(cpl[[j]], "DataFrame")
            chrs@.Data[i, j][[1L]] <- tmp
        }
    }
    pb$tick()
    ## Process features - that is not perfect.
    if (hasFeatures(object)) {
        message("Processing features")
        pb <- progress_bar$new(format = paste0("[:bar] :current/:",
                                               "total (:percent) in ",
                                               ":elapsed"),
                               total = nrow(chrs) + 1L, clear = FALSE)
        pks_sub <- chromPeaks(chrs)
        fts <- lapply(seq_len(nrow(chrs)), function(r) {
            fdev <- featureDefinitions(object, mz = fd[r, mzc], rt = fd[r, rtc])
            pb$tick()
            if (nrow(fdev)) {
                fdev$row <- r
                .subset_features_on_chrom_peaks(
                    fdev, .chromPeaks(object), pks_sub)
            } else data.frame()
        })
        chrs@featureDefinitions <- DataFrame(do.call(rbind, fts))
        pb$tick()
    }
    chrs@.processHistory <- object@processHistory
    chrs
}

#' @rdname XcmsExperiment
featureArea <- function(object, mzmin = min, mzmax = max, rtmin = min,
                        rtmax = max, features = character()) {
    if (!hasFeatures(object))
        stop("No correspondence results available. Please run ",
             "'groupChromPeaks' first.")
    if (!length(features))
        features <- rownames(featureDefinitions(object))
    .features_ms_region(object, mzmin = mzmin, mzmax = mzmax, rtmin = rtmin,
                        rtmax = rtmax, features = features)
}

#' @title Define MS regions for features
#'
#' @param x `XcmsExperiment` or `XCMSnExp`.
#'
#' @param mzmin, mzmax, rtmin, rtmax `function` to be applied to the values
#'     (rtmin, ...) of the chrom peaks. Defaults to `median` but would also
#'     work with `mean` etc.
#'
#' @param features `character` with the IDs of the features. Mandatory!
#'
#' @return `matrix` with columns `"mzmin"`, `"mzmax"`, `"rtmin"`, `"rtmax"`
#'     defining the range of
#'
#' @author Johannes Rainer
#'
#' @noRd
.features_ms_region <- function(x, mzmin = median, mzmax = median,
                                rtmin = median, rtmax = median,
                                features = character()) {
    features <- .i2index(features, ids = rownames(featureDefinitions(x)),
                         "features")
    pks <- .chromPeaks(x)[, c("mzmin", "mzmax", "rtmin", "rtmax")]
    res <- do.call(
        rbind, lapply(featureDefinitions(x)$peakidx[features],
                      function(i) {
                          ## maybe consider/drop gap-filled peaks?
                          c(mzmin(pks[i, "mzmin"]),
                            mzmax(pks[i, "mzmax"]),
                            rtmin(pks[i, "rtmin"]),
                            rtmax(pks[i, "rtmax"]))
                      }))
    rownames(res) <- rownames(featureDefinitions(x))[features]
    colnames(res) <- c("mzmin", "mzmax", "rtmin", "rtmax")
    res
}

#' *Reconstruct* MS2 spectra for DIA data:
#' For each MS1 chromatographic peak:
#'
#' - find all MS2 chrom peaks from the same isolation window
#' - reduce to MS2 chrom peaks with an rt similar to the one from the MS1
#'   chrom peak
#' - remove EICs with 2 or less data points
#' - create an MS2 spectrum from all MS2 chrom peaks with peak shape
#'   correlation > `minCor`.
#'
#' @param object `XcmsExperiment` with data from a **single** file.
#'
#' @author Johannes Rainer, Michael Witting
#'
#' @noRd
.mse_reconstruct_dia_ms2 <-
    function(object, expandRt = 2, diffRt = 5, minCor = 0.8, fromFile = 1L,
             column = "maxo",
             peakId = rownames(chromPeaks(object, msLevel = 1L))) {
        if (hasAdjustedRtime(object))
            object@spectra$rtime <- spectra(object)$rtime_adjusted
        message("Reconstructing MS2 spectra for ", length(peakId),
                " chrom peaks ...", appendLF = FALSE)
        pks <- .chromPeaks(object)[, c("mz", "mzmin", "mzmax", "rt", "rtmin",
                                       "rtmax", column)]
        pks[, "rtmin"] <- pks[, "rtmin"] - expandRt
        pks[, "rtmax"] <- pks[, "rtmax"] + expandRt
        ord <- order(pks[, "mz"])       # m/z need to be ordered in a Spectra
        pks <- pks[ord, ]
        ilmz <- .chromPeakData(object)$isolationWindowLowerMz[ord]
        iumz <- .chromPeakData(object)$isolationWindowUpperMz[ord]
        ## Get EICs for all chrom peaks (all MS levels)
        object <- filterRt(object, rt = range(pks[, c("rtmin", "rtmax")]))
        chrs <- .chromatograms_for_peaks(
            peaksData(object@spectra, f = factor()),
            rt = rtime(object@spectra),
            msl = msLevel(object@spectra), file_idx = fromFile,
            tmz = isolationWindowTargetMz(object@spectra), pks = pks,
            pks_msl = .chromPeakData(object)$ms_level[ord],
            pks_tmz = .chromPeakData(object)$isolationWindow[ord])
        idx <- match(peakId, rownames(pks)) # MS1 peaks to loop over
        res <- data.frame(
            peak_id = peakId, precursorMz = pks[idx, "mz"],
            rtime = pks[idx, "rt"], msLevel = 2L,
            polarity = polarity(object@spectra)[1L],
            precursorIntensity = pks[idx, column],
            fromFile = fromFile)
        ms2_peak_id <- lapply(idx, function(z) character())
        mzs <- ints <- ms2_peak_cor <-
            lapply(seq_len(nrow(res)), function(z) numeric())
        for (i in seq_along(idx)) {
            ii <- idx[i]
            imz <- .which_mz_in_range(pks[ii, "mz"], ilmz, iumz)
            irt <- .which_chrom_peak_diff_rt(pks[ii, c("rt", "rtmax")],
                                             pks, diffRt = diffRt)
            ix <- intersect(imz, irt)
            if (!length(ix))
                next
            ## Filter empty or sparse chromatograms
            ix <- ix[vapply(
                chrs@.Data[ix], function(z) sum(!is.na(z@intensity)),
                integer(1)) > 2]
            if (!length(ix))
                next
            ## Correlate
            cors <- vapply(
                chrs@.Data[ix], compareChromatograms, numeric(1),
                y = chrs[[ii]], ALIGNFUNARGS = list(method = "approx"))
            keep <- which(cors >= minCor)
            if (length(keep)) {
                ix <- ix[keep]
                res$rtime[i] <- median(pks[ix, "rt"])
                ms2_peak_id[[i]] <- rownames(pks)[ix]
                ms2_peak_cor[[i]] <- unname(cors[keep])
                mzs[[i]] <- unname(pks[ix, "mz"])
                ints[[i]] <- unname(pks[ix, column])
            }
        }
        res$mz <- mzs
        res$intensity <- ints
        res$ms2_peak_id <- ms2_peak_id
        res$ms2_peak_cor <- ms2_peak_cor
        message(" OK")
        .require_spectra()
        Spectra::Spectra(res)
    }

#' helper function to get indices of chromatographic peaks for eventual
#' subsetting with `rt`, `mz`, `msLevel` and `type`.
#'
#' @param object `XcmsExperiment` or `XCMSnExp` object with `chromPeaks`.
#'
#' @return `integer` vector with the indices of the `chromPeaks` matching the
#'     requested filtering.
#'
#' @noRd
.index_chrom_peaks <- function(object, rt = numeric(),
                               mz = numeric(), ppm = 0,
                               msLevel = integer(),
                               type = c("any", "within",
                                        "apex_within")) {
    type <- match.arg(type)
    pks <- .chromPeaks(object)
    keep <- rep(TRUE, nrow(pks))
    if (length(msLevel))
        keep <- keep &
            chromPeakData(
                object, return.type = "data.frame")$ms_level %in% msLevel
    ## Select peaks within rt range.
    if (length(rt)) {
        rt <- range(as.numeric(rt))
        if (type == "any")
            keep <- keep & pks[, "rtmin"] <= rt[2L] & pks[, "rtmax"] >= rt[1L]
        if (type == "within")
            keep <- keep & pks[, "rtmin"] >= rt[1L] & pks[, "rtmax"] <= rt[2L]
        if (type == "apex_within")
            keep <- keep & pks[, "rt"] >= rt[1L] & pks[, "rt"] <= rt[2L]
    }
    ## Select peaks within mz range, considering also ppm
    if (length(mz)) {
        mz <- range(as.numeric(mz))
        if (is.finite(mz[1L]))
            mz[1L] <- mz[1L] - mz[1L] * ppm / 1e6
        if (is.finite(mz[2L]))
            mz[2L] <- mz[2L] + mz[2L] * ppm / 1e6
        if (type == "any")
            keep <- keep & pks[, "mzmin"] <= mz[2L] & pks[, "mzmax"] >= mz[1L]
        if (type == "within")
            keep <- keep & pks[, "mzmin"] >= mz[1L] & pks[, "mzmax"] <= mz[2L]
        if (type == "apex_within")
            keep <- keep & pks[, "mz"] >= mz[1L] & pks[, "mz"] <= mz[2L]
    }
    which(keep)
}

#' Helper function to return the chromPeakData as-is (as a data.frame) from
#' a `XcmsExperiment` object.
#'
#' @noRd
.chromPeakData <- function(object, msLevel = integer()) {
    if (length(msLevel))
        object@chromPeakData[object@chromPeakData$ms_level %in% msLevel, ]
    else object@chromPeakData
}

#' Helper function to return the **full**  chromPeaks matrix without any
#' filtering from either a `XcmsExperiment` or `XCMSnExp` object
#'
#' @noRd
.chromPeaks <- function(object) {
    if (inherits(object, "XcmsExperiment"))
        object@chromPeaks
    else chromPeaks(object@msFeatureData)
}

#' function to create an empty `XcmsExperiment` object
XcmsExperiment <- function() {
    as(MsExperiment(), "XcmsExperiment")
}
