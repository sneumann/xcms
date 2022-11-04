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
                                    ignoreHistory = FALSE, ...) {
    ## if (!length(i))
    ##     return(x)
    i <- MsCoreUtils::i2index(i, length(x))
    ## This is a special case that would make life (=performance) miserable
    if (length(i) != length(unique(i)))
        stop("Duplicated indices are not (yet) supported for ",
             "'[,XcmsExperiment'", call. = FALSE)
    ## if (keepAdjustedRtime && hasAdjustedRtime(x)) {
    ##     ## Keep only adjusted rtimes for the selected samples.
    ## }
    ## if (keepFeatures && hasFeatures(x)) {
    ## }
    if (hasChromPeaks(x)) {
        if (keepChromPeaks) {
            keep <- x@chromPeaks[, "sample"] %in% i
            x@chromPeaks <- x@chromPeaks[keep, , drop = FALSE]
            x@chromPeakData <- x@chromPeakData[keep, , drop = FALSE]
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
