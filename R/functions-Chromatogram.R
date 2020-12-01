#' @description
#'
#' *Align* chromatogram `x` to chromatogram `y`.
#' The retention times of the first `Chromatogram` (`x`) will be replaced
#' by the retention times of the second (`y`) estimated by the `approx`
#' function.
#'
#' This type of alignment should be used if the chromatograms were measured in
#' the same run (e.g. in SWATH experiments).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Michael Witting
#'
#' @noRd
.align_chromatogram_approx <- function(x, y, na.value = NA_real_, ...) {
    x_aligned <- approx(x@rtime, x@intensity, y@rtime)
    if (!is.na(na.value)) {
        x_aligned$y[is.na(x_aligned$y)] <- na.value
    }
    ## correct rtime and int in chromatogram object
    x@rtime <- x_aligned$x
    x@intensity <- x_aligned$y
    x
}

#' @description
#'
#' *Align* chromatogram `x` against chromatogram `y` matching each data point
#' in `x` to the data point in `y` with the smallest difference between their
#' retention times, but only if the difference in retention times is `<=` the
#' minimal average difference between retention times in `x` or `y` (otherwise
#' the data point will be discarded).
#'
#' @param x [Chromatogram()] object that will be aligned against `y`
#'
#' @param y [Chromatogram()] object.
#'
#' @param na.value optional parameter allowing to specify the value with which
#'     `NA`s in the aligned chromatogram `x` will be replaced.
#'
#' @param ... optional parameters to be passed along to the `.match_closest`
#'     function.
#'
#' @return aligned `x` `Chromatogram` with the same number of data points and
#'     same retention times than `y`.
#'
#' @author Johannes Rainer
#'
#' @noRd
.align_chromatogram_match_rtime <- function(x, y, na.value = NA_real_, ...) {
    idx <- closest(x@rtime, y@rtime,
                   tolerance = min(mean(diff(x@rtime)), mean(diff(y@rtime))))
    ## idx <- .match_closest(x@rtime, y@rtime, ...)
    not_na <- !is.na(idx)
    x@rtime <- rtime(y)
    new_int <- rep(na.value, length(y))
    if (any(not_na))
        new_int[idx[not_na]] <- x@intensity[not_na]
    x@intensity <- new_int
    x
}

.align_chromatogram_none <- function(x, y, na.value = NA_real_, ...) {
    idx <- match(x@rtime, y@rtime)
    not_na <- !is.na(idx)
    x@rtime <- y@rtime
    new_int <- rep(na.value, length(y))
    if (any(not_na))
        new_int[idx[not_na]] <- x@intensity[not_na]
    x@intensity <- new_int
    x
}

.align_chromatogram <- function(x, y, method = c("matchRtime", "approx", "none"),
                                na.value = NA_real_, ...) {
    method <- match.arg(method)
    switch(method,
           matchRtime = .align_chromatogram_match_rtime(
               x, y, na.value = na.value, ...),
           approx = .align_chromatogram_approx(x, y, na.value = na.value, ...),
           none = .align_chromatogram_none(x, y, na.value = na.value, ...))
}

.filter_intensity_chromatogram <- function(x, intensity = 0, ...) {
    if (is.numeric(intensity)) {
        keep <- x@intensity >= intensity[1]
    } else if (is.function(intensity)) {
        keep <- intensity(x, ...)
    } else stop("'intensity' should be either a numeric value or a function.")
    if (!is.logical(keep) | length(keep) != length(x@intensity))
        stop("The filter function seems to not return the expected result.")
    keep <- which(keep)
    x@intensity <- x@intensity[keep]
    x@rtime <- x@rtime[keep]
    x
}

#' @title Correlate chromatograms
#'
#' @description
#'
#' Correlate intensities of two chromatograms with each other. If the two
#' `Chromatogram` objects have different retention times they are first
#' *aligned* to match data points in the first to data points in the second
#' chromatogram. See [align()] for more details.
#'
#' @param x [Chromatogram()] object.
#'
#' @param y [Chromatogram()] object.
#'
#' @param use `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @param method `character(1)` passed to the `cor` function. See [cor()] for
#'     details.
#'
#' @param align `character(1)` defining the alignment method to be used. See
#'     [align()] for details. The value of this parameter is passed to the
#'     `method` parameter of `align`.
#'
#' @param ... optional parameters passed along to the `align` method.
#'
#' @return `numeric(1)` with the correlation coefficient.
#'
#' @author Michael Witting, Johannes Rainer
#'
#' @noRd
.correlate_chromatogram <- function(x, y, use = "pairwise.complete.obs",
                                    method = c("pearson", "kendall",
                                               "spearman"),
                                    align = c("matchRtime", "approx", "none"),
                                    ...) {
    align <- match.arg(align)
    if(length(x) != length(y) || !all(rtime(x) %in% rtime(y)))
        x <- .align_chromatogram(x, y, method = align, ...)
    cor(x@intensity, y@intensity, use = use, method = method)
}

#' @title Merge neighboring peaks in chromatogram
#'
#' @description
#'
#' Peak detection sometimes fails to identify a chromatographic peak correctly,
#' especially for broad peaks and if the peak shape is irregular (mostly for
#' HILIC data). In such cases several smaller peaks are reported. This function
#' tries to combine such peaks again considering their distance in retention
#' time dimension and the measured intensity between them.
#'
#' In detail, the function evaluates if the peaks are close enough, i.e. checks
#' if the difference between the `"rtmax"` and `"rtmin"` of consecutive peaks
#' is smaller than `diffRt`. If so, the average intensity of the 3 data points
#' at half way between them is calculated. If this average is larger than
#' `minProp` of the smaller `"maxo"` value of both, the peaks are merged.
#' In other words, if the (average) intensity between the two peaks is larger
#' than `minProp` times the smaller maximal peak intensity of both peaks, they
#' are joined. Note that peaks are **not** joined if all 3 data points in the
#' middle between the peaks are `NA`.
#' The joined peaks get the `"mz"`, `"rt"`, `"sn"` and `"maxo"` values from
#' the peak with the largest signal (`"maxo"`) as well as its row in the
#' *chrom peak data* `pkd`. The `"rtmin"`, `"rtmax"` are updated and
#' `"into"` is recalculated based on all the signal between `"rtmin"` and
#' `"rtmax"` of the new merged peak. The smallest and largest m/z value of all
#' data points within the provided extracted ion chromatogram is used as
#' `"mzmin"` and `"mzmax`" of the merged peak.
#'
#' Note that the `"maxo"` of the merged peak is updated (the maximum of the two
#' merged peaks is used) and used in any further merging. If for example two
#' peaks were merged, the maxo of these merged peak is used in the evaluation
#' of any additional peak that should be merged.
#'
#' @param x `Chromatogram` object with the extracted ion chromatogram containing
#'     the signal for the `pks`.
#'
#' @param pks `matrix` representing a peaks matrix. Columns `"rtmin"`,
#'     `"rtmax"` and `"maxo"` are required. It is supposed that these peaks are
#'     close enough on retention time to potentially represent signal from the
#'     same compound.
#'
#' @param pkd `DataFrame` representing the peak data (as returned by the
#'     [chromPeakData()] function.
#'
#' @param minProp `numeric(1)` representing the proportion of intensity to
#'     be required for peaks to be joined. See description for more details.
#'
#' @param diffRt `numeric(1)` representing the maximal difference between
#'     `"rtmax"` and `"rtmin"` of consecutive peaks to be considered for
#'     merging.
#'
#' @return `list` with element `"chromPeaks"`, that contains the peaks `matrix`
#'     containing newly merged peaks and original peaks if they could not be
#'     merged and `"chromPeakData"` that represents the `DataFrame` with the
#'     corresponding metadata information. The merged peaks will have a row
#'     name of `NA`.
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @noRd
#'
#' @examples
#'
#' xd <- readMSData(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
#'     mode = "onDisk")
#' chr <- chromatogram(xd, mz = c(-0.5, 0.5) + 453.2)
#' xchr <- findChromPeaks(chr, param = CentWaveParam(snthresh = 0))
#' plot(xchr)
#'
#' res <- .chrom_merge_neighboring_peaks(chr[[1]], chromPeaks(xchr[[1]]),
#'     chromPeakData(xchr[[1]]), diffRt = 2)
#' res
#' rect(res[, "rtmin"], 0, res[, "rtmax"], res[, "maxo"], border = "red")
.chrom_merge_neighboring_peaks <- function(x, pks, pkd, minProp = 0.75,
                                           diffRt = 0) {
    if (nrow(pks) < 2 || all(is.na(intensity(x))))
        return(list(chromPeaks = pks, chromPeakData = pkd))
    ## x <- clean(x, all = TRUE)
    idx <- order(pks[, "rtmin"])
    pks <- pks[idx, , drop = FALSE]
    pkd <- extractROWS(pkd, idx)
    cns <- colnames(pks)
    if (!any(cns == "mz"))
        pks <- cbind(pks, mz = NA_real_, mzmin = NA_real_, mzmax = NA_real_)
    if (is.null(rownames(pks)))
        rownames(pks) <- rownames(pkd) <- seq_len(nrow(pks))
    pks_new <- pks
    pks_new[ , ] <- NA_real_
    rownames(pks_new) <- rep(NA_character_, nrow(pks))
    pks_new[1, ] <- pks[1, ]
    rownames(pks_new)[1] <- rownames(pks)[1]
    current_peak <- 1 # point always to the current *new* (merged) peak.
    drop_cols <- !(colnames(pks_new) %in% c("mz", "mzmin", "mzmax", "rt",
                                            "rtmin", "rtmax", "into",
                                            "maxo", "sn", "sample"))
    for (i in 2:nrow(pks)) {
        if ((pks[i, "rtmin"] - pks_new[current_peak, "rtmax"]) < diffRt) {
            ## skip if second peak contained within first
            if (pks[i, "rtmin"] >= pks_new[current_peak, "rtmin"] &
                pks[i, "rtmax"] <= pks_new[current_peak, "rtmax"])
                next
            rt_mid <- (pks[i, "rtmin"] + pks_new[current_peak, "rtmax"]) / 2
            ## If rt_mid is NOT between the peaks, take the midpoint between
            ## the apexes instead.
            apexes <- range(c(pks[i, "rt"], pks[current_peak, "rt"]))
            if (rt_mid < apexes[1] || rt_mid > apexes[2])
                rt_mid <- sum(apexes) / 2
            ## Calculate the mean of the 3 data points closest to rt_mid. Skip
            ## if all of them are `NA`.
            mid_vals <- intensity(x)[order(abs(rtime(x) - rt_mid))[1:3]]
            if (!all(is.na(mid_vals)) &&
                mean(mid_vals, na.rm = TRUE) >
                min(pks_new[current_peak, "maxo"], pks[i, "maxo"]) * minProp) {
                rownames(pks_new)[current_peak] <- NA_character_
                pks_new[current_peak, drop_cols] <- NA_real_
                if (pks[i, "rtmax"] > pks_new[current_peak, "rtmax"])
                    pks_new[current_peak, "rtmax"] <- pks[i, "rtmax"]
                if (!is.na(pks[i, "mzmin"])) {
                    ## Use the mzmin and mzmax of the actual EIC - since we are
                    ## also integrating the intensities from there.
                    pks_new[current_peak, c("mzmin", "mzmax")] <-
                        mz(x, filter = FALSE)
                }
                idx_min <- which.min(
                    abs(rtime(x) - pks_new[current_peak, "rtmin"]))
                idx_max <- which.min(
                    abs(rtime(x) - pks_new[current_peak, "rtmax"]))
                peak_width <- (pks_new[current_peak, "rtmax"] -
                               pks_new[current_peak, "rtmin"]) / (idx_max - idx_min)
                ## Calculate into as done in centWave.
                pks_new[current_peak, "into"] <-
                    sum(intensity(x)[rtime(x) >= pks_new[current_peak, "rtmin"] &
                                     rtime(x) <= pks_new[current_peak, "rtmax"]],
                        na.rm = TRUE) * peak_width
                if (pks[i, "maxo"] > pks_new[current_peak, "maxo"]) {
                    pks_new[current_peak, c("mz", "rt", "maxo", "sn")] <-
                        pks[i, c("mz", "rt", "maxo", "sn")]
                    pkd[current_peak, ] <- extractROWS(pkd, i) # replace peak data with new
                }
            } else {
                current_peak <- current_peak + 1
                pks_new[current_peak, ] <- pks[i, ]
                rownames(pks_new)[current_peak] <- rownames(pks)[i]
                pkd[current_peak, ] <- extractROWS(pkd, i)
            }
        } else {
            current_peak <- current_peak + 1
            pks_new[current_peak, ] <- pks[i, ]
            rownames(pks_new)[current_peak] <- rownames(pks)[i]
            pkd[current_peak, ] <- extractROWS(pkd, i)
        }
    }
    keep <- which(!is.na(pks_new[, "rt"]))
    list(chromPeaks = pks_new[keep, cns, drop = FALSE],
         chromPeakData = extractROWS(pkd, keep))
}

.normalize_chromatogram <- function(x, method = "max") {
    ref <- do.call(method, list(x = x@intensity, na.rm = TRUE))
    x@intensity <- x@intensity / ref
    x
}
