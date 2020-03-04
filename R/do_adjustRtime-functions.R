## Retention time correction methods.
#' @include DataClasses.R functions-MsFeatureData.R

#' @title Align spectrum retention times across samples using peak groups
#' found in most samples
#'
#' @description
#'
#' The function performs retention time correction by assessing
#' the retention time deviation across all samples using peak groups
#' (features) containg chromatographic peaks present in most/all samples.
#' The retention time deviation for these features in each sample is
#' described by fitting either a polynomial (\code{smooth = "loess"}) or
#' a linear (\code{smooth = "linear"}) model to the data points. The
#' models are subsequently used to adjust the retention time for each
#' spectrum in each sample.
#'
#' @note The method ensures that returned adjusted retention times are
#'     increasingly ordered, just as the raw retention times.
#'
#' @details
#'
#' The alignment bases on the presence of compounds that can be found
#' in all/most samples of an experiment. The retention times of individual
#' spectra are then adjusted based on the alignment of the features
#' corresponding to these \emph{house keeping compounds}. The paraneters
#' \code{minFraction} and \code{extraPeaks} can be used to fine tune which
#' features should be used for the alignment (i.e. which features
#' most likely correspond to the above mentioned house keeping compounds).
#'
#' Parameter \code{subset} allows to define a subset of samples within the
#' experiment that should be aligned. All samples not being part of the subset
#' will be aligned based on the adjustment of the closest sample within the
#' subset. This allows to e.g. exclude blank samples from the alignment process
#' with their retention times being still adjusted based on the alignment
#' results of the \emph{real} samples.
#'
#' @inheritParams adjustRtime-peakGroups
#'
#' @param peaks a \code{matrix} or \code{data.frame} with the identified
#'     chromatographic peaks in the samples.
#'
#' @param peakIndex a \code{list} of indices that provides the grouping
#'     information of the chromatographic peaks (across and within samples).
#'
#' @param rtime a \code{list} of \code{numeric} vectors with the retention
#'     times per file/sample.
#'
#' @param peakGroupsMatrix optional \code{matrix} of (raw) retention times for
#'     peak groups on which the alignment should be performed. Each column
#'     represents a sample, each row a feature/peak group. If not provided,
#'     this matrix will be determined depending on parameters
#'     \code{minFraction} and \code{extraPeaks}. If provided,
#'     \code{minFraction} and \code{extraPeaks} will be ignored.
#'
#' @return A \code{list} with \code{numeric} vectors with the adjusted
#'     retention times grouped by sample.
#'
#' @family core retention time correction algorithms
#'
#' @author Colin Smith, Johannes Rainer
#'
#' @references
#' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
#' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
#' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
#' \emph{Anal. Chem.} 2006, 78:779-787.
do_adjustRtime_peakGroups <-
    function(peaks, peakIndex, rtime, minFraction = 0.9, extraPeaks = 1,
             smooth = c("loess", "linear"), span = 0.2,
             family = c("gaussian", "symmetric"),
             peakGroupsMatrix = matrix(ncol = 0, nrow = 0),
             subset = integer(), subsetAdjust = c("average", "previous"))
{
    subsetAdjust <- match.arg(subsetAdjust)
    ## Check input.
    if (missing(peaks) | missing(peakIndex) | missing(rtime))
        stop("Arguments 'peaks', 'peakIndex' and 'rtime' are required!")
    smooth <- match.arg(smooth)
    family <- match.arg(family)
    ## minFraction
    if (any(minFraction > 1) | any(minFraction < 0))
        stop("'minFraction' has to be between 0 and 1!")
    ## Check peaks:
    OK <- xcms:::.validChromPeaksMatrix(peaks)
    if (is.character(OK))
        stop(OK)
    ## Check peakIndex:
    if (any(!(unique(unlist(peakIndex)) %in% seq_len(nrow(peaks)))))
        stop("Some indices listed in 'peakIndex' are outside of ",
             "1:nrow(peaks)!")
    ## Check rtime:
    if (!is.list(rtime))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (!all(unlist(lapply(rtime, is.numeric), use.names = FALSE)))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (length(rtime) != max(peaks[, "sample"]))
        stop("The length of 'rtime' does not match with the total number of ",
             "samples according to the 'peaks' matrix!")
    total_samples <- length(rtime)
    if (length(subset)) {
        if (!is.numeric(subset))
            stop("If provided, 'subset' is expected to be an integer")
        if (!all(subset %in% seq_len(total_samples)))
            stop("One or more indices in 'subset' are out of range.")
        if (length(subset) < 2)
            stop("Length of 'subset' too small: minimum required samples for ",
                 "alignment is 2.")
    } else subset <- seq_len(total_samples)
    ## Translate minFraction to number of allowed missing samples.
    nSamples <- length(subset)
    missingSample <- nSamples - (nSamples * minFraction)
    ## Remove peaks not present in "subset" from the peakIndex
    peaks_in_subset <- which(peaks[, "sample"] %in% subset)
    peakIndex <- lapply(peakIndex, function(z) z[z %in% peaks_in_subset])
    ## Check if we've got a valid peakGroupsMatrix
    ## o Same number of samples.
    ## o range of rt values is within the rtime.
    if (nrow(peakGroupsMatrix)) {
        if (ncol(peakGroupsMatrix) != nSamples)
            stop("'peakGroupsMatrix' has to have the same number of columns ",
                 "as there are samples!")
        pg_range <- range(peakGroupsMatrix, na.rm = TRUE)
        rt_range <- range(rtime)
        if (!(pg_range[1] >= rt_range[1] & pg_range[2] <= rt_range[2]))
            stop("The retention times in 'peakGroupsMatrix' have to be within",
                 " the retention time range of the experiment!")
        rt <- peakGroupsMatrix
    } else
        rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, subset,
                                     missingSample, extraPeaks)
    ## Fix for issue #175
    if (length(rt) == 0)
        stop("No peak groups found in the data for the provided settings")
    if (ncol(rt) != length(subset))
        stop("Length of 'subset' and number of columns of the peak group ",
             "matrix do not match.")
    message("Performing retention time correction using ", nrow(rt),
            " peak groups.")

    ## Calculate the deviation of each peak group in each sample from its
    ## median
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)

    if (smooth == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            smooth <- "linear"
            warning("Too few peak groups for 'loess', reverting to linear",
                    " method")
        } else if (mingroups * span < 4) {
            span <- 4 / mingroups
            warning("Span too small for 'loess' and the available number of ",
                    "peak groups, resetting to ", round(span, 2))
        }
    }

    rtdevsmo <- vector("list", nSamples)

    ## Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE
    warn.tweak.rt <- FALSE

    rtime_adj <- rtime
    ## Adjust samples in subset.
    for (i in seq_along(subset)) {
        i_all <- subset[i]              # Index of sample in whole dataset.
        pts <- na.omit(data.frame(rt = rt[, i], rtdev = rtdev[, i]))

        ## order the data.frame such that rt and rtdev are increasingly ordered.
        pk_idx <- order(pts$rt, pts$rtdev)
        pts <- pts[pk_idx, ]
        if (smooth == "loess") {
            lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span,
                                         degree = 1, family = family))

            rtdevsmo[[i]] <- xcms:::na.flatfill(
                                        predict(lo, data.frame(rt = rtime[[i_all]])))
            ## Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                          quantile(abs(rtdevsmo[[i]]), 0.9,
                                   na.rm = TRUE) * 2] <- NA
            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(
                    approx(na.omit(data.frame(rtime[[i_all]], rtdevsmo[[i]])),
                           xout = rtime[[i_all]][naidx], rule = 2)$y
                )

            ## Check if there are adjusted retention times that are not ordered
            ## increasingly. If there are, search for each first unordered rt
            ## the next rt that is larger and linearly interpolate the values
            ## in between (see issue #146 for an illustration).
            while (length(decidx <- which(diff(rtime[[i_all]] - rtdevsmo[[i]]) < 0))) {
                warn.tweak.rt <- TRUE  ## Warn that we had to tweak the rts.
                rtadj <- rtime[[i_all]] - rtdevsmo[[i]]
                rtadj_start <- rtadj[decidx[1]] ## start interpolating from here
                ## Define the
                next_larger <- which(rtadj > rtadj[decidx[1]])
                if (length(next_larger) == 0) {
                    ## Fix if there is no larger adjusted rt up to the end.
                    next_larger <- length(rtadj) + 1
                    rtadj_end <- rtadj_start
                } else {
                    next_larger <- min(next_larger)
                    rtadj_end <- rtadj[next_larger]
                }
                ## linearly interpolate the values in between.
                adj_idxs <- (decidx[1] + 1):(next_larger - 1)
                incr <- (rtadj_end - rtadj_start) / length(adj_idxs)
                rtdevsmo[[i]][adj_idxs] <- rtime[[i_all]][adj_idxs] -
                    (rtadj_start + (1:length(adj_idxs)) * incr)
            }

            rtdevsmorange <- range(rtdevsmo[[i]])
            if (any(rtdevsmorange / rtdevrange > 2))
                warn.overcorrect <- TRUE
        } else {
            if (nrow(pts) < 2) {
                stop("Not enough peak groups even for linear smoothing ",
                     "available!")
            }
            ## Use lm instead?
            fit <- lsfit(pts$rt, pts$rtdev)
            rtdevsmo[[i]] <- rtime[[i_all]] * fit$coef[2] + fit$coef[1]
            ptsrange <- range(pts$rt)
            minidx <- rtime[[i_all]] < ptsrange[1]
            maxidx <- rtime[[i_all]] > ptsrange[2]
            rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
            rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
        }
        ## Finally applying the correction
        rtime_adj[[i_all]] <- rtime[[i_all]] - rtdevsmo[[i]]
    }
    ## Adjust the remaining samples.
    rtime_adj <- adjustRtimeSubset(rtime, rtime_adj, subset = subset,
                                   method = subsetAdjust)
    if (warn.overcorrect) {
        warning("Fitted retention time deviation curves exceed points by more",
                " than 2x. This is dangerous and the algorithm is probably ",
                "overcorrecting your data. Consider increasing the span ",
                "parameter or switching to the linear smoothing method.")
    }

    if (warn.tweak.rt) {
        warning(call. = FALSE, "Adjusted retention times had to be ",
                "re-adjusted for some files to ensure them being in the same",
                " order than the raw retention times. A call to ",
                "'dropAdjustedRtime' might thus fail to restore retention ",
                "times of chromatographic peaks to their original values. ",
                "Eventually consider to increase the value of the 'span' ",
                "parameter.")
    }
    rtime_adj
}
## That's the original code that fails to fix unsorted adjusted retention times
## (see issue #146).
do_adjustRtime_peakGroups_orig <- function(peaks, peakIndex, rtime,
                                           minFraction = 0.9, extraPeaks = 1,
                                           smooth = c("loess", "linear"), span = 0.2,
                                           family = c("gaussian", "symmetric")) {
    ## Check input.
    if (missing(peaks) | missing(peakIndex) | missing(rtime))
        stop("Arguments 'peaks', 'peakIndex' and 'rtime' are required!")
    smooth <- match.arg(smooth)
    family <- match.arg(family)
    ## minFraction
    if (any(minFraction > 1) | any(minFraction < 0))
        stop("'minFraction' has to be between 0 and 1!")

    ## Check peaks:
    OK <- .validChromPeaksMatrix(peaks)
    if (is.character(OK))
        stop(OK)
    ## Check peakIndex:
    if (any(!(unique(unlist(peakIndex)) %in% seq_len(nrow(peaks)))))
        stop("Some indices listed in 'peakIndex' are outside of ",
             "1:nrow(peaks)!")
    ## Check rtime: in line with the total number of samples we've got in
    ## peaks?
    if (!is.list(rtime))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (!all(unlist(lapply(rtime, is.numeric), use.names = FALSE)))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (length(rtime) != max(peaks[, "sample"]))
        stop("The length of 'rtime' does not match with the total number of ",
             "samples according to the 'peaks' matrix!")

    nSamples <- length(rtime)
    ## Translate minFraction to number of allowed missing samples.
    missingSample <- nSamples - (nSamples * minFraction)

    rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, seq_len(nSamples),
                                 missingSample, extraPeaks)

    message("Performing retention time correction using ", nrow(rt),
            " peak groups.")

    ## Calculate the deviation of each peak group in each sample from its
    ## median
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)

    if (smooth == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            smooth <- "linear"
            warning("Too few peak groups for 'loess', reverting to linear",
                    " method")
        } else if (mingroups * span < 4) {
            span <- 4 / mingroups
            warning("Span too small for 'loess' and the available number of ",
                    "peak groups, resetting to ", round(span, 2))
        }
    }

    rtdevsmo <- vector("list", nSamples)

    ## Code for checking to see if retention time correction is overcorrecting
    rtdevrange <- range(rtdev, na.rm = TRUE)
    warn.overcorrect <- FALSE

    for (i in 1:nSamples) {
        pts <- na.omit(data.frame(rt = rt[, i], rtdev = rtdev[, i]))

        if (smooth == "loess") {
            lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span,
                                         degree = 1, family = family))

            rtdevsmo[[i]] <- na.flatfill(predict(lo, data.frame(rt = rtime[[i]])))
            ## Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                          quantile(abs(rtdevsmo[[i]]), 0.9) * 2] <- NA
            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(
                    approx(na.omit(data.frame(rtime[[i]], rtdevsmo[[i]])),
                           xout = rtime[[i]][naidx], rule = 2)$y
                )
            ## That's to ensure that the adjusted retention times are in the
            ## same order than the raw retention times - I guess...
            ## Check if adjustment changes the order of the adjusted retention
            ## times. If yes, the difference between consecutive retention times
            ## will be negative.
            ## What does this code do:
            ## o move the last adjusted retention time to the left by its
            ##   difference to the next one.
            while (length(decidx <- which(diff(rtime[[i]] - rtdevsmo[[i]]) < 0))) {
                d <- diff(rtime[[i]] - rtdevsmo[[i]])[tail(decidx, 1)]
                rtdevsmo[[i]][tail(decidx, 1)] <- rtdevsmo[[i]][tail(decidx, 1)] - d
                if (abs(d) <= 1e-06)
                    break
            }

            rtdevsmorange <- range(rtdevsmo[[i]])
            if (any(rtdevsmorange / rtdevrange > 2))
                warn.overcorrect <- TRUE
        } else {
            if (nrow(pts) < 2) {
                stop("Not enough peak groups even for linear smoothing ",
                     "available!")
            }
            ## Use lm instead?
            fit <- lsfit(pts$rt, pts$rtdev)
            rtdevsmo[[i]] <- rtime[[i]] * fit$coef[2] + fit$coef[1]
            ptsrange <- range(pts$rt)
            minidx <- rtime[[i]] < ptsrange[1]
            maxidx <- rtime[[i]] > ptsrange[2]
            rtdevsmo[[i]][minidx] <- rtdevsmo[[i]][head(which(!minidx), n = 1)]
            rtdevsmo[[i]][maxidx] <- rtdevsmo[[i]][tail(which(!maxidx), n = 1)]
        }
        ## Finally applying the correction
        rtime[[i]] <- rtime[[i]] - rtdevsmo[[i]]
    }

    if (warn.overcorrect) {
        warning("Fitted retention time deviation curves exceed points by more",
                " than 2x. This is dangerous and the algorithm is probably ",
                "overcorrecting your data. Consider increasing the span ",
                "parameter or switching to the linear smoothing method.")
    }

    return(rtime)
}

#' This function adjusts retentin times in the vector/matrix `x` given the
#' provided `numeric` vectors `rtraw` and `rtadj`.
#'
#' @note
#'
#' Values in `x` that are outside of the range of `rtraw` are linearly shifted
#' by the difference of the first/last adjusted value.
#'
#' @details
#'
#' The function uses `stepfun` to adjust `x` given adjustment from `rtraw`
#' to `rtadj`. It is possible to perform or to revert retention time
#' correction in `x` depending on whether raw or adjusted retention times are
#' provided with `rtraw` or `rtadj`, respectively.
#' See examples for details.
#'
#' @param x A `numeric` or `matrix` with retention time values that should be
#'     adjusted.
#'
#' @param rtraw `numeric` with raw retention times.
#'
#' @param rtadj `numeric` with adjusted retention times.
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @md
#'
#' @examples
#'
#' ## Perform retention time correction:
#' ## feats is supposed to be the peaks matrix FOR A SINGLE SAMPLE, rtr and
#' ## rtc the raw and adjusted retention times of the spectras from the same
#' ## samples:
#' ## adjFts <- feats
#' ## adjFts[, c("rt", "rtmin", "rtmax")] <- .applyRtAdjustment(feats[, c("rt", "rtmin", "rtmax")], rtr, rtc)
#'
#' ## To revert the adjustment: just switch the order of rtr and rtc
.applyRtAdjustment <- function(x, rtraw, rtadj) {
    ## re-order everything if rtraw is not sorted; issue #146
    if (is.unsorted(rtraw)) {
        idx <- order(rtraw)
        rtraw <- rtraw[idx]
        rtadj <- rtadj[idx]
    }
    adjFun <- stepfun(rtraw[-1] - diff(rtraw) / 2, rtadj)
    res <- adjFun(x)
    ## Fix margins.
    idx_low <- which(x < rtraw[1])
    if (length(idx_low)) {
        first_adj <- idx_low[length(idx_low)] + 1
        res[idx_low] <- x[idx_low] + res[first_adj] - x[first_adj]
    }
    idx_high <- which(x > rtraw[length(rtraw)])
    if (length(idx_high)) {
        last_adj <- idx_high[1] - 1
        res[idx_high] <- x[idx_high] + res[last_adj] - x[last_adj]
    }
    if (is.null(dim(res)))
        names(res) <- names(x)
    res
}

#' Helper function to apply retention time adjustment to already identified
#' peaks in the peaks matrix of an XCMSnExp (or peaks matrix of an
#' xcmsSet).
#'
#' @noRd
.applyRtAdjToChromPeaks <- function(x, rtraw, rtadj) {
    if (!is.list(rtraw) | !is.list(rtadj))
        stop("'rtraw' and 'rtadj' are supposed to be lists!")
    if (length(rtraw) != length(rtadj))
        stop("'rtraw' and 'rtadj' have to have the same length!")
    ## Going to adjust the columns rt, rtmin and rtmax in x.
    ## Using a for loop here.
    for (i in 1:length(rtraw)) {
        whichSample <- which(x[, "sample"] == i)
        if (length(whichSample)) {
            x[whichSample, c("rt", "rtmin", "rtmax")] <-
                .applyRtAdjustment(x[whichSample, c("rt", "rtmin", "rtmax")],
                                   rtraw = rtraw[[i]], rtadj = rtadj[[i]])
        }
    }
    x
}

#' Simple helper function to create a matrix with retention times for well
#' aligned peak groups, each row containing the rt of a peak group,
#' columns being samples.
#'
#' @details This function is called internally by the
#'     do_adjustRtime_peakGroups function and the retcor.peakgroups method.
#'
#' @noRd
.getPeakGroupsRtMatrix <- function(peaks, peakIndex, sampleIndex,
                                   missingSample, extraPeaks) {
    ## For each feature:
    ## o extract the retention time of the peak with the highest intensity.
    ## o skip peak groups if they are not assigned a peak in at least a
    ##   minimum number of samples OR if have too many peaks from the same
    ##   sample assigned to it.
    nSamples <- length(sampleIndex)
    rt <- lapply(peakIndex, function(z) {
        cur_fts <- peaks[z, c("rt", "into", "sample"), drop = FALSE]
        ## Return NULL if we've got less samples that required or is the total
        ## number of peaks is larger than a certain threshold.
        ## Note that the original implementation is not completely correct!
        ## nsamp > nsamp + extraPeaks might be correct.
        nsamp <- length(unique(cur_fts[, "sample"]))
        if (nsamp < (nSamples - missingSample) |
            nrow(cur_fts) > (nsamp + extraPeaks))
            return(NULL)
        cur_fts[] <- cur_fts[order(cur_fts[, 2], decreasing = TRUE), , drop = FALSE]
        cur_fts[match(sampleIndex, cur_fts[, 3]), 1]
    })
    rt <- do.call(rbind, rt)
    ## Order them by median retention time. NOTE: this is different from the
    ## original code, in which the peak groups are ordered by the median
    ## retention time that is calculated over ALL peaks within the peak
    ## group, not only to one peak selected for each sample (for multi
    ## peak per sample assignments).
    ## Fix for issue #175
    if (is(rt, "matrix")) {
        rt <- rt[order(rowMedians(rt, na.rm = TRUE)), , drop = FALSE]
    }
    rt
}

#' For a given index `x` return one from `idx` that is the *closest*, can be
#' either simply the `"next"`, or the `"closet"` (smallest difference to `x`).
#' For `"next"`: if there is no *next* index, it takes the previous.
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' .get_closest_index(3, c(2, 4, 6, 8))
.get_closest_index <- function(x, idx, method = c("next", "previous",
                                                  "closest")) {
    method <- match.arg(method)
    switch(method,
           `next` = {
               nxt <- idx > x
               if (any(nxt))
                   idx[nxt][1]
               else idx[!nxt][sum(!nxt)]
           },
           `previous` = {
               prv <- idx < x
               if (any(prv))
                   idx[prv][sum(prv)]
               else idx[!prv][1]
           },
           closest = {
               dst <- abs(idx - x)
               idx[which.min(dst)]
           })
}

#' *align* two vectors with each other. If they have a different length they
#' will be trimmed to have the same length starting from the end with the
#' smaller average difference (between 5 datapoints).
#'
#' @param x `list` of two `numeric` vectors with potentially different length.
#'
#' @return `list` of two `numeric` vectors with equal length.
#'
#' @noRd
#'
#' @author Johannes Rainer
#'
#' @examples
#'
#' x <- list(a = 1:10, b = 3:10)
#' .match_trim_vectors(x)
#'
#' x <- list(a = 1:10, b = 2:15)
#' .match_trim_vectors(x)
#'
#' x <- list(a = 1:20, b = 1:5)
#' .match_trim_vectors(x)
.match_trim_vector_index <- function(x, n = 5) {
    lens <- lengths(x)
    min_len <- min(lens)
    if (length(unique(lens)) == 1)
        replicate(n = length(lens), seq_len(min_len))
    hd <- vapply(x, function(z) mean(head(z, n = n)), numeric(1))
    tl <- vapply(x, function(z) mean(tail(z, n = n)), numeric(1))
    if (diff(range(hd)) <= diff(range(tl)))
        replicate(n = length(x), 1:min_len, FALSE)
    else
        lapply(lens, function(z) (z - min_len + 1):z)
}

.match_trim_vectors <- function(x, n = 5, idxs) {
    if (missing(idxs))
        idxs <- .match_trim_vector_index(x, n = n)
    mapply(x, idxs, FUN = function(z, idx) z[idx], SIMPLIFY = FALSE)
}

#' @title Adjust retention times based on alignment result from subset
#'
#' @description
#'
#' This function adjusts retention times based on the alignment results on a
#' subset of samples from an experiment. Specifically, the samples **not** part
#' of `subset` are adjusted based on the adjusted retention times of the
#' *closest* `subset` samples:
#'
#' How the retention times will be adjusted depends on `method`:
#' - `"previous"`: adjusted retention times will match the adjusted retention
#'   times of the closest previous subset-sample.
#' - `"average"`: adjusted retention times for a non-subset sample is calculated
#'   based on the average retention times of the closest previous and following
#'   subset-sample. The average is calculated based on an weighted average
#'   with weights representing the distance of the non-subset sample to the
#'   closest subset samples.
#'
#' @param rtraw `list` of raw retention times, one element/vector per sample.
#'
#' @param rtadj `list` of adjusted retention times, one element/vector per
#'     sample. Has to have the same length than `rtraw`.
#'
#' @param subset `integer` with the indices of the `subset` on which the
#'     alignment has been performed.
#'
#' @param method `character` specifying the method with which the non-subset
#'     samples are adjusted: either `"previous"` or `"average"`. See details.
#'
#' @return `list` of adjusted retention times.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @md
adjustRtimeSubset <- function(rtraw, rtadj, subset,
                              method = c("average", "previous")) {
    method <- match.arg(method)
    if (length(rtraw) != length(rtadj))
        stop("Lengths of 'rtraw' and 'rtadj' have to match.")
    if (missing(subset))
        subset <- seq_along(rtraw)
    if (!all(subset %in% seq_along(rtraw)))
        stop("'subset' is out of bounds.")
    ## if (length(subset) == length(rtraw)) {
    ##     cat("return rtadj\n")
    ##     return(rtadj)
    ## }
    no_subset <- seq_len(length(rtraw))[-subset]
    for (i in no_subset) {
        message("Aligning sample number ", i, " against subset ... ",
                appendLF = FALSE)
        if (method == "previous") {
            i_adj <- xcms:::.get_closest_index(i, subset, method = "previous")
            rtadj[[i]] <- .applyRtAdjustment(rtraw[[i]], rtraw[[i_adj]],
                                                 rtadj[[i_adj]])
        }
        if (method == "average") {
            i_ref <- c(xcms:::.get_closest_index(i, subset, method = "previous"),
                       xcms:::.get_closest_index(i, subset, method = "next"))
            trim_idx <- .match_trim_vector_index(rtraw[i_ref])
            rt_raw_ref <- do.call(
                cbind, .match_trim_vectors(rtraw[i_ref], idxs = trim_idx))
            rt_adj_ref <- do.call(
                cbind, .match_trim_vectors(rtadj[i_ref], idxs = trim_idx))
            wghts <- 1 / abs(i_ref - i) # weights depending on distance to i
            rt_raw_ref <- apply(rt_raw_ref, 1, weighted.mean, w = wghts)
            rt_adj_ref <- apply(rt_adj_ref, 1, weighted.mean, w = wghts)
            rtadj[[i]] <- .applyRtAdjustment(rtraw[[i]], rt_raw_ref,
                                             rt_adj_ref)
        }
        message("OK")
    }
    rtadj
}
