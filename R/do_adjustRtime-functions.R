## Retention time correction methods.
#' @include DataClasses.R functions-MsFeatureData.R

#' @title Align spectrum retention times across samples using peak groups
#' found in most samples
#'
#' @description The function performs retention time correction by assessing
#'     the retention time deviation across all samples using peak groups
#'     (features) containg chromatographic peaks present in most/all samples.
#'     The retention time deviation for these features in each sample is
#'     described by fitting either a polynomial (\code{smooth = "loess"}) or
#'     a linear (\code{smooth = "linear"}) model to the data points. The
#'     models are subsequently used to adjust the retention time for each
#'     spectrum in each sample. 
#'
#' @note The method ensures that returned adjusted retention times are
#'     increasingly ordered, just as the raw retention times.
#' 
#' @details The alignment bases on the presence of compounds that can be found
#'     in all/most samples of an experiment. The retention times of individual
#'     spectra are then adjusted based on the alignment of the features 
#'     corresponding to these \emph{house keeping compounds}. The paraneters
#'      \code{minFraction} and \code{extraPeaks} can be used to fine tune which
#'     features should be used for the alignment (i.e. which features 
#'     most likely correspond to the above mentioned house keeping compounds).
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
             peakGroupsMatrix = matrix(ncol = 0, nrow = 0))
{
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
        rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, nSamples,
                                     missingSample, extraPeaks)
    ## Fix for issue #175
    if (length(rt) == 0)
        stop("No peak groups found in the data for the provided settings")
    ## ## Check if we have peak groups with almost the same retention time. If yes
    ## ## select the best matching peaks among these.
    ## rtmeds <- rowMedians(rt, na.rm = TRUE)
    ## sim_rt <- which(diff(rtmeds) < 1e-6)
    ## if (length(sim_rt)) {
    ##     pk_grps <- list()
    ##     current_idxs <- NULL
    ##     last_idx <- -1
    ##     for (current_idx in sim_rt) {
    ##         if ((current_idx - last_idx) > 1) {
    ##             if (!is.null(current_idxs))
    ##                 pk_grps <- c(pk_grps, list(current_idxs))
    ##             current_idxs <- c(current_idx - 1, current_idx)
    ##         } else {
    ##             ## Just add the index.
    ##             current_idxs <- c(current_idxs, current_idx)
    ##         }
    ##         last_idx <- current_idx
    ##     }
    ##     pk_grps <- c(pk_grps, list(current_idxs))
    ##     ## Now, for each of these select one present in most samples.
    ##     sel_idx <- unlist(lapply(pk_grps, function(z) {
    ##         tmp <- rt[z, , drop = FALSE]
    ##         z[which.max(apply(tmp, MARGIN = 1, function(zz) sum(!is.na(zz))))]
    ##     }))
    ##     ## Define the other peaks that we can keep as.is
    ##     if (any(!(1:nrow(rt) %in% unique(unlist(pk_grps)))))
    ##         spec_idx <- (1:nrow(rt))[-unique(unlist(pk_grps))]
    ##     else spec_idx <- NULL
    ##     sel_idx <- sort(c(spec_idx, sel_idx))
    ##     rt <- rt[sel_idx, , drop = FALSE]
    ## }
    
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
    
    for (i in 1:nSamples) {
        pts <- na.omit(data.frame(rt = rt[, i], rtdev = rtdev[, i]))

        ## order the data.frame such that rt and rtdev are increasingly ordered.
        pk_idx <- order(pts$rt, pts$rtdev)
        pts <- pts[pk_idx, ]
        if (smooth == "loess") {
            lo <- suppressWarnings(loess(rtdev ~ rt, pts, span = span,
                                         degree = 1, family = family))
            
            rtdevsmo[[i]] <- na.flatfill(predict(lo, data.frame(rt = rtime[[i]])))
            ## Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                          quantile(abs(rtdevsmo[[i]]), 0.9,
                                   na.rm = TRUE) * 2] <- NA
            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(
                    approx(na.omit(data.frame(rtime[[i]], rtdevsmo[[i]])),
                           xout = rtime[[i]][naidx], rule = 2)$y
                )

            ## Check if there are adjusted retention times that are not ordered
            ## increasingly. If there are, search for each first unordered rt
            ## the next rt that is larger and linearly interpolate the values
            ## in between (see issue #146 for an illustration).
            while (length(decidx <- which(diff(rtime[[i]] - rtdevsmo[[i]]) < 0))) {
                warn.tweak.rt <- TRUE  ## Warn that we had to tweak the rts.
                rtadj <- rtime[[i]] - rtdevsmo[[i]]
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
                rtdevsmo[[i]][adj_idxs] <- rtime[[i]][adj_idxs] -
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

    if (warn.tweak.rt) {
        warning(call. = FALSE, "Adjusted retention times had to be ",
                "re-adjusted for some files to ensure them being in the same",
                " order than the raw retention times. A call to ",
                "'dropAdjustedRtime' might thus fail to restore retention ",
                "times of chromatographic peaks to their original values. ",
                "Eventually consider to increase the value of the 'span' ",
                "parameter.")
    }

    return(rtime)
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
   
    rt <- .getPeakGroupsRtMatrix(peaks, peakIndex, nSamples,
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

#' This function adjusts retentin times in the vector/matrix \code{x} given the
#' provided \code{numeric} vectors \code{rtraw} and \code{rtadj}.
#'
#' @details The function uses the \code{stepfun} to adjust \code{x} and adjusts
#'     it given \code{rtraw} towards \code{rtadj}. Hence it is possible to
#'     perform or to revert retention time correction in \code{x} depending
#'     on what is provided with parameters \code{rtraw} and \code{rtadj}.
#'     See examples for details.
#' 
#' @param x A numeric or matrix with retention time values that should be
#'     adjusted.
#' 
#' @noRd
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
        adjFun <- stepfun(rtraw[idx][-1] - diff(rtraw[idx]) / 2, rtadj[idx])
        ## if (!is.null(dim(x)))
        ##     return(adjFun(x[idx, ]))
        ## else 
        ##     return(adjFun(x[idx]))
    } else {
        adjFun <- stepfun(rtraw[-1] - diff(rtraw) / 2, rtadj)
    }
    adjFun(x)
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
#' @noRd
.getPeakGroupsRtMatrix <- function(peaks, peakIndex, nSamples,
                                   missingSample, extraPeaks) {
    ## For each feature:
    ## o extract the retention time of the peak with the highest intensity.
    ## o skip peak groups if they are not assigned a peak in at least a
    ##   minimum number of samples OR if have too many peaks from the same
    ##   sample assigned to it.
    seq_samp <- seq_len(nSamples)
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
        cur_fts[] <- cur_fts[order(cur_fts[, 2], decreasing = TRUE), ]
        cur_fts[match(seq_samp, cur_fts[, 3]), 1]
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
