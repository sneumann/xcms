## Retention time correction methods.
#' @include functions-MsFeatureData.R

## Needs:
## o features
## o feature groups
## o rtime

##' @title Align spectrum retention times across samples using feature groups
##' found in most samples
##' 
##' @description The function performs retention time correction by assessing
##' the retention time deviation across all samples using feature groups containg
##' features present in most/all samples. The retention time deviation
##' for these feature groups in each sample is described by fitting either a
##' polynomial (\code{smooth = "loess"}) or a linear (\code{smooth = "linear"})
##' model to the data points. The models are subsequently used to adjust the
##' retention time for each spectrum in each sample. 
##'
##' @details The alignment bases on the presence of compounds that can be found
##' in all/most samples of an experiment. The retention times of individual
##' spectra are then adjusted based on the alignment of the feature groups
##' corresponding to these \emph{house keeping compounds}. The paraneters
##' \code{minFraction} and \code{extraFeatures} can be used to fine tune which
##' feature groups should be used for the alignment (i.e. which feature groups
##' most likely correspond to the above mentioned house keeping compounds).
##' 
##' @param features a \code{matrix} or \code{data.frame} with the identified
##' features in the samples.
##'
##' @param featureIndex a \code{list} of indices that provides the grouping
##' information of the features (across and within samples).
##' 
##' @param rtime A \code{list} of \code{numeric} vectors with the retention times
##' per file/sample.
##' 
##' @param minFraction numeric(1) between 0 and 1 defining the minimum required
##' fraction of samples in which features for the feature group were identified.
##' Feature groups passing this criteria will aligned across samples and retention
##' times of individual spectra will be adjusted based on this alignment. For
##' \code{minFraction = 1} the feature group has to contain features in all
##' samples of the experiment.
##' 
##' @param extraFeatures numeric(1) defining the maximal number of additional
##' features for all samples to be assigned to a feature group for retention time
##' correction. For a data set with 6 samples, \code{extraFeatures = 1} uses all
##' feature groups with a total feature count \code{<= 6 + 1}. The total feature
##' count is the total number of features being assigned to a feature group and
##' considers also multiple features within a sample being assigned to the group.
##'
##' @param smooth character defining the function to be used, to interpolate
##' corrected retention times for all feature groups. Either \code{"loess"} or
##' \code{"linear"}.
##'
##' @param span numeric(1) defining the degree of smoothing (if
##' \code{smooth = "loess"}). This parameter is passed to the internal call
##' to \code{\link{loess}}.
##'
##' @param family character defining the method to be used for loess smoothing.
##' Allowed values are \code{"gaussian"} and \code{"symmetric"}.See
##' \code{\link{loess}} for more information.
##' 
##' @return A \code{numeric} vector with the adjusted retention times.
##'
##' @family core retention time correction algorithms
##'
##' @author Colin Smith, Johannes Rainer
##'
##' @references
##' Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
##' Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
##' Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
##' \emph{Anal. Chem.} 2006, 78:779-787.
do_adjustRtime_featureGroups <- function(features, featureIndex, rtime,
                                         minFraction = 0.9, extraFeatures = 1,
                                         smooth = c("loess", "linear"), span = 0.2,
                                         family = c("gaussian", "symmetric")) {
    ## Check input.
    if (missing(features) | missing(featureIndex) | missing(rtime))
        stop("Arguments 'features', 'featureIndex' and 'rtime' are required!")
    smooth <- match.arg(smooth)
    family <- match.arg(family)
    ## minFraction
    if (any(minFraction > 1) | any(minFraction < 0))
        stop("'minFraction' has to be between 0 and 1!")
    
    ## Check features:
    OK <- .validFeatureMatrix(features)
    if (is.character(OK))
        stop(OK)
    ## Check featureIndex:
    if (any(!(unique(unlist(featureIndex)) %in% seq_len(nrow(features)))))
        stop("Some indices listed in 'featureIndex' are outside of ",
             "1:nrow(features)!")
    ## Check rtime: in line with the total number of samples we've got in
    ## features?
    if (!is.list(rtime))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (!all(unlist(lapply(rtime, is.numeric), use.names = FALSE)))
        stop("'rtime' should be a list of numeric vectors with the retention ",
             "times of the spectra per sample!")
    if (length(rtime) != max(features[, "sample"]))
        stop("The length of 'rtime' does not match with the total number of ",
             "samples according to the 'features' matrix!")
    
    nSamples <- length(rtime)
    ## Translate minFraction to number of allowed missing samples.
    missingSample <- nSamples - (nSamples * minFraction)
   
    ## For each feature group:
    ## o extract the retention time of the feature with the highest intensity.
    ## o skip feature groups if they are not assigned a feature in at least a
    ##   minimum number of samples OR if have too many features from the same
    ##   sample assigned to it.
    seq_samp <- seq_len(nSamples)
    rt <- lapply(featureIndex, function(z) {
        cur_fts <- features[z, c("rt", "into", "sample"), drop = FALSE]
        ## Return NULL if we've got less samples that required or is the total
        ## number of features is larger than a certain threshold.
        ## Note that the original implementation is not completely correct!
        ## nsamp > nsamp + extraFeatures might be correct.
        nsamp <- length(unique(cur_fts[, "sample"]))
        if (nsamp < (nSamples - missingSample) |
            nrow(cur_fts) > (nsamp + extraFeatures))
            return(NULL)
        cur_fts[] <- cur_fts[order(cur_fts[, 2], decreasing = TRUE), ]
        cur_fts[match(seq_samp, cur_fts[, 3]), 1]
    })
    rt <- do.call(rbind, rt)
    ## Order them by median retention time. NOTE: this is different from the
    ## original code, in which the feature groups are ordered by the median
    ## retention time that is calculated over ALL features within the feature
    ## group, not only to one feature selected for each sample (for multi
    ## feature per sample assignments).
    rt <- rt[order(rowMedians(rt, na.rm = TRUE)), , drop = FALSE]
    
    message("Performing retention time correction using ", nrow(rt),
            " feature groups.")

    ## Calculate the deviation of each feature group in each sample from its
    ## median
    rtdev <- rt - apply(rt, 1, median, na.rm = TRUE)

    if (smooth == "loess") {
        mingroups <- min(colSums(!is.na(rt)))
        if (mingroups < 4) {
            smooth <- "linear"
            warning("Too few feature groups for 'loess', reverting to linear",
                    " method")
        } else if (mingroups * span < 4) {
            span <- 4 / mingroups
            warning("Span too small for 'loess' and the available number of ",
                    "feature groups, resetting to ", round(span, 2))
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
            
            rtdevsmo[[i]] <- xcms:::na.flatfill(predict(lo,
                                                        data.frame(rt = rtime[[i]])))
            ## Remove singularities from the loess function
            rtdevsmo[[i]][abs(rtdevsmo[[i]]) >
                          quantile(abs(rtdevsmo[[i]]), 0.9) * 2] <- NA
            if (length(naidx <- which(is.na(rtdevsmo[[i]]))))
                rtdevsmo[[i]][naidx] <- suppressWarnings(
                    approx(na.omit(data.frame(rtime[[i]], rtdevsmo[[i]])),
                           xout = rtime[[i]][naidx], rule = 2)$y
                )
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
                stop("Not enough feature groups even for linear smoothing ",
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

    ## Oh my god! We're really adjusting the retention times for the identified
    ## features here!!! This should be done in the xcmsSet or the XCMSnExp
    ## object!!!
    ## for (i in 1:n) {

    ##     cfun <- stepfun(rtcor[[i]][-1] - diff(rtcor[[i]])/2, rtcor[[i]] - rtdevsmo[[i]])
    ##     rtcor[[i]] <- rtcor[[i]] - rtdevsmo[[i]]

    ##     sidx <- which(corpeaks[,"sample"] == i)
    ##     corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    ## }

    ## object@rt$corrected <- rtcor
    ## peaks(object) <- corpeaks
    ## groups(object) <- matrix(nrow = 0, ncol = 0)
    ## groupidx(object) <- list()
    ## invisible(object)
    return(rtime)
}

do_adjustRtime_obiwarp <- function() {
}

##' This function adjusts retentin times in the vector/matrix \code{x} given the
##' provided \code{numeric} vectors \code{rtraw} and \code{rtadj}.
##'
##' @details The function uses the \code{stepfun} to adjust \code{x} and adjusts
##' it given \code{rtraw} towards \code{rtadj}. Hence it is possible to perform
##' or to revert retention time correction in \code{x} depending on what is
##' provided with parameters \code{rtraw} and \code{rtadj}. See examples for
##' details.
##' 
##' @param x A numeric or matrix with retention time values that should be
##' adjusted.
##' @noRd
##' @examples
##'
##' ## Perform retention time correction:
##' ## feats is supposed to be the features matrix FOR A SINGLE SAMPLE, rtr and
##' ## rtc the raw and adjusted retention times of the spectras from the same
##' ## samples:
##' ## adjFts <- feats
##' ## adjFts[, c("rt", "rtmin", "rtmax")] <- .applyRtAdjustment(feats[, c("rt", "rtmin", "rtmax")], rtr, rtc)
##'
##' ## To revert the adjustment: just switch the order of rtr and rtc
.applyRtAdjustment <- function(x, rtraw, rtadj) {
    adjFun <- stepfun(rtraw[-1] - diff(rtraw) / 2, rtadj)
    return(adjFun(x))
}

##' Helper function to apply retention time adjustment to already identified
##' features in the features matrix of an XCMSnExp (or peaks matrix of an
##' xcmsSet).
##' 
##' @noRd
.applyRtAdjustmentToFeatures <- function(x, rtraw, rtadj) {
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
    return(x)
}
