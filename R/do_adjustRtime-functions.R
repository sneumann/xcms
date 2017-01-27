## Retention time correction methods.
#' @include DataClasses.R functions-MsFeatureData.R

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
##' @inheritParams adjustRtime-featureGroups
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
##' @return A \code{list} with \code{numeric} vectors with the adjusted retention
##' times grouped by sample.
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
   
    rt <- .getFeatureGroupsRtMatrix(features, featureIndex, nSamples,
                                    missingSample, extraFeatures)

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

    return(rtime)
}

## Creating a "do" function for this is inefficient, since we need the profile
## matrix for all samples.
do_adjustRtime_obiwarp <- function(features, featureIndex, rtime,
                                   profStep = 1, centerSample = NULL,
                                   response = 1,
                                   distFun = "cor_opt", gapInit = NULL,
                                   gapExtend = NULL, factorDiag = 2,
                                   factorGap = 1, localAlignment = FALSE,
                                   initPenalty = 0) {
    ## Check input.
    if (missing(features) | missing(featureIndex) | missing(rtime))
        stop("Arguments 'features', 'featureIndex' and 'rtime' are required!")    
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
    ## centerSample
    if (!is.null(centerSample)) {
        if (length(centerSample) > 1 | !(centerSample %in% 1:nSamples))
            stop("'centerSample' has to be a single integer between 1 and ",
                 nSamples, "!")
    }
    ## Note: localAlignment has to be converted to a number.

    ## set default parameter.
    if (is.null(gapInit)) {
        if (distFunc == "cor") gapInit = 0.3
        if (distFunc == "cor_opt") gapInit = 0.3
        if (distFunc == "cov") gapInit = 0.0
        if (distFunc == "euc") gapInit = 0.9
        if (distFunc == "prd") gapInit = 0.0
    }
    if (is.null(gapExtend)) {
        if (distFunc == "cor") gapExtend = 2.4
        if (distFunc == "cor_opt") gapExtend = 2.4
        if (distFunc == "cov") gapExtend = 11.7
        if (distFunc == "euc") gapExtend = 1.8
        if (distFunc == "prd") gapExtend = 7.8
    }

    ## peakmat: features
    ## rtcor: rtime
    ## samples
    
    ## peakmat <- peaks(object)
    ## samples <- sampnames(object)
    ## classlabel <- as.vector(unclass(sampclass(object)))
    ## N <- length(samples)
    ## corpeaks <- peakmat
    ## plottype <- match.arg(plottype)

    ## Initialize result lists.
    rtimecor <- vector("list", nSamples)
    rtdevsmo <- vector("list", nSamples)

    ## Determine the center fun
    if (is.null(centerSample))
        centerSample <- which.max(table(features[, "sample"]))[1]

    ## Now for XCMSnExp (or OnDiskMSnExp) objects:
    ## o Have to ensure that the scan range of all objects fits with the center
    ##   sample - let's to that within the lapply though - just like in the
    ##   original function.
    ## o split it by file
    ## o Extract the center sample, remove it from the list
    ## o do an lapply over them:
    
    message("Sample number ", centerSample, " used as center sample.")
    
    cat("center sample: ", samples[center], "\nProcessing: ")
    idx <- which(seq(1,N) != center)
    obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0)

	## added t automatically find the correct scan range from the xcmsSet object
	if(length(obj1@scantime) != length(object@rt$raw[[center]])){
		##figure out the scan time range
		scantime.start	<-object@rt$raw[[center]][1]
		scantime.end	<-object@rt$raw[[center]][length(object@rt$raw[[center]])]

		scanrange.start	<-which.min(abs(obj1@scantime - scantime.start))
		scanrange.end	<-which.min(abs(obj1@scantime - scantime.end))
		scanrange<-c(scanrange.start, scanrange.end)
		obj1 <- xcmsRaw(object@filepaths[center], profmethod="bin", profstep=0, scanrange=scanrange)
	} else{
		scanrange<-NULL
	}

    for (si in 1:length(idx)) {
        s <- idx[si]
        cat(samples[s], " ")

        ##
        ## Might be better to just get the profile matrix from the center object
        ## outside of the for loop and then modifying a internal variable within
        ## the loop - avoids creation of two profile matrices in each iteration.
        profStepPad(obj1) <- profStep ## (re-)generate profile matrix, since it might have been modified during previous iteration
		if(is.null(scanrange)){
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0)
		} else{
			obj2 <- xcmsRaw(object@filepaths[s], profmethod="bin", profstep=0, scanrange=scanrange)
		}
        profStepPad(obj2) <- profStep ## generate profile matrix

        mzmin <-  min(obj1@mzrange[1], obj2@mzrange[1])
        mzmax <-  max(obj1@mzrange[2], obj2@mzrange[2])

        mz <- seq(mzmin,mzmax, by=profStep)
        mz <- as.double(mz)
        mzval <- length(mz)

        scantime1 <- obj1@scantime
        scantime2 <- obj2@scantime

        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        rtup1 <- c(1:length(scantime1))
        rtup2 <- c(1:length(scantime2))

        ## ???
        mst1 <- which(diff(scantime1)>5*mstdiff)[1]
        if(!is.na(mst1)) {
            rtup1 <- which(rtup1<=mst1)
            cat("Found gaps: cut scantime-vector at ", scantime1[mst1],"seconds", "\n")
        }

        mst2 <- which(diff(scantime2)>5*mstdiff)[1]
        if(!is.na(mst2)) {
            rtup2 <- which(rtup2<=mst2)
            cat("Found gaps: cut scantime-vector at ", scantime2[mst2],"seconds", "\n")
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]

        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        if(rtmaxdiff>(5*mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            rtup1 <- which(scantime1<=rtmax)
            rtup2 <- which(scantime2<=rtmax)
        }

        scantime1 <- scantime1[rtup1]
        scantime2 <- scantime2[rtup2]
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)

        if(length(obj1@scantime)>valscantime1) {
            obj1@env$profile <- obj1@env$profile[,-c((valscantime1+1):length(obj1@scantime))]
        }
        if(length(obj2@scantime)>valscantime2) {
            obj2@env$profile <- obj2@env$profile[,-c((valscantime2+1):length(obj2@scantime))]
        }

        if(mzmin < obj1@mzrange[1]) {
            seqlen <- length(seq(mzmin, obj1@mzrange[1], profStep))-1
            x <- matrix(0, seqlen,dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(x, obj1@env$profile)
        }
        if(mzmax > obj1@mzrange[2]){
            seqlen <- length(seq(obj1@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj1@env$profile)[2])
            obj1@env$profile <- rbind(obj1@env$profile, x)
        }
        if(mzmin < obj2@mzrange[1]){
            seqlen <- length(seq(mzmin, obj2@mzrange[1], profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(x, obj2@env$profile)
        }
        if(mzmax > obj2@mzrange[2]){
            seqlen <- length(seq(obj2@mzrange[2], mzmax, profStep))-1
            x <- matrix(0, seqlen, dim(obj2@env$profile)[2])
            obj2@env$profile <- rbind(obj2@env$profile, x)
        }

        intensity1 <- obj1@env$profile
        intensity2 <- obj2@env$profile

        if ((mzval * valscantime1 != length(intensity1)) ||  (mzval * valscantime2 != length(intensity2)))
            stop("Dimensions of profile matrices do not match !\n")

        ## Would it be possible to supply non-binned data too???
        rtimecor[[s]] <-.Call("R_set_from_xcms",
                              valscantime1,scantime1,mzval,mz,intensity1,
                              valscantime2,scantime2,mzval,mz,intensity2,
                              response, distFunc,
                              gapInit, gapExtend,
                              factorDiag, factorGap,
                              localAlignment, initPenalty)

        if(length(obj2@scantime) > valscantime2) {
            object@rt$corrected[[s]] <- c(rtimecor[[s]],
                                          obj2@scantime[(max(rtup2)+1):length(obj2@scantime)])
        } else {
            object@rt$corrected[[s]] <- rtimecor[[s]]
        }

        rtdevsmo[[s]] <- round(rtime[[s]]-object@rt$corrected[[s]],2)

        rm(obj2)
        gc()

        ## updateProgressInfo
        object@progressInfo$retcor.obiwarp <-  si / length(idx)
        xcms:::progressInfoUpdate(object)

    }

    cat("\n")
    rtdevsmo[[center]] <- round(rtime[[center]] - object@rt$corrected[[center]], 2)


    for (i in 1:N) {
        cfun <- stepfun(rtime[[i]][-1] - diff(rtime[[i]])/2, rtime[[i]] - rtdevsmo[[i]])
        rtime[[i]] <- rtime[[i]] - rtdevsmo[[i]]

        sidx <- which(corpeaks[,"sample"] == i)
        corpeaks[sidx, c("rt", "rtmin", "rtmax")] <- cfun(corpeaks[sidx, c("rt", "rtmin", "rtmax")])
    }

    peaks(object) <- corpeaks
    groups(object) <- matrix(nrow = 0, ncol = 0)
    groupidx(object) <- list()
    invisible(object)
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
.applyRtAdjToFeatures <- function(x, rtraw, rtadj) {
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

##' Simple helper function to create a matrix with retention times for well
##' aligned feature groups, each row containing the rt of a feature group,
##' columns being samples.
##' @details This function is called internally by the
##' do_adjustRtime_featureGroups function and the retcor.peakgroups method.
##' @noRd
.getFeatureGroupsRtMatrix <- function(features, featureIndex, nSamples,
                                      missingSample, extraFeatures) {
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
    return(rt[order(rowMedians(rt, na.rm = TRUE)), , drop = FALSE])
}
