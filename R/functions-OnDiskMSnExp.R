## Functions for MSnbase's OnDiskMSnExp objects
#' @include do_findChromPeaks-functions.R DataClasses.R


#' @param x an OnDiskMSnExp representing the whole experiment.
#'
#' @param method The (chromatographic) peak detection method to be used. Can be
#'     "centWave" etc.
#'
#' @param param A class extending Param containing all parameters for the
#'     peak detection method.
#'
#' @return a list of length 2, \code{peaks} containing a matrix with the
#'     identified peaks and \code{date} the time stamp when the peak detection
#'     was started.
#'
#' @noRd
findChromPeaks_OnDiskMSnExp <- function(object, method = "centWave",
                                        param) {
    require("xcms", quietly = TRUE, character.only = TRUE)
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    ## Since we're calling this function already with bplapply ensure that
    ## the spectra call is not firing its own parallel processing!
    findChromPeaks_Spectrum_list(x = spectra(object, BPPARAM = SerialParam()),
                                 method = method, param = param,
                                 rt = rtime(object))
}


#' Run the peak detection on a list of Spectrum1 objects from the same
#' file
#'
#' @param x A list of Spectrum1 objects of a sample.
#'
#' @param method The peak detection method to be used. Can be "centWave" etc.
#'
#' @param param A class extending Param containing all parameters for the
#'     peak detection method.
#'
#' @param rt Numeric with the retention times for the spectra. If not provided
#'     it is extracted from the spectra.
#'
#' @return a list of length 2, \code{peaks} containing a matrix with the
#'     identified peaks and \code{date} the time stamp when the peak detection
#'     was started.
#'
#' @author Johannes Rainer
#'
#' @noRd
findChromPeaks_Spectrum_list <- function(x, method = "centWave", param, rt) {
    method <- match.arg(method, c("centWave", "massifquant", "matchedFilter",
                                  "MSW", "centWaveWithPredIsoROIs"))
    method <- paste0("do_findChromPeaks_", method)
    if (method == "do_findChromPeaks_MSW")
        method <- "do_findPeaks_MSW"
    if (method == "do_findChromPeaks_matchedFilter") {
        ## Issue #325: empty spectra is not supported
        x <- lapply(x, function(z) {
            if (!length(z@mz)) {
                z@mz <- 0.0
                z@intensity <- 0.0
            }
            z
        })
    }
    if (missing(param))
        stop("'param' has to be specified!")
    if (missing(rt))
        rt <- unlist(lapply(x, rtime), use.names = FALSE)
    if (is.unsorted(rt))
        stop("Spectra are not ordered by retention time!")
    mzs <- lapply(x, mz)
    vals_per_spect <- lengths(mzs, FALSE)
    if (any(vals_per_spect == 0))
        warning("Found empty spectra. Please run 'filterEmptySpectra' first.",
                call. = FALSE)
    procDat <- date()
    res <- do.call(
        method, args = c(list(mz = unlist(mzs, use.names = FALSE),
                              int = unlist(lapply(x, intensity),
                                           use.names = FALSE),
                              valsPerSpect = vals_per_spect,
                              scantime = rt), as(param, "list")))
    ## Ensure that we call the garbage collector to eventually clean unused stuff
    rm(mzs)
    rm(x)
    rm(rt)
    gc()
    list(peaks = res, date = procDat)
}

## That's a special case since we don't expect to have rt available for this.
findPeaks_MSW_OnDiskMSnExp <- function(object, method = "MSW",
                                            param) {
    require("xcms", quietly = TRUE, character.only = TRUE)
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    findPeaks_MSW_Spectrum_list(x = spectra(object, BPPARAM = SerialParam()),
                                method = method, param = param)
}
findPeaks_MSW_Spectrum_list <- function(x, method = "MSW", param) {
    method <- match.arg(method, c("MSW"))
    method <- paste0("do_findPeaks_", method)
    if (missing(param))
        stop("'param' has to be specified!")
    mzs <- lapply(x, mz)
    procDat <- date()
    list(peaks = do.call(
             method, args = c(list(mz = unlist(mzs, use.names = FALSE),
                                   int = unlist(lapply(x, intensity),
                                                use.names = FALSE)),
                              as(param, "list"))), date = procDat)
}

#' Fast way to split an `XCMSnExp` object by file:
#' - subsets the object to MS level specified with `msLevel`.
#' - subsets feature data to the smallest possible set of columns (if
#'   `selectFeatureData = TRUE`.
#' - returns by default an `OnDiskMSnExp`, unless `to_class = "XCMSnExp"`, in
#'   which case also potentially present chromatographic peaks are preserved.
#'
#' @param keep_sample_idx if column "sample" should be kept as it is.
#'
#' @note
#'
#' This function needs a considerable amount of memory if
#' `to_class = "XCMSnExp"` because, for efficiency reasons, it first splits
#' the `chromPeaks` and `chromPeakData` per file.
#'
#' @noRd
.split_by_file <- function(x, msLevel. = unique(msLevel(x)),
                           subsetFeatureData = TRUE,
                           to_class = "OnDiskMSnExp",
                           keep_sample_idx = FALSE) {
    if (is(x, "XCMSnExp") && hasAdjustedRtime(x))
        x@featureData$retentionTime <- adjustedRtime(x)
    if (subsetFeatureData) {
        fcs <- intersect(c(MSnbase:::.MSnExpReqFvarLabels, "centroided",
                           "polarity", "seqNum"), colnames(.fdata(x)))
        x <- selectFeatureData(x, fcol = fcs)
    }
    procd <- x@processingData
    expd <- new(
        "MIAPE",
        instrumentManufacturer = x@experimentData@instrumentManufacturer[1],
        instrumentModel = x@experimentData@instrumentModel[1],
        ionSource = x@experimentData@ionSource[1],
        analyser = x@experimentData@analyser[1],
        detectorType = x@experimentData@detectorType[1])
    create_object <- function(x, i, to_class) {
        a <- new(to_class)
        slot(procd, "files", check = FALSE) <- x@processingData@files[i]
        slot(a, "processingData", check = FALSE) <- procd
        slot(a, "featureData", check = FALSE) <- extractROWS(
            x@featureData, which(x@featureData$msLevel %in% msLevel. &
                                 x@featureData$fileIdx == i))
        if (!nrow(a@featureData))
            stop("No MS level ", msLevel., " spectra present.", call. = FALSE)
        a@featureData$fileIdx <- 1L
        slot(a, "experimentData", check = FALSE) <- expd
        slot(a, "spectraProcessingQueue", check = FALSE) <-
            x@spectraProcessingQueue
        slot(a, "phenoData", check = FALSE) <- x@phenoData[i, , drop = FALSE]
        a
    }
    if (to_class == "XCMSnExp" && is(x, "XCMSnExp") && hasChromPeaks(x)) {
        if (any(colnames(.chrom_peak_data(x@msFeatureData)) == "ms_level")) {
            pk_idx <- which(
                .chrom_peak_data(x@msFeatureData)$ms_level %in% msLevel.)
        } else pk_idx <- seq_len(nrow(chromPeaks(x@msFeatureData)))
        fct <- as.factor(
            as.integer(chromPeaks(x@msFeatureData)[pk_idx, "sample"]))
        pksl <- split.data.frame(
            chromPeaks(x@msFeatureData)[pk_idx, , drop = FALSE], fct)
        pkdl <- split.data.frame(
            extractROWS(.chrom_peak_data(x@msFeatureData), pk_idx), fct)
        res <- vector("list", length(fileNames(x)))
        for (i in seq_along(res)) {
            a <- create_object(x, i, to_class)
            newFd <- new("MsFeatureData")
            pks <- pksl[[as.character(i)]]
            if (!is.null(pks) && nrow(pks)) {
                if (!keep_sample_idx)
                    pks[, "sample"] <- 1
                chromPeaks(newFd) <- pks
                chromPeakData(newFd) <- pkdl[[as.character(i)]]
            } else {
                chromPeaks(newFd) <- chromPeaks(x@msFeatureData)[0, ]
                chromPeakData(newFd) <- .chrom_peak_data(x@msFeatureData)[0, ]
            }
            lockEnvironment(newFd, bindings = TRUE)
            slot(a, "msFeatureData", check = FALSE) <- newFd
            res[[i]] <- a
        }
        res
    } else {
        lapply(seq_along(fileNames(x)), function(z) {
            create_object(x, z, to_class)
        })
    }
}

#' Same as `.split_by_file` but *faster* because it splits the chrom peaks
#' matrix too - and requires thus more memory.
#'
#' @author Johannes Rainer
#'
#' @noRd
.split_by_file2 <- function(x, msLevel. = unique(msLevel(x)),
                            subsetFeatureData = FALSE,
                            to_class = "OnDiskMSnExp",
                            keep_sample_idx = FALSE) {
    if (is(x, "XCMSnExp") && hasAdjustedRtime(x))
        x@featureData$retentionTime <- adjustedRtime(x)
    if (subsetFeatureData) {
        fcs <- intersect(c(MSnbase:::.MSnExpReqFvarLabels, "centroided",
                           "polarity", "seqNum"), colnames(.fdata(x)))
        x <- selectFeatureData(x, fcol = fcs)
    }
    fdl <- split.data.frame(x@featureData, as.factor(fromFile(x)))
    procd <- x@processingData
    expd <- new(
        "MIAPE",
        instrumentManufacturer = x@experimentData@instrumentManufacturer[1],
        instrumentModel = x@experimentData@instrumentModel[1],
        ionSource = x@experimentData@ionSource[1],
        analyser = x@experimentData@analyser[1],
        detectorType = x@experimentData@detectorType[1])
    create_object <- function(i, fd, x, to_class) {
        a <- new(to_class)
        slot(procd, "files", check = FALSE) <- x@processingData@files[i]
        slot(a, "processingData", check = FALSE) <- procd
        slot(a, "featureData", check = FALSE) <-
            extractROWS(fd, which(fd$msLevel %in% msLevel.))
        if (!nrow(a@featureData))
            warning("No MS level ", msLevel., " spectra present in file ",
                    basename(x@processingData@files[i]), call. = FALSE)
        else a@featureData$fileIdx <- 1L
        slot(a, "experimentData", check = FALSE) <- expd
        slot(a, "spectraProcessingQueue", check = FALSE) <-
            x@spectraProcessingQueue
        slot(a, "phenoData", check = FALSE) <- x@phenoData[i, , drop = FALSE]
        a
    }
    if (to_class == "XCMSnExp" && is(x, "XCMSnExp") && hasChromPeaks(x)) {
        if (any(colnames(.chrom_peak_data(x@msFeatureData)) == "ms_level")) {
            pk_idx <- which(
                .chrom_peak_data(x@msFeatureData)$ms_level %in% msLevel.)
        } else pk_idx <- seq_len(nrow(chromPeaks(x@msFeatureData)))
        fct <- as.factor(
            as.integer(chromPeaks(x@msFeatureData)[pk_idx, "sample"]))
        pksl <- split.data.frame(
            chromPeaks(x@msFeatureData)[pk_idx, , drop = FALSE], fct)
        pkdl <- split.data.frame(
            extractROWS(.chrom_peak_data(x@msFeatureData), pk_idx), fct)
        res <- vector("list", length(fileNames(x)))
        for (i in seq_along(res)) {
            a <- create_object(i, fdl[[i]], x, to_class)
            newFd <- new("MsFeatureData")
            pks <- pksl[[as.character(i)]]
            if (!is.null(pks) && nrow(pks)) {
                if (!keep_sample_idx)
                    pks[, "sample"] <- 1
                chromPeaks(newFd) <- pks
                chromPeakData(newFd) <- pkdl[[as.character(i)]]
            } else {
                chromPeaks(newFd) <- chromPeaks(x@msFeatureData)[0, ]
                chromPeakData(newFd) <- .chrom_peak_data(x@msFeatureData)[0, ]
            }
            lockEnvironment(newFd, bindings = TRUE)
            slot(a, "msFeatureData", check = FALSE) <- newFd
            res[[i]] <- a
        }
        res
    } else
        mapply(seq_along(fileNames(x)), fdl, FUN = create_object,
               MoreArgs = list(x = x, to_class = to_class))
}


############################################################
#' @description Fill some settings and data from an OnDiskMSnExp or pSet into an
#' xcmsSet:
#' o fileNames
#' o phenoData
#' o rt
#' o mslevel
#' o scanrange: this should enable to subset the raw data again by scan index
#' (which should be equivalent to acquisitionNum/scanIndex)
#'
#' @param pset The pSet from which data should be extracted
#'
#' @noRd
.pSet2xcmsSet <- function(pset) {
    object <- new("xcmsSet")
    filepaths(object) <- fileNames(pset)
    phenoData(object) <- pData(pset)
    ## rt
    rt <- split(unname(rtime(pset)), f = as.factor(fromFile(pset)))
    object@rt <- list(raw = rt, corrected = rt)
    ## mslevel
    mslevel(object) <- unique(msLevel(pset))
    ## scanrange
    ## scanrange(object) <- range(acquisitionNum(pset))
    ## Using scanIndex instead of acquisitionNum as the scan index should
    ## represent the index of the spectrum within the file, while acquisitionNum
    ## might be different e.g. if the mzML was subsetted.
    scanrange(object) <- range(scanIndex(pset))
    object
}

#' @description Processes the result list returned by an lapply/bplapply to
#'     findChromPeaks_Spectrum_list or findChromPeaks_OnDiskMSnExp and returns a
#'     list with two elements: \code{$peaks} the peaks matrix of identified
#'     peaks and \code{$procHist} a list of ProcessHistory objects (empty if
#'     \code{getProcHist = FALSE}).
#'
#' @param x See description above.
#'
#' @param getProcHist Wheter ProcessHistory objects should be returned too.
#'
#' @param fnames The file names.
#'
#' @return See description above.
#'
#' @author Johannes Rainer
#'
#' @noRd
.processResultList <- function(x, getProcHist = TRUE, fnames) {
    nSamps <- length(x)
    pks <- vector("list", nSamps)
    phList <- vector("list", nSamps)
    for (i in 1:nSamps) {
        n_pks <- nrow(x[[i]]$peaks)
        if (is.null(n_pks))
            n_pks <- 0
        if (n_pks == 0) {
            pks[[i]] <- NULL
            warning("No peaks found in sample number ", i, ".")
        } else {
            pks[[i]] <- cbind(x[[i]]$peaks, sample = rep.int(i, n_pks))
        }
        if (getProcHist)
            phList[[i]] <- ProcessHistory(
                info. = paste0("Chromatographic peak detection in '",
                               basename(fnames[i]), "': ", n_pks,
                               " peaks identified."),
                date. = x[[i]]$date,
                type. = .PROCSTEP.PEAK.DETECTION,
                fileIndex. = i
            )
    }
    list(peaks = pks, procHist = phList)
}


#' Calculate adjusted retention times by aligning each sample against a center
#' sample.
#'
#' @note Adjustment should be performed only on spectra from the same MS level!
#'     It's up to the calling function to ensure that.
#'
#' @param object An \code{OnDiskMSnExp}.
#'
#' @param param An \code{ObiwarpParam}.
#'
#' @param msLevel \code{integer} defining the MS level on which the adjustment
#'     should be performed.
#'
#' @return The function returns a \code{list} of adjusted retention times
#'     grouped by file.
#'
#' @noRd
.obiwarp <- function(object, param) {
    if (missing(object))
        stop("'object' is mandatory!")
    if (missing(param))
        param <- ObiwarpParam()
    subs <- subset(param)
    if (!length(subs))
        subs <- seq_along(fileNames(object))
    total_samples <- length(fileNames(object))
    nSamples <- length(subs)
    if (nSamples <= 1)
        stop("Can not perform a retention time correction on less than two",
             " files.")

    ## centerSample
    if (length(centerSample(param))) {
        if (!(centerSample(param) %in% 1:nSamples))
            stop("'centerSample' has to be a single integer between 1 and ",
                 nSamples, "!")
    } else
        centerSample(param) <- floor(median(1:nSamples))
    message("Sample number ", centerSample(param), " used as center sample.")

    rtraw <- split(rtime(object), as.factor(fromFile(object)))
    object <- filterFile(object, file = subs)
    ## Get the profile matrix of the center sample:
    ## Using the (hidden) parameter returnBreaks to return also the breaks of
    ## the bins of the profile matrix. I can use them to align the matrices
    ## later.
    ## NOTE: it might be event better to just re-use the breaks from the center
    ## sample for the profile matrix generation of all following samples.
    suppressMessages(
        profCtr <- profMat(object, method = "bin",
                           step = binSize(param),
                           fileIndex = centerSample(param),
                           returnBreaks = TRUE)[[1]]
    )
    ## Now split the object by file
    objL <- .split_by_file2(object, msLevel. = 1)
    objL <- objL[-centerSample(param)]
    centerObject <- filterFile(object, file = centerSample(param))
    ## Now we can bplapply here!
    res <- bplapply(objL, function(z, cntr, cntrPr, parms) {
        message("Aligning ", basename(fileNames(z)), " against ",
                basename(fileNames(cntr)), " ... ", appendLF = FALSE)
        ## Get the profile matrix for the current file.
        suppressMessages(
            curP <- profMat(z, method = "bin", step = binSize(parms),
                            returnBreaks = TRUE)[[1]]
        )
        rtadj <- .obiwarp_bare(scantime1 = unname(rtime(cntr)),
                               scantime2 = unname(rtime(z)),
                               ref_pr = cntrPr, other_pr = curP,
                               param = parms)
        message("OK")
        rtadj
    }, cntr = centerObject, cntrPr = profCtr, parms = param)
    ## Create result
    adjRt <- vector("list", total_samples)
    adjRt[subs[centerSample(param)]] <- list(unname(rtime(centerObject)))
    adjRt[subs[-centerSample(param)]] <- res
    adjustRtimeSubset(rtraw, adjRt, subset = subs, method = subsetAdjust(param))
}

.obiwarp_bare <- function(scantime1, scantime2, ref_pr, other_pr, param) {
    ## ---------------------------------------
    ## 1)Check the scan times of both objects:
    scantime1_orig <- scantime1         # reference
    scantime2_orig <- scantime2         # sample
    scantime1_diff <- diff(scantime1)
    scantime2_diff <- diff(scantime2)
    ## median difference between spectras' scan time.
    mstdiff <- median(c(scantime1_diff, scantime2_diff), na.rm = TRUE)

    mst1 <- which(scantime1_diff > param@rtimeDifferenceThreshold * mstdiff)[1]
    if (!is.na(mst1)) {
        warning("Found gaps in scan times of the center sample: cut ",
                "scantime-vector at ", scantime1[mst1]," seconds.")
        scantime1 <- scantime1[seq_len(max(2, (mst1 - 1)))]
    }
    mst2 <- which(scantime2_diff > param@rtimeDifferenceThreshold * mstdiff)[1]
    if (!is.na(mst2)) {
        warning("Found gaps in scan time. Cutting scantime-vector at ",
                scantime2[mst2]," seconds.")
        scantime2 <- scantime2[seq_len(max(2, (mst2 - 1)))]
    }
    ## Drift of measured scan times - expected to be largest at the end.
    rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                            scantime2[length(scantime2)])))
    ## If the drift is larger than the threshold, cut the matrix up to the
    ## max allowed difference.
    if (rtmaxdiff > (5 * mstdiff)) {
        rtmax <- min(scantime1[length(scantime1)],
                     scantime2[length(scantime2)])
        scantime1 <- scantime1[scantime1 <= rtmax]
        scantime2 <- scantime2[scantime2 <= rtmax]
    }
    valscantime1 <- length(scantime1)
    valscantime2 <- length(scantime2)
    ## Ensure we have the same number of scans.
    if (valscantime1 != valscantime2) {
        min_number <- min(valscantime1, valscantime2)
        diffs <- abs(range(scantime1) - range(scantime2))
        ## Cut at the start or at the end, depending on where we have the
        ## larger difference
        if (diffs[2] > diffs[1]) {
            scantime1 <- scantime1[1:min_number]
            scantime2 <- scantime2[1:min_number]
        } else {
            scantime1 <- rev(rev(scantime1)[1:min_number])
            scantime2 <- rev(rev(scantime2)[1:min_number])
        }
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)
    }
    ## Finally, restrict the profile matrix to the restricted data
    if (ncol(ref_pr$profMat) != valscantime1) {
        ## Find out whether we were cutting at the start or end.
        start_idx <- which(scantime1[1] == scantime1_orig)
        end_idx <- which(scantime1[length(scantime1)] == scantime1_orig)
        ref_pr$profMat <- ref_pr$profMat[, start_idx:end_idx]
    }
    if(ncol(other_pr$profMat) != valscantime2) {
        start_idx <- which(scantime2[1] == scantime2_orig)
        end_idx <- which(scantime2[length(scantime2)] == scantime2_orig)
        other_pr$profMat <- other_pr$profMat[, start_idx:end_idx]
    }
    ## ---------------------------------
    ## 2) Now match the breaks/mz range.
    ##    The -1 below is because the breaks define the upper and lower
    ##    boundary. Have to do it that way to be in line with the orignal
    ##    code... would be better to use the breaks as is.
    mzr1 <- c(ref_pr$breaks[1], ref_pr$breaks[length(ref_pr$breaks) - 1])
    mzr2 <- c(other_pr$breaks[1], other_pr$breaks[length(other_pr$breaks) - 1])
    mzmin <- min(c(mzr1[1], mzr2[1]))
    mzmax <- max(c(mzr1[2], mzr2[2]))
    mzs <- seq(mzmin, mzmax, by = binSize(param))
    ## Eventually add empty rows at the beginning
    if (mzmin < mzr1[1]) {
        tmp <- matrix(0, (length(seq(mzmin, mzr1[1], binSize(param))) - 1),
                      ncol = ncol(ref_pr$profMat))
        ref_pr$profMat <- rbind(tmp, ref_pr$profMat)
    }
    ## Eventually add empty rows at the end
    if (mzmax > mzr1[2]) {
        tmp <- matrix(0, (length(seq(mzr1[2], mzmax, binSize(param))) - 1),
                      ncol = ncol(ref_pr$profMat))
        ref_pr$profMat <- rbind(ref_pr$profMat, tmp)
    }
    ## Eventually add empty rows at the beginning
    if (mzmin < mzr2[1]) {
        tmp <- matrix(0, (length(seq(mzmin, mzr2[1], binSize(param))) - 1),
                      ncol = ncol(other_pr$profMat))
        other_pr$profMat <- rbind(tmp, other_pr$profMat)
    }
    ## Eventually add empty rows at the end
    if (mzmax > mzr2[2]) {
        tmp <- matrix(0, (length(seq(mzr2[2], mzmax, binSize(param))) - 1),
                      ncol = ncol(other_pr$profMat))
        other_pr$profMat <- rbind(other_pr$profMat, tmp)
    }
    ## A final check of the data.
    mzvals <- length(mzs)
    cntrVals <- length(ref_pr$profMat)
    curVals <- length(other_pr$profMat)
    if ((mzvals * valscantime1) != cntrVals | (mzvals * valscantime2) != curVals)
        stop("Dimensions of profile matrices do not match!")
    ## Done with preparatory stuff - now I can perform the alignment.
    rtadj <- .Call("R_set_from_xcms", valscantime1, scantime1, mzvals, mzs,
                   ref_pr$profMat, valscantime2, scantime2, mzvals, mzs,
                   other_pr$profMat, response(param), distFun(param),
                   gapInit(param), gapExtend(param), factorDiag(param),
                   factorGap(param), as.numeric(localAlignment(param)),
                   initPenalty(param))
    if (length(scantime2_orig) != valscantime2) {
        nrt <- length(scantime2_orig)
        adj_starts_at <- which(scantime2_orig == scantime2[1])
        adj_ends_at <- which(scantime2_orig == scantime2[length(scantime2)])
        if (adj_ends_at < nrt)
            rtadj <- c(
                rtadj, rtadj[length(rtadj)] +
                       cumsum(diff(scantime2_orig[adj_ends_at:nrt])))
        if (adj_starts_at > 1)
            rtadj <- c(
                rtadj[1] +
                rev(cumsum(diff(scantime2_orig[adj_starts_at:1]))), rtadj)
    }
    unname(rtadj)
}


.concatenate_OnDiskMSnExp <- function(...) {
    x <- list(...)
    if (length(x) == 0)
        return(NULL)
    if (length(x) == 1)
        return(x[[1]])
    ## Check that all are XCMSnExp objects.
    if (!all(unlist(lapply(x, function(z) is(z, "OnDiskMSnExp")))))
        stop("All passed objects should be 'OnDiskMSnExp' objects")
    ## Check processingQueue
    procQ <- lapply(x, function(z) z@spectraProcessingQueue)
    new_procQ <- procQ[[1]]
    is_ok <- unlist(lapply(procQ, function(z)
        !is.character(all.equal(new_procQ, z))
        ))
    if (any(!is_ok)) {
        warning("Processing queues from the submitted objects differ! ",
                "Dropping the processing queue.")
        new_procQ <- list()
    }
    ## processingData
    fls <- lapply(x, function(z) z@processingData@files)
    startidx <- cumsum(lengths(fls))
    ## featureData
    featd <- lapply(x, fData)
    ## Have to update the file index and the spectrum names.
    for (i in 2:length(featd)) {
        featd[[i]]$fileIdx <- featd[[i]]$fileIdx + startidx[i - 1]
        rownames(featd[[i]]) <- MSnbase:::formatFileSpectrumNames(
                                              fileIds = featd[[i]]$fileIdx,
                                              spectrumIds = featd[[i]]$spIdx,
                                              nSpectra = nrow(featd[[i]]),
                                              nFiles = length(unlist(fls))
                                          )
    }
    featd <- do.call(rbind, featd)
    featd$spectrum <- 1:nrow(featd)
    ## experimentData
    expdata <- lapply(x, function(z) {
        ed <- z@experimentData
        data.frame(instrumentManufacturer = ed@instrumentManufacturer,
                   instrumentModel = ed@instrumentModel,
                   ionSource = ed@ionSource,
                   analyser = ed@analyser,
                   detectorType = ed@detectorType,
                   stringsAsFactors = FALSE)
    })
    expdata <- do.call(rbind, expdata)
    expdata <- new("MIAPE",
                   instrumentManufacturer = expdata$instrumentManufacturer,
                   instrumentModel = expdata$instrumentModel,
                   ionSource = expdata$ionSource,
                   analyser = expdata$analyser,
                   detectorType = expdata$detectorType)

    ## protocolData
    protodata <- lapply(x, function(z) z@protocolData)
    if (any(unlist(lapply(protodata, nrow)) > 0))
        warning("Found non-empty protocol data, but merging protocol data is",
                " currently not supported. Skipped.")
    ## phenoData
    pdata <- do.call(rbind, lapply(x, pData))
    res <- new(
        "OnDiskMSnExp",
        phenoData = new("NAnnotatedDataFrame", data = pdata),
        featureData = new("AnnotatedDataFrame", featd),
        processingData = new("MSnProcess",
                             processing = paste0("Concatenated [", date(), "]"),
                             files = unlist(fls), smoothed = NA),
        experimentData = expdata,
        spectraProcessingQueue = new_procQ)
    if (validObject(res))
        res
}

#' @title Change the file path of an `OnDiskMSnExp` object
#'
#' @aliases dirname dirname,OnDiskMSnExp-method
#'
#' @name dirname
#'
#' @description
#'
#' `dirname` allows to get and set the path to the directory containing the
#' source files of the [OnDiskMSnExp-class] (or [XCMSnExp-class]) object.
#'
#' @param path [OnDiskMSnExp-class].
#'
#' @param value `character` of length 1 or length equal to the number of files
#'     defining the new path to the files.
#'
#' @md
#'
#' @author Johannes Rainer
setMethod("dirname", "OnDiskMSnExp", function(path) {
    dirname(fileNames(path))
})
#' @rdname dirname
setReplaceMethod("dirname", "OnDiskMSnExp", function(path, value) {
    flnms <- fileNames(path)
    if (length(value) == 1)
        value <- rep(value, length(flnms))
    new_flnms <- normalizePath(paste0(value, .Platform$file.sep,
                                      basename(flnms)))
    do_exist <- file.exists(new_flnms)
    if (any(!do_exist))
        stop("The following files do not exist: ",
             paste(new_flnms[!do_exist], ", "))
    path@processingData@files <- new_flnms
    validObject(path)
    path
})

#' @description
#'
#' Estimate the precursor intensity for an MS2 spectrum based on interpolation
#' using the intensity of the respective m/z peak from the previous and
#' following MS1 spectrum.
#'
#' @param x `OnDiskMSnExp` for a single file
#'
#' @param ppm `numeric(1)` with acceptable difference to MS1 m/z.
#'
#' @param method `character(1)` specifying how the precursor intensity should
#'     be estimated, either based on the previous MS1 scan
#'     (`method = "previous"`, the default) or using an interpolation between
#'     the previous and the subsequent MS1 scan (`method = "interpolation"`)
#'     considering also their retention time.
#'
#' @return `numeric` same length than `x` with the estimated precursor intensity
#'     or `NA` for MS1 spectra.
#'
#' @author Johannes Rainer
#'
#' @noRd
#'
#' @examples
#'
#' library(MSnbase)
#' fl <- system.file("TripleTOF-SWATH", "PestMix1_DDA.mzML", package = "msdata")
#' pest_dda <- readMSData(fl, mode = "onDisk")
#' res <- .estimate_prec_intensity(pest_dda)
#' fData(pest_dda)$precursorIntensity <- res
#'
#' ms2 <- filterMsLevel(pest_dda, msLevel = 2)
#' tic <- vapply(intensity(ms2), function(z) sum(z, na.rm = TRUE), numeric(1))
#' plot(tic, precursorIntensity(ms2))  # not that nice...
#'
#' fl <- proteomics(full.names = TRUE)[4]
#' tmt <- readMSData(fl, mode = "onDisk")
#'
#' res <- .estimate_prec_intensity(tmt, method = "interpolation")
#' res_2 <- .estimate_prec_intensity(tmt, method = "previous")
#'
#' par(mfrow = c(1, 2))
#' plot(res, precursorIntensity(tmt))
#' plot(res_2, precursorIntensity(tmt))
.estimate_prec_intensity <- function(x, ppm = 10,
                                     method = c("previous", "interpolation")) {
    method <- match.arg(method)
    pmz <- precursorMz(x)
    pmi <- rep(NA_real_, length(pmz))
    idx <- which(!is.na(pmz))
    x_ms1 <- filterMsLevel(x, msLevel = 1L)
    ms1_rt <- rtime(x_ms1)
    sps <- spectra(x_ms1)
    if (method == "previous") {
        for (i in idx) {
            ms2_rt <- .fdata(x)$retentionTime[i]
            ## Find the closest rtime before and the closest rtime after.
            before_idx <- which(ms1_rt < ms2_rt)
            before_int <- numeric()
            if (length(before_idx)) {
                sp <- sps[[before_idx[length(before_idx)]]]
                before_idx <- closest(pmz[i], sp@mz, ppm = ppm, tolerance = 0,
                                      duplicates = "closest")
                if (!is.na(before_idx)) {
                    before_rt <- sp@rt
                    before_int <- sp@intensity[before_idx]
                    before_int <- before_int[!is.na(before_int)]
                }
            }
            if (length(before_int))
                pmi[i] <- before_int
        }
    } else {
        for (i in idx) {
            ms2_rt <- .fdata(x)$retentionTime[i]
            ## Find the closest rtime before and the closest rtime after.
            before_idx <- which(ms1_rt < ms2_rt)
            before_int <- numeric()
            if (length(before_idx)) {
                sp <- sps[[before_idx[length(before_idx)]]]
                before_idx <- closest(pmz[i], sp@mz, ppm = ppm, tolerance = 0,
                                      duplicates = "closest")
                if (!is.na(before_idx)) {
                    before_rt <- sp@rt
                    before_int <- sp@intensity[before_idx]
                    before_int <- before_int[!is.na(before_int)]
                }
            }
            after_idx <- which(ms1_rt > ms2_rt)
            after_int <- numeric()
            if (length(after_idx)) {
                sp <- sps[[after_idx[1L]]]
                after_idx <- closest(pmz[i], sp@mz, ppm = ppm, tolerance = 0,
                                     duplicates = "closest")
                if (!is.na(after_idx)) {
                    after_rt <- sp@rt
                    after_int <- sp@intensity[after_idx]
                    after_int <- after_int[!is.na(after_int)]
                }
            }
            ## Check if we have before and after value
            if (length(before_int) && length(after_int)) {
                pmi[i] <- approx(c(before_rt, after_rt),
                                 c(before_int, after_int),
                                 xout = ms2_rt)$y
            } else {
                if (length(before_int))
                    pmi[i] <- before_int
                if (length(after_int))
                    pmi[i] <- after_int
            }
        }
    }
    pmi
}

#' @title Estimate precursor intensity for MS level 2 spectra
#'
#' @description
#'
#' `estimatePrecursorIntensity` determines the precursor intensity for a MS 2
#' spectrum based on the intensity of the respective signal from the
#' neighboring MS 1 spectra (i.e. based on the peak with the m/z matching the
#' precursor m/z of the MS 2 spectrum). Based on parameter `method` either the
#' intensity of the peak from the previous MS 1 scan is used
#' (`method = "previous"`) or an interpolation between the intensity from the
#' previous and subsequent MS1 scan is used (`method = "interpolation"`, which
#' considers also the retention times of the two MS1 scans and the retention
#' time of the MS2 spectrum).
#'
#' @param x `OnDiskMSnExp` or `XCMSnExp` object.
#'
#' @param ppm `numeric(1)` defining the maximal acceptable difference (in ppm)
#'     of the precursor m/z and the m/z of the corresponding peak in the MS 1
#'     scan.
#'
#' @param method `character(1)` defining the method how the precursor intensity
#'     should be determined (see description above for details). Defaults to
#'     `method = "previous"`.
#'
#' @param BPPARAM parallel processing setup. See [bpparam()] for details.
#'
#' @return `numeric` with length equal to the number of spectra in `x`. `NA` is
#'     returned for MS 1 spectra or if no matching peak in a MS 1 scan can be
#'     found for an MS 2 spectrum
#'
#' @author Johannes Rainer
#'
#' @md
estimatePrecursorIntensity <- function(x, ppm = 10,
                                       method = c("previous", "interpolation"),
                                       BPPARAM = bpparam()) {
    method <- match.arg(method)
    unlist(bplapply(.split_by_file2(x, subsetFeatureData = FALSE),
                    .estimate_prec_intensity, ppm = ppm, method = method,
                    BPPARAM = BPPARAM), use.names = FALSE)
}

#' Helper function to convert an OnDiskMSnExp to a Spectra object. This will
#' only convert the spectra data, but no sample information.
#'
#' @noRd
.OnDiskMSnExp2MsBackendMzR <- function(x) {
    .fData2MsBackendMzR(.fdata(x), fileNames(x))
}

.fData2MsBackendMzR <- function(x, filenames, res = new("MsBackendMzR")) {
    x$dataStorage <- x$dataOrigin <- filenames[x$fileIdx]
    colnames(x)[colnames(x) == "retentionTime"] <- "rtime"
    colnames(x)[colnames(x) == "seqNum"] <- "scanIndex"
    colnames(x)[colnames(x) == "precursorScanNum"] <- "precScanNum"
    colnames(x)[colnames(x) == "precursorMZ"] <- "precursorMz"
    colnames(x)[colnames(x) == "isolationWindowTargetMZ"] <-
        "isolationWindowTargetMz"
    colnames(x)[colnames(x) == "fileIdx"] <- "fromFile"
    x$isolationWindowLowerMz <- x$isolationWindowTargetMz -
        x$isolationWindowLowerOffset
    x$isolationWindowUpperMz <- x$isolationWindowTargetMz +
        x$isolationWindowUpperOffset
    x$isolationWindowUpperOffset <- NULL
    x$isolationWindowLowerOffset <- NULL
    slot(res, "spectraData", check = FALSE) <- as(x, "DataFrame")
    res
}
