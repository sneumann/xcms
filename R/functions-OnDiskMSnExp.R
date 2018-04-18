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
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    ## Since we're calling this function already with bplapply ensure that
    ## the spectra call is not firing its own parallel processing!
    return(findChromPeaks_Spectrum_list(x = spectra(object,
                                                    BPPARAM = SerialParam()),
                                        method = method,
                                        param = param, rt = rtime(object)))
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
    if (method == "MSW")
        method <- paste0("do_findPeaks_", method)
    if (missing(param))
        stop("'param' has to be specified!")
    ## Check if the spectra are orderd by rt.
    if (missing(rt))
        rt <- unlist(lapply(x, rtime), use.names = FALSE)
    if (is.unsorted(rt))
        stop("Spectra are not ordered by retention time!")
    mzs <- lapply(x, mz)
    vals_per_spect <- lengths(mzs, FALSE)
    procDat <- date()
    res <- do.call(method, args = c(list(mz = unlist(mzs,
                                                     use.names = FALSE),
                                         int = unlist(lapply(x, intensity),
                                                      use.names = FALSE),
                                         valsPerSpect = vals_per_spect,
                                         scantime = rt),
                                    as(param, "list")))
    ## Ensure that we call the garbage collector to eventually clean unused stuff
    rm(mzs)
    rm(x)
    rm(rt)
    gc()
    return(list(peaks = res, date = procDat))
}

## That's a special case since we don't expect to have rt available for this.
findPeaks_MSW_OnDiskMSnExp <- function(object, method = "MSW",
                                            param) {
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    return(findPeaks_MSW_Spectrum_list(x = spectra(object,
                                                   BPPARAM = SerialParam()),
                                       method = method,
                                       param = param))
}
findPeaks_MSW_Spectrum_list <- function(x, method = "MSW", param) {
    method <- match.arg(method, c("MSW"))
    method <- paste0("do_findPeaks_", method)
    if (missing(param))
        stop("'param' has to be specified!")
    mzs <- lapply(x, mz)
    procDat <- date()
    return(list(peaks = do.call(method,
                                args = c(list(mz = unlist(mzs,
                                                          use.names = FALSE),
                                              int = unlist(lapply(x, intensity),
                                                           use.names = FALSE)
                                              ),
                                         as(param, "list"))),
                date = procDat))
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
    rt <- split(unname(rtime(pset)), f = fromFile(pset))
    object@rt <- list(raw = rt, corrected = rt)
    ## mslevel
    mslevel(object) <- unique(msLevel(pset))
    ## scanrange
    ## scanrange(object) <- range(acquisitionNum(pset))
    ## Using scanIndex instead of acquisitionNum as the scan index should
    ## represent the index of the spectrum within the file, while acquisitionNum
    ## might be different e.g. if the mzML was subsetted.
    scanrange(object) <- range(scanIndex(pset))
    return(object)
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
    return(list(peaks = pks, procHist = phList))
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
    nSamples <- length(fileNames(object))
    if (nSamples <= 1)
        stop("Can not perform a retention time correction on less than to",
             " files.")
    
    ## centerSample
    if (length(centerSample(param))) {
        if (!(centerSample(param) %in% 1:nSamples))
            stop("'centerSample' has to be a single integer between 1 and ",
                 nSamples, "!")
    } else {
        centerSample(param) <- floor(median(1:nSamples))
    }
    message("Sample number ", centerSample(param), " used as center sample.")

    ## Get the profile matrix of the center sample:
    ## Using the (hidden) parameter returnBreaks to return also the breaks of
    ## the bins of the profile matrix. I can use them to align the matrices
    ## later.
    ## NOTE: it might be event better to just re-use the breaks from the center
    ## sample for the profile matrix generation of all following samples.
    suppressMessages(
        profCtr <- profMat(object, method = "bin", step = binSize(param),
                           fileIndex = centerSample(param),
                           returnBreaks = TRUE)[[1]]
    )
    ## Now split the object by file
    objL <- splitByFile(object, f = factor(seq_len(nSamples)))
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
        ## ---------------------------------------
        ## 1)Check the scan times of both objects:
        scantime1 <- unname(rtime(cntr))
        scantime2 <- unname(rtime(z))
        ## median difference between spectras' scan time.
        mstdiff <- median(c(diff(scantime1), diff(scantime2)))

        ## rtup1 <- seq_along(scantime1)
        ## rtup2 <- seq_along(scantime2)

        mst1 <- which(diff(scantime1) > 5 * mstdiff)[1]
        if (!is.na(mst1)) {
            scantime1 <- scantime1[seq_len((mst1 - 1))]
            message("Found gaps in scan times of the center sample: cut ",
                    "scantime-vector at ", scantime1[mst1]," seconds.")
        }
        mst2 <- which(diff(scantime2) > 5 * mstdiff)[1]
        if(!is.na(mst2)) {
            scantime2 <- scantime2[seq_len((mst2 - 1))]
            message("Found gaps in scan time of file ", basename(fileNames(z)),
                    ": cut scantime-vector at ", scantime2[mst2]," seconds.")
        }
        ## Drift of measured scan times - expected to be largest at the end.
        rtmaxdiff <- abs(diff(c(scantime1[length(scantime1)],
                                scantime2[length(scantime2)])))
        ## If the drift is larger than the threshold, cut the matrix up to the
        ## max allowed difference.
        if(rtmaxdiff > (5 * mstdiff)){
            rtmax <- min(scantime1[length(scantime1)],
                         scantime2[length(scantime2)])
            scantime1 <- scantime1[scantime1 <= rtmax]
            scantime2 <- scantime2[scantime2 <= rtmax]
        }
        valscantime1 <- length(scantime1)
        valscantime2 <- length(scantime2)
        ## Finally, restrict the profile matrix to columns 1:valscantime
        if (ncol(cntrPr$profMat) > valscantime1) {
            cntrPr$profMat <- cntrPr$profMat[, -c((valscantime1 + 1):
                                                  ncol(cntrPr$profMat))]
        }
        if(ncol(curP$profMat) > valscantime2) {
            curP$profMat <- curP$profMat[, -c((valscantime2 + 1):
                                              ncol(curP$profMat))]
        }
        ## ---------------------------------
        ## 2) Now match the breaks/mz range.
        ##    The -1 below is because the breaks define the upper and lower
        ##    boundary. Have to do it that way to be in line with the orignal
        ##    code... would be better to use the breaks as is.
        mzr1 <- c(cntrPr$breaks[1], cntrPr$breaks[length(cntrPr$breaks) - 1])
        mzr2 <- c(curP$breaks[1], curP$breaks[length(curP$breaks) - 1])
        mzmin <- min(c(mzr1[1], mzr2[1]))
        mzmax <- max(c(mzr1[2], mzr2[2]))
        mzs <- seq(mzmin, mzmax, by = binSize(parms))
        ## Eventually add empty rows at the beginning
        if (mzmin < mzr1[1]) {
            tmp <- matrix(0, (length(seq(mzmin, mzr1[1], binSize(parms))) - 1),
                          ncol = ncol(cntrPr$profMat))
            cntrPr$profMat <- rbind(tmp, cntrPr$profMat)
        }
        ## Eventually add empty rows at the end
        if (mzmax > mzr1[2]) {
            tmp <- matrix(0, (length(seq(mzr1[2], mzmax, binSize(parms))) - 1),
                          ncol = ncol(cntrPr$profMat))
            cntrPr$profMat <- rbind(cntrPr$profMat, tmp)
        }
        ## Eventually add empty rows at the beginning
        if (mzmin < mzr2[1]) {
            tmp <- matrix(0, (length(seq(mzmin, mzr2[1], binSize(parms))) - 1),
                          ncol = ncol(curP$profMat))
            curP$profMat <- rbind(tmp, curP$profMat)
        }
        ## Eventually add empty rows at the end
        if (mzmax > mzr2[2]) {
            tmp <- matrix(0, (length(seq(mzr2[2], mzmax, binSize(parms))) - 1),
                          ncol = ncol(curP$profMat))
            curP$profMat <- rbind(curP$profMat, tmp)
        }
        ## A final check of the data.
        mzvals <- length(mzs)
        cntrVals <- length(cntrPr$profMat)
        curVals <- length(curP$profMat)
        if ((mzvals * valscantime1) != cntrVals | (mzvals * valscantime2) != curVals)
            ## Here the question is if we REALLY need to have the same numbers
            ## of values in both. This caused the problems in issue #196
            ## | cntrVals != curVals)
            stop("Dimensions of profile matrices of files ",
                 basename(fileNames(cntr)), " and ", basename(fileNames(z)),
                 " do not match!")
        ## Done with preparatory stuff - now I can perform the alignment.
        rtadj <- .Call("R_set_from_xcms", valscantime1, scantime1, mzvals, mzs,
                       cntrPr$profMat, valscantime2, scantime2, mzvals, mzs,
                       curP$profMat, response(parms), distFun(parms),
                       gapInit(parms), gapExtend(parms), factorDiag(parms),
                       factorGap(parms), as.numeric(localAlignment(parms)),
                       initPenalty(parms))
        if (length(rtime(z)) > valscantime2) {
            ## Adding the raw retention times if we were unable to align all of
            ## them.
            rtadj <- c(rtadj, rtime(z)[(valscantime2 + 1):length(rtime(z))])
            warning(basename(fileNames(z)), " :could only align up to a ",
                    "retention time of ", rtime(z)[valscantime2], " seconds. ",
                    "After that raw retention times are reported.")
        }
        message("OK")
        return(rtadj)
        ## Related to issue #122: try to resemble the rounding done in the
        ## recor.obiwarp method.
        ## return(round(rtadj, 2))
    }, cntr = centerObject, cntrPr = profCtr, parms = param)
    ## Add also the rtime of the center sample:
    adjRt <- vector("list", nSamples)
    adjRt[centerSample(param)] <- list(unname(rtime(centerObject)))
    ## Add the result.
    idxs <- 1:nSamples
    idxs <- idxs[idxs != centerSample(param)]
    adjRt[idxs] <- res
    return(adjRt)
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
