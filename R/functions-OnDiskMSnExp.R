## Functions for MSnbase's OnDiskMSnExp objects
#' @include do_detectFeatures-functions.R DataClasses.R


##' @param x an OnDiskMSnExp representing the whole experiment.
##' @param method The feature detection method to be used. Can be "centWave" etc.
##' @param param A class extending Param containing all parameters for the
##' feature detection method.
##'
##' @return a list of length 2, \code{peaks} containing a matrix with the
##' identified peaks and \code{date} the time stamp when the feature detection
##' was started.
##' @noRd
detectFeatures_OnDiskMSnExp <- function(object, method = "centWave",
                                        param) {
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    return(detectFeatures_Spectrum_list(x = spectra(object), method = method,
                                        param = param, rt = rtime(object)))
}


##' Run the feature detection on a list of Spectrum1 objects from the same
##' file
##'
##' @param x A list of Spectrum1 objects of a sample.
##' @param method The feature detection method to be used. Can be "centWave" etc.
##' @param param A class extending Param containing all parameters for the
##' feature detection method.
##' @param rt Numeric with the retention times for the spectra. If not provided
##' it is extracted from the spectra.
##' @return a list of length 2, \code{peaks} containing a matrix with the
##' identified peaks and \code{date} the time stamp when the feature detection
##' was started.
##' @author Johannes Rainer
##' @noRd
detectFeatures_Spectrum_list <- function(x, method = "centWave", param, rt) {
    method <- match.arg(method, c("centWave", "massifquant", "matchedFilter",
                                  "MSW", "centWaveWithPredIsoROIs"))
    method <- paste0("do_detectFeatures_", method)
    if (missing(param))
        stop("'param' has to be specified!")
    ## Check if the spectra are orderd by rt.
    if (missing(rt))
        rt <- unlist(lapply(x, rtime), use.names = FALSE)
    if (is.unsorted(rt))
        stop("Spectra are not ordered by retention time!")
    mzs <- lapply(x, mz)
    procDat <- date()
    return(list(peaks = do.call(method,
                                args = c(list(mz = unlist(mzs,
                                                          use.names = FALSE),
                                              int = unlist(lapply(x, intensity),
                                                           use.names = FALSE),
                                              valsPerSpect = lengths(mzs, FALSE),
                                              scantime = rt),
                                         as(param, "list"))),
                date = procDat))
}

## That's a special case since we don't expect to have rt available for this.
detectFeatures_MSW_OnDiskMSnExp <- function(object, method = "MSW",
                                            param) {
    if (missing(param))
        stop("'param' has to be specified!")
    ## pass the spectra to the _Spectrum_list function
    return(detectFeatures_MSW_Spectrum_list(x = spectra(object), method = method,
                                            param = param))
}
detectFeatures_MSW_Spectrum_list <- function(x, method = "MSW", param) {
    method <- match.arg(method, c("MSW"))
    method <- paste0("do_detectFeatures_", method)
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
##' @description Fill some settings and data from an OnDiskMSnExp or pSet into an
##' xcmsSet:
##' o fileNames
##' o phenoData
##' o rt
##' o mslevel
##' o scanrange: this should enable to subset the raw data again by scan index
##' (which should be equivalent to acquisitionNum/scanIndex)
##' @param pset The pSet from which data should be extracted
##' @noRd
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

##' @description Processes the result list returned by an lapply/bplapply to
##' detectFeatures_Spectrum_list or detectFeatures_OnDiskMSnExp and returns a
##' list with two elements: \code{$peaks} the peaks matrix of identified
##' features and \code{$procHist} a list of ProcessHistory objects (empty if
##' \code{getProcHist = FALSE}).
##' @param x See description above.
##' @param getProcHist Wheter ProcessHistory objects should be returned too.
##' @param fnames The file names.
##' @return See description above.
##' @author Johannes Rainer
##' @noRd
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
            phList[[i]] <- ProcessHistory(info. = paste0("Feature detection in '",
                                                         basename(fnames[i]),
                                                         "': ", n_pks,
                                                         " features identified."),
                                          date. = x[[i]]$date,
                                          type. = .PROCSTEP.FEATURE.DETECTION,
                                          fileIndex. = i
                                          )
    }
    return(list(peaks = pks, procHist = phList))
}
