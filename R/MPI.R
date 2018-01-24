##
## findPeaks slave function for parallel execution
##

############################################################
## findPeaksPar
##
## This should at some point be replaced by a call that does not need
## parameter lists and getting options from the environment.
findPeaksPar <- function(arg) {
    require(xcms)

    procDate <- date()
    params <- arg$params
    myID <- arg$id
    if (is.null(params$method))
        params$method <- getOption("BioC")$xcms$findPeaks.method
    method <- match.arg(params$method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep=".")

    ## What about using the getXcmsRaw call here???
    xRaw <- xcmsRaw(arg$file, profmethod=params$profmethod,
                    profparam=params$profparam, profstep = 0,
                    includeMSn=params$includeMSn, mslevel=params$mslevel,
                    scanrange=params$scanrange)
    if(params$lockMassFreq == TRUE){
        xRaw <- stitch(xRaw, AutoLockMass(xRaw))
    }

    ## remove parameters which are not used by method() from the parameter list
    params["method"] <- NULL
    params["id"] <- NULL
    params["profmethod"] <- NULL
    params["profparam"] <- NULL
    params["includeMSn"] <- NULL
    params["lockMassFreq"] <- NULL
    params["mslevel"] <- NULL
    params["scanrange"] <- NULL ## avoid filtering scanrange twice.

    peaks <- do.call(method, args = c(list(object = xRaw), params))

    ## Ensure to remove data to avoid memory accumulation.
    scanT <- xRaw@scantime
    rm(xRaw)
    gc()
    peaks <- cbind(peaks, sample = rep.int(myID, nrow(peaks)))
    ## Ensure that last column is named "sample" even if we didn't find
    ## anything (issue #220)
    colnames(peaks)[ncol(peaks)] <- "sample"
    list(scantime = scanT,
         peaks = peaks,
         date = procDate)
}

############################################################
## findChromPeaks
##
## Same as findPeaksPar but without the need to pass argument lists
## and read settings from the global options.
## args should be a list with arguments
## o file: the file name
## o readParams: parameter class to read the file; actually we would only
##   need the scanrange, the includeMSn and the lockMassFreq here.
## o detectParams: parameter class for the peak detection.
findChromPeaksInFile <- function(args) {
    ## Placeholder
}


##
## findPeaks slave function for parallel execution
##

fillPeaksChromPar <- function(arg) {

    require(xcms)

    params <- arg$params
    myID <- arg$id
    cat(arg$file, "\n")

    prof <- params$prof
    rtcor <- params$rtcor
    peakrange <- params$peakrange
    expand.mz <- params$expand.mz
    expand.rt <- params$expand.rt
    gvals <- params$gvals$gvals

    lcraw <- xcmsRaw(arg$file, profmethod=params$prof$method, profstep = 0)

    if (length(params$dataCorrection) > 1) {
        ## Note: dataCorrection (as set in the xcmsSet function) is either
        ## 1 for all or for none.
        if (any(params$dataCorrection == 1))
            lcraw <- stitch(lcraw, AutoLockMass(lcraw))
    }

    if (exists("params$polarity") && length(params$polarity) >0) {
        if (length(params$polarity) > 0) {
            ## Retain wanted polarity only
            lcraws <- split(lcraw, lcraw@polarity, DROP=TRUE)
            lcraw <- lcraws[[params$polarity]]
        }
    }

    if (length(prof) > 2)
        lcraw@profparam <- prof[seq(3, length(prof))]
    if (length(rtcor) == length(lcraw@scantime) ) {
        lcraw@scantime <- rtcor
    } else {
        warning("(corrected) retention time vector length mismatch for ", basename(arg$file))
    }


    ## Expanding the peakrange
    incrMz <- (peakrange[, "mzmax"] - peakrange[, "mzmin"]) / 2 * (expand.mz - 1)
    peakrange[, "mzmax"] <- peakrange[, "mzmax"] + incrMz
    peakrange[, "mzmin"] <- peakrange[, "mzmin"] - incrMz
    incrRt <- (peakrange[, "rtmax"] - peakrange[, "rtmin"]) / 2 * (expand.rt - 1)
    peakrange[, "rtmax"] <- peakrange[, "rtmax"] + incrRt
    peakrange[, "rtmin"] <- peakrange[, "rtmin"] - incrRt
    
    naidx <- which(is.na(gvals[,myID]))

    newpeaks <- getPeaks(lcraw, peakrange[naidx,,drop=FALSE], step = prof$step)

    list(myID=myID, newpeaks=cbind(newpeaks, sample=myID))
}



msgfun.featureDetection <- function(x,i) {
    message("Detecting features in file #",i,":",basename(x[[i]]$file))
    flush.console();
}

msgfunGeneric <- function(x, i) {
    message(i,":",basename(x[[i]]$file))
    flush.console();
}
