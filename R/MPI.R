
##
## findPeaks slave function for parallel execution
##

findPeaksPar <- function(arg) {
    require(xcms)

    params <- arg$params
    myID <- arg$id
    if (is.null(params$method))
        params$method <- getOption("BioC")$xcms$findPeaks.method
    method <- match.arg(params$method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep=".")

    ## What about using the getXcmsRaw call here???
    xRaw <- xcmsRaw(arg$file, profmethod=params$profmethod, profparam=params$profparam, profstep = 0,
                    includeMSn=params$includeMSn, mslevel=params$mslevel, scanrange=params$scanrange)
    if(params$lockMassFreq == TRUE){
        xRaw<-stitch(xRaw, AutoLockMass(xRaw))
    }
    ## params["object"] <- xRaw

    ## remove parameters which are not used by method() from the parameter list
    params["method"] <- params["id"] <- params["profmethod"] <- params["profparam"] <- params["includeMSn"] <- params["lockMassFreq"] <-  params["mslevel"] <- NULL

    params["scanrange"] <- NULL ## avoid filtering scanrange twice, first in xRaw then in findPeaks

    peaks <- do.call(method, args = c(list(object = xRaw), params))

    list(scantime=xRaw@scantime, peaks=cbind(peaks, sample = rep.int(myID, nrow(peaks))))
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
        if (any(params$dataCorrection) == 1)
            lcraw <- stitch(lcraw, AutoLockMass(lcraw))
    }

    if (exists("params$polarity") && length(params$polarity) >0) {
        if (length(params$polarity) >0) {
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


                                        # Expanding the peakrange
    peakrange[,"mzmax"]  <-  peakrange[,"mzmax"]   +    (   (peakrange[,"mzmax"]-peakrange[,"mzmin"])/2    )*(expand.mz-1)
    peakrange[,"mzmin"]  <-  peakrange[,"mzmin"]   -    (   (peakrange[,"mzmax"]-peakrange[,"mzmin"])/2    )*(expand.mz-1)
    peakrange[,"rtmax"]  <-  peakrange[,"rtmax"]   +    (   (peakrange[,"rtmax"]-peakrange[,"rtmin"])/2    )*(expand.rt-1)
    peakrange[,"rtmin"]  <-  peakrange[,"rtmin"]   -    (   (peakrange[,"rtmax"]-peakrange[,"rtmin"])/2    )*(expand.rt-1)




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
