findPeaksMPI <- function(arg) {
    require(xcms)
    
    params <- arg$params
    myID <- arg$id
    if (is.null(params$method)) 
        params$method <- getOption("BioC")$xcms$findPeaks.method
    method <- match.arg(params$method, getOption("BioC")$xcms$findPeaks.methods)
    if (is.na(method))
        stop("unknown method : ", method)
    method <- paste("findPeaks", method, sep=".")    
    
    xRaw <- xcmsRaw(arg$file,profmethod=params$profmethod, profparam=params$profparam, profstep = 0)
    params["object"] <- xRaw
    
    ## remove parameters which are not used by method() from the parameter list
    params["method"] <- params["id"] <- params["profmethod"] <- params["profparam"] <- NULL
    
    peaks <- do.call(method, params)

    list(scantime=xRaw@scantime, peaks=cbind(peaks, sample = rep.int(myID, nrow(peaks))))
}
