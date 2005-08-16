mzXMLIsFile <- function(filename) {

    if (!file.exists(filename))
        return(FALSE)
    text <- readChar(filename, 500)
    
    length(text) > 0
}

mzXMLOpen <- function(filename) {

    result <- .C("MzXMLOpen",
                 as.character(filename),
                 mxid = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")
    
    if (result$status)
        return(result$status)
    
    return(result$mxid)
}

mzXMLClose <- function(mxid) {

    result <- .C("MzXMLClose",
                 as.integer(mxid),
                 status = integer(1),
                 PACKAGE = "xcms")
    
    result$status
}

mzXMLNumScans <- function(mxid) {

    result <- .C("MzXMLNumScans",
                 as.integer(mxid),
                 numscans = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")
    
    if (result$status)
        return(result$status)
    
    return(result$numscans)
}

mzXMLRawData <- function(mxid) {

    len <- mzXMLNumScans(mxid)
    if (len < 1)
        stop("Couldn't read any scans")
    
    result1 <- .C("MzXMLRTTICPeaks",
                  as.integer(mxid),
                  rt = double(len),
                  tic = double(len),
                  peaksTotal = integer(1),
                  status = integer(1),
                  PACKAGE = "xcms")
    
    if (result1$status)
        stop("Couldn't read retention times/total ion current")
    
    result2 <- .C("MzXMLSIPeaks",
                  as.integer(mxid),
                  scanindex = integer(len),
                  mz = double(result1$peaksTotal),
                  intensity = double(result1$peaksTotal),
                  status = integer(1),
                  DUP = FALSE, PACKAGE = "xcms")
    
    if (result2$status)
        stop("Couldn't read mass/intensity values")
    
    return(list(rt = result1$rt, tic = result1$tic, 
                scanindex = result2$scanindex, mz = result2$mz,
                intensity = result2$intensity))
}
