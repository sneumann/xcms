rampInit <- function() {

    result <- .C("RampRInit",
                 PACKAGE = "xcms")
}

rampPrintFiles <- function() {

    result <- .C("RampRPrintFiles",
                 PACKAGE = "xcms")
}

rampIsFile <- function(filename) {

    # The C version doesn't do anything extra
    #.C("RampRIsFile",
    #   as.character(filename),
    #   isfile = logical(1),
    #   PACKAGE = "xcms")$isfile

    if (!file.exists(filename))
        return(FALSE)
    text <- readChar(filename, 1024)

    length(text) > 0
}

rampOpen <- function(filename) {

    result <- .C("RampROpen",
                 as.character(filename),
                 rampid = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(result$status)

    return(result$rampid)
}

rampClose <- function(rampid) {

    result <- .C("RampRClose",
                 as.integer(rampid),
                 PACKAGE = "xcms")
}

rampCloseAll <- function() {

    result <- .C("RampRCloseAll",
                 PACKAGE = "xcms")
}

rampNumScans <- function(rampid) {

    result <- .C("RampRNumScans",
                 as.integer(rampid),
                 numscans = integer(1),
                 status = integer(1),
                 PACKAGE = "xcms")

    if (result$status)
        return(NA)

    return(result$numscans)
}

rampScanHeaders <- function(rampid) {

    .Call("RampRScanHeaders",
          as.integer(rampid),
          PACKAGE = "xcms")
}

rampSIPeaks <- function(rampid, seqNum, peaksCount) {

    if (!is.integer(seqNum))
        seqNum <- as.integer(seqNum)
    if (!is.integer(peaksCount))
        peaksCount <- as.integer(peaksCount)
    .Call("RampRSIPeaks",
          as.integer(rampid),
          seqNum,
          peaksCount,
          PACKAGE = "xcms")
}

rampRawData <- function(rampid) {

    scanHeaders <- rampScanHeaders(rampid)

    # Some of these checks work around buggy RAMP indexing code
    scans <- scanHeaders$msLevel == 1 & scanHeaders$seqNum > 0 &
             !duplicated(scanHeaders$acquisitionNum) &
             scanHeaders$peaksCount > 0
    if ("Full" %in% levels(scanHeaders$scanType))
        scans <- scans & scanHeaders$scanType == "Full"

    scans <- which(scans)

    sipeaks <- rampSIPeaks(rampid, scans, scanHeaders$peaksCount[scans])

    return(list(rt = scanHeaders$retentionTime[scans],
                acquisitionNum = scanHeaders$acquisitionNum[scans],
                tic = scanHeaders$totIonCurrent[scans],
                scanindex = sipeaks$scanindex, mz = sipeaks$mz,
                intensity = sipeaks$intensity))
}

rampRawDataMSn <- function(rampid) {

    # Check if we have MSn at all
    scanHeaders <- rampScanHeaders(rampid)
    if (max(scanHeaders[,"msLevel"]) < 2) {
        warning("MSn spectra requested but not found")
        return (NULL);
    }

    # Some of these checks work around buggy RAMP indexing code
    scans <- ( scanHeaders$msLevel >= 2 & scanHeaders$seqNum > 0
              & !duplicated(scanHeaders$acquisitionNum)
              & scanHeaders$peaksCount > 0)

    scans <- which(scans)

    sipeaks <- rampSIPeaks(rampid, scans, scanHeaders$peaksCount[scans])

    retdata <- list(rt = scanHeaders$retentionTime[scans],
                    acquisitionNum = scanHeaders$acquisitionNum[scans],
                    precursorNum=scanHeaders$precursorScanNum[scans],
                    precursorMZ = scanHeaders$precursorMZ[scans],
                    precursorIntensity = scanHeaders$precursorIntensity[scans],
                    peaksCount=scanHeaders$peaksCount[scans],
                    msLevel = scanHeaders$msLevel[scans],
                    precursorCharge = scanHeaders$precursorCharge[scans],
                    scanindex = sipeaks$scanindex, collisionEnergy = scanHeaders$collisionEnergy[scans],
                    mz = sipeaks$mz,
                    intensity =sipeaks$intensity);

    return(retdata)
}
