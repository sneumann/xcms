test.filenotfound <- function() {
    filename <- "rhabarber.mzData"
    checkTrue(xcms:::rampOpen(filename) < 0)
}

test.mzData <- function() {
    filename <- system.file('iontrap/extracted.mzData', package = "msdata")

    rampid <- xcms:::rampOpen(filename)
    checkTrue(rampid >= 0)

    rawdata <- xcms:::rampRawData(rampid)
    checkTrue(length(rawdata$rt) > 0)

    xcms:::rampClose(rampid)
}


