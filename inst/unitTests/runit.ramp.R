#test.filenotfound <- function() {
#    filename <- "rhabarber.mzData"
#    checkTrue(xcms:::rampOpen(filename) < 0)
#}

#test.mzData <- function() {
#    filename <- system.file('iontrap/extracted.mzData', package = "msdata")
#
#    result <- .C("RampROpen", filename,
#       rampid = integer(1), status = integer(1),
#       PACKAGE = "xcms")
#
#    checkTrue(result[[1]] == filename)
#    checkTrue(result$rampid >= 0)
#    checkEqualsNumeric(result$status, 0)
#
#    xcms:::rampClose(result$rampid)
#}
