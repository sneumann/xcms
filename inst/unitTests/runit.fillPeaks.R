testFilledFlag <- function() {
    xsg <- group(faahko)
    xsgf <- fillPeaks(xsg, method="chrom")

    checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled), nrow(peaks(xsgf)))
}

testFillPeaksPar <- function() {
    xsg <- group(faahko)

    xsgfSerial <- fillPeaks(xsg, method="chrom")
    xsgfParallel <- fillPeaks(xsg, method="chrom", nSlaves=2)

    checkEqualsNumeric(nrow(peaks(xsgfSerial)),
                       nrow(peaks(xsgfParallel)))
}

## testFilledFlagMSW <- function() {

##   xsg <- group(ham)
##   xsgf <- fillPeaks(xsg)

##   checkEqualsNumeric(nrow(peaks(xsg)) + length(xsgf@filled, nrow(peaks(xsgf))) )
## }
