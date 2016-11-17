## Test Param and Subclasses.

############################################################
## CentWaveParam
test_CentWaveParam <- function() {
    p <- new("CentWaveParam")
    checkDefaultValues <- function(x) {
        checkEquals(x@ppm, 25)
        checkEquals(x@peakwidth, c(20, 50))
        checkEquals(x@snthresh, 10)
        checkEquals(x@prefilter, c(3, 100))
        checkEquals(x@mzCenterFun, "wMean")
        checkEquals(x@integrate, 1)
        checkEquals(x@mzdiff, -0.001)
        checkEquals(x@fitgauss, FALSE)
        checkEquals(x@noise, 0)
        checkEquals(x@verboseColumns, FALSE)
        checkEquals(x@roiList, list())
        checkEquals(x@firstBaselineCheck, TRUE)
        checkEquals(x@roiScales, numeric())
    }
    checkDefaultValues(p)
    p <- new("CentWaveParam", fitgauss = TRUE)
    checkEquals(p@fitgauss, TRUE)
    p@fitgauss <- FALSE
    checkDefaultValues(p)
    p <- CentWaveParam()
    checkDefaultValues(p)
    p <- CentWaveParam(verboseColumns = TRUE, ppm = 30)
    checkEquals(p@verboseColumns, TRUE)
    checkEquals(p@ppm, 30)
    checkEquals(p@snthresh, 10)
    ## Check errors and validity
    checkException(CentWaveParam(ppm = -1))
    checkException(new("CentWaveParam", prefilter = c(3, 4, 5)))

    p@mzCenterFun <- "other"
    checkTrue(is.character(validObject(p)))
    ##plist <- as.list(p)
}

