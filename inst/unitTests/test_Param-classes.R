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

    ## Check getter/setter methods:
    ppm(p) <- 12
    checkEquals(ppm(p), 12)
    checkException(ppm(p) <- "a")
    peakwidth(p) <- c(5, 100)
    checkEquals(peakwidth(p), c(5, 100))
    checkException(peakwidth(p) <- c(1, 2, 3))
    snthresh(p) <- 300
    checkEquals(snthresh(p), 300)
    prefilter(p) <- c(1, 1000)
    checkEquals(prefilter(p), c(1, 1000))
    checkException(prefilter(p) <- 3)
    mzCenterFun(p) <- "meanApex3"
    checkEquals(mzCenterFun(p), "meanApex3")
    checkException(mzCenterFun(p) <- "other")
    integrate(p) <- 2
    checkEquals(integrate(p), 2L)
    checkException(integrate(p) <- 3L)
    mzdiff(p) <- 123.3
    checkEquals(mzdiff(p), 123.3)
    checkException(mzdiff(p) <- 1:4)
    fitgauss(p) <- TRUE
    checkEquals(fitgauss(p), TRUE)
    checkException(fitgauss(p) <- c(TRUE, TRUE))
    noise(p) <- 5
    checkEquals(noise(p), 5)
    checkException(noise(p) <- 1:4)
    verboseColumns(p) <- TRUE
    checkEquals(verboseColumns(p), TRUE)
    checkException(verboseColumns(p) c(TRUE, TRUE))
    firstBaselineCheck(b) <- TRUE
    checkEquals(firstBaselineCheck(b), TRUE)
    checkException(firstBaselineCheck(p) c(TRUE, TRUE))
}

############################################################
## MatchedFilterParam
test_matchedFilterParam <- function() {
    p <- new("MatchedFilterParam")
    checkDefaultValues <- function(x) {
        checkEquals(x@binSize, 0.1)
        checkEquals(x@impute, "none")
        checkEquals(x@baseValue, numeric())
        checkEquals(x@distance, numeric())
        checkEquals(x@fwhm, 30)
        checkEquals(x@sigma, 12.73994)
        checkEquals(x@max, 5)
        checkEquals(x@snthresh, 10)
        checkEquals(x@steps, 2)
        checkEquals(x@mzdiff, 0.6)
        checkEquals(x@index, FALSE)
    }
    checkDefaultValues(p)
    p <- new("MatchedFilterParam", binSize = 0.5)
    checkEquals(binSize(p), 0.5)
    ## binSize
    binSize(p) <- 0.6
    checkEquals(binSize(p), 0.6)
    p <- MatchedFilterParam(binSize = 0.6)
    checkEquals(binSize(p), 0.6)
    checkException(binSize(p) <- 1:5)
    ## impute
    impute(p) <- "lin"
    checkEquals(impute(p), "lin")
    p <- MatchedFilterParam(impute = "lin")
    checkEquals(impute(p), "lin")
    checkException(impute(p) <- "other")
    ## baseValue
    baseValue(p) <- 3
    checkEquals(baseValue(p), 3)
    p <- MatchedFilterParam(baseValue = 3)
    checkEquals(baseValue(p), 3)
    checkException(baseValue(p) <- 1:3)
    ## distance
    distance(p) <- 3
    checkEquals(distance(p), 3)
    p <- MatchedFilterParam(distance = 3)
    checkEquals(distance(p), 3)
    checkException(distance(p) <- 1:3)
    ## fwhm
    fwhm(p) <- 11
    checkEquals(fwhm(p), 11)
    p <- MatchedFilterParam(fwhm = 11)
    checkEquals(fwhm(p), 11)
    checkException(fwhm(p) <- 1:3)
    ## sigma
    sigma(p) <- 11
    checkEquals(sigma(p), 11)
    p <- MatchedFilterParam(sigma = 11)
    checkEquals(sigma(p), 11)
    checkException(sigma(p) <- 1:3)
    ## max
    max(p) <- 11
    checkEquals(max(p), 11)
    p <- MatchedFilterParam(max = 11)
    checkEquals(max(p), 11)
    checkException(max(p) <- 1:3)
    ## snthresh
    snthresh(p) <- 11
    checkEquals(snthresh(p), 11)
    p <- MatchedFilterParam(snthresh = 11)
    checkEquals(snthresh(p), 11)
    checkException(snthresh(p) <- 1:3)
    ## steps
    steps(p) <- 11
    checkEquals(steps(p), 11)
    p <- MatchedFilterParam(steps = 11)
    checkEquals(steps(p), 11)
    checkException(steps(p) <- 1:3)
    ## mzdiff
    mzdiff(p) <- 11
    checkEquals(mzdiff(p), 11)
    p <- MatchedFilterParam(mzdiff = 11)
    checkEquals(mzdiff(p), 11)
    checkException(mzdiff(p) <- 1:3)
    ## index
    index(p) <- TRUE
    checkEquals(index(p), TRUE)
    p <- MatchedFilterParam(index = TRUE)
    checkEquals(index(p), TRUE)
    checkException(index(p) <- 1:3)
}
