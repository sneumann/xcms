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
    ppm(p) <- 1
    checkEquals(ppm(p), 1)
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
    checkException(verboseColumns(p) <- c(TRUE, TRUE))
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
    distance(p) <- 4
    checkEquals(distance(p), 4)
    p <- MatchedFilterParam(distance = 4)
    checkEquals(distance(p), 4)
    checkException(distance(p) <- 1:3)
    ## fwhm
    fwhm(p) <- 11
    checkEquals(fwhm(p), 11)
    p <- MatchedFilterParam(fwhm = 11)
    checkEquals(fwhm(p), 11)
    checkException(fwhm(p) <- 1:3)
    ## sigma
    sigma(p) <- 12
    checkEquals(sigma(p), 12)
    p <- MatchedFilterParam(sigma = 12)
    checkEquals(sigma(p), 12)
    checkException(sigma(p) <- 1:3)
    ## max
    max(p) <- 13
    checkEquals(max(p), 13)
    p <- MatchedFilterParam(max = 13)
    checkEquals(max(p), 13)
    checkException(max(p) <- 1:3)
    ## snthresh
    snthresh(p) <- 14
    checkEquals(snthresh(p), 14)
    p <- MatchedFilterParam(snthresh = 14)
    checkEquals(snthresh(p), 14)
    checkException(snthresh(p) <- 1:3)
    ## steps
    steps(p) <- 15
    checkEquals(steps(p), 15)
    p <- MatchedFilterParam(steps = 15)
    checkEquals(steps(p), 15)
    checkException(steps(p) <- 1:3)
    ## mzdiff
    mzdiff(p) <- 16
    checkEquals(mzdiff(p), 16)
    p <- MatchedFilterParam(mzdiff = 16)
    checkEquals(mzdiff(p), 16)
    checkException(mzdiff(p) <- 1:3)
    ## index
    index(p) <- TRUE
    checkEquals(index(p), TRUE)
    p <- MatchedFilterParam(index = TRUE)
    checkEquals(index(p), TRUE)
    checkException(index(p) <- 1:3)
}

test_MassifquantParam <- function() {
    ## Check getter/setter methods:
    p <- new("MassifquantParam")
    ppm(p) <- 1
    checkEquals(ppm(p), 1)
    p <- MassifquantParam(ppm = 1)
    checkEquals(ppm(p), 1)
    checkException(ppm(p) <- 1:4)
    peakwidth(p) <- c(5, 100)
    checkEquals(peakwidth(p), c(5, 100))
    p <- MassifquantParam(peakwidth = c(5, 100))
    checkEquals(peakwidth(p), c(5, 100))
    checkException(peakwidth(p) <- c(1, 2, 3))
    snthresh(p) <- 300
    checkEquals(snthresh(p), 300)
    p <- MassifquantParam(snthresh = 300)
    checkEquals(snthresh(p), 300)
    checkException(snthresh(p) <- 2:4)
    prefilter(p) <- c(1, 1000)
    checkEquals(prefilter(p), c(1, 1000))
    p <- MassifquantParam(prefilter = c(1, 1000))
    checkException(prefilter(p) <- 3)
    mzCenterFun(p) <- "meanApex3"
    checkEquals(mzCenterFun(p), "meanApex3")
    p <- MassifquantParam(mzCenterFun = "meanApex3")
    checkEquals(mzCenterFun(p), "meanApex3")
    checkException(mzCenterFun(p) <- "other")
    integrate(p) <- 2
    checkEquals(integrate(p), 2L)
    p <- MassifquantParam(integrate = 2)
    checkEquals(integrate(p), 2L)
    checkException(integrate(p) <- 3L)
    mzdiff(p) <- 123.3
    checkEquals(mzdiff(p), 123.3)
    p <- MassifquantParam(mzdiff = 123.3)
    checkEquals(mzdiff(p), 123.3)
    checkException(mzdiff(p) <- 1:4)
    fitgauss(p) <- TRUE
    checkEquals(fitgauss(p), TRUE)
    p <- MassifquantParam(fitgauss = TRUE)
    checkEquals(fitgauss(p), TRUE)
    checkException(fitgauss(p) <- c(TRUE, TRUE))
    noise(p) <- 5
    checkEquals(noise(p), 5)
    p <- MassifquantParam(noise = 5)
    checkEquals(noise(p), 5)
    checkException(noise(p) <- 1:4)
    verboseColumns(p) <- TRUE
    checkEquals(verboseColumns(p), TRUE)
    p <- MassifquantParam(verboseColumns = TRUE)
    checkEquals(verboseColumns(p), TRUE)
    checkException(verboseColumns(p) <- c(TRUE, TRUE))
    criticalValue(p) <- 34
    checkEquals(criticalValue(p), 34)
    p <- MassifquantParam(criticalValue = 34)
    checkEquals(criticalValue(p), 34)
    checkException(criticalValue(p) <- as.numeric(3:5))
    consecMissedLimit(p) <- 4
    checkEquals(consecMissedLimit(p), 4L)
    p <- MassifquantParam(consecMissedLimit = 4)
    checkEquals(consecMissedLimit(p), 4L)
    checkException(consecMissedLimit(p) <- 3:5)
    ## unions, checkBack, withWave
    unions(p) <- 0
    checkEquals(unions(p), 0L)
    p <- MassifquantParam(unions = 1)
    checkEquals(unions(p), 1L)
    checkException(unions(p) <- 4)
    checkBack(p) <- 0
    checkEquals(checkBack(p), 0L)
    p <- MassifquantParam(checkBack = 1)
    checkEquals(checkBack(p), 1L)
    checkException(checkBack(p) <- 4)
    withWave(p) <- TRUE
    checkEquals(withWave(p), TRUE)
    p <- MassifquantParam(withWave = TRUE)
    checkEquals(withWave(p), TRUE)
    checkException(withWave(p) <- c(TRUE, FALSE))
}
