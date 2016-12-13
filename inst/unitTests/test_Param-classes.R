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


test_MSWParam <- function() {
    ## Check getter/setter methods:
    p <- new("MSWParam", snthresh = 14)
    checkEquals(snthresh(p), 14)
    snthresh(p) <- 13
    checkEquals(snthresh(p), 13)
    p <- MSWParam(snthresh = 21)
    checkEquals(snthresh(p), 21)
    checkException(MSWParam(snthresh = -4))

    p <- new("MSWParam", verboseColumns = TRUE)
    checkEquals(verboseColumns(p), TRUE)
    verboseColumns(p) <- FALSE
    checkEquals(verboseColumns(p), FALSE)
    p <- MSWParam(verboseColumns = TRUE)
    checkEquals(verboseColumns(p), TRUE)
    checkException(MSWParam(verboseColumns = c(FALSE, TRUE)))

    p <- new("MSWParam", scales = 1:5)
    checkEquals(scales(p), 1:5)
    scales(p) <- 6:12
    checkEquals(scales(p), 6:12)
    p <- MSWParam(scales = 1:3)
    checkEquals(scales(p), 1:3)

    p <- new("MSWParam", nearbyPeak = FALSE)
    checkEquals(nearbyPeak(p), FALSE)
    nearbyPeak(p) <- TRUE
    checkEquals(nearbyPeak(p), TRUE)
    p <- MSWParam(nearbyPeak = FALSE)
    checkEquals(nearbyPeak(p), FALSE)
    checkException(nearbyPeak(p) <- c(TRUE, FALSE))

    p <- new("MSWParam", peakScaleRange = 8)
    checkEquals(peakScaleRange(p), 8)
    peakScaleRange(p) <- 7
    checkEquals(peakScaleRange(p), 7)
    p <- MSWParam(peakScaleRange = 9)
    checkEquals(peakScaleRange(p), 9)
    checkException(peakScaleRange(p) <- 4:6)

    p <- new("MSWParam", ampTh = 0.4)
    checkEquals(ampTh(p), 0.4)
    ampTh(p) <- 7
    checkEquals(ampTh(p), 7)
    p <- MSWParam(ampTh = 9)
    checkEquals(ampTh(p), 9)
    checkException(ampTh(p) <- 4:6)

    p <- new("MSWParam", minNoiseLevel = 0.4)
    checkEquals(minNoiseLevel(p), 0.4)
    minNoiseLevel(p) <- 7
    checkEquals(minNoiseLevel(p), 7)
    p <- MSWParam(minNoiseLevel = 9)
    checkEquals(minNoiseLevel(p), 9)
    checkException(minNoiseLevel(p) <- -4)

    p <- new("MSWParam", ridgeLength = 14)
    checkEquals(ridgeLength(p), 14)
    ridgeLength(p) <- 7
    checkEquals(ridgeLength(p), 7)
    p <- MSWParam(ridgeLength = 9)
    checkEquals(ridgeLength(p), 9)
    checkException(ridgeLength(p) <- -4)

    p <- new("MSWParam", peakThr = 14)
    checkEquals(peakThr(p), 14)
    peakThr(p) <- 7
    checkEquals(peakThr(p), 7)
    p <- MSWParam(peakThr = 9)
    checkEquals(peakThr(p), 9)
    checkException(peakThr(p) <- 3:5)

    p <- new("MSWParam", tuneIn = FALSE)
    checkEquals(tuneIn(p), FALSE)
    tuneIn(p) <- TRUE
    checkEquals(tuneIn(p), TRUE)
    p <- MSWParam(tuneIn = FALSE)
    checkEquals(tuneIn(p), FALSE)
    checkException(tuneIn(p) <- c(TRUE, TRUE))

    p <- new("MSWParam", addParams = list(a = 3, b = 4))
    checkEquals(addParams(p), list(a = 3, b = 4))
    addParams(p) <- list(d = "a", e = 2:4)
    checkEquals(addParams(p), list(d = "a", e = 2:4))
    p <- MSWParam(otherp = 1:4, z = "a")
    checkEquals(addParams(p), list(otherp = 1:4, z = "a"))

    ## Check the .param2list method:
    p <- new("MSWParam", addParams = list(z = "z", bla = 1:4))
    L <- xcms:::.param2list(p)
    checkEquals(L$z, "z")
    checkEquals(L$bla, 1:4)
    checkEquals(L$snthresh, 3)
    p <- new("MSWParam")
    L <- xcms:::.param2list(p)
    checkTrue(!any(names(L) == "addParams"))
    checkEquals(L$snthresh, 3)
}

test_CentWavePredIsoParam <- function() {
    ## Check getter/setter methods:
    p <- new("CentWavePredIsoParam", ppm = 14)
    checkEquals(ppm(p), 14)
    ppm(p) <- 13
    checkEquals(ppm(p), 13)
    p <- CentWavePredIsoParam(ppm = 21)
    checkEquals(ppm(p), 21)
    checkException(CentWavePredIsoParam(ppm = -4))

    p <- new("CentWavePredIsoParam", peakwidth = c(1, 2))
    checkEquals(peakwidth(p), c(1, 2))
    peakwidth(p) <- c(2, 3)
    checkEquals(peakwidth(p), c(2, 3))
    p <- CentWavePredIsoParam(peakwidth = c(3, 4))
    checkEquals(peakwidth(p), c(3, 4))
    checkException(CentWavePredIsoParam(peakwidth = 1:3))

    p <- new("CentWavePredIsoParam", snthresh = 5)
    checkEquals(snthresh(p), 5)
    snthresh(p) <- 6
    checkEquals(snthresh(p), 6)
    p <- CentWavePredIsoParam(snthresh = 7)
    checkEquals(snthresh(p), 7)
    checkException(CentWavePredIsoParam(snthresh = 1:3))

    p <- new("CentWavePredIsoParam", prefilter = c(1, 2))
    checkEquals(prefilter(p), c(1, 2))
    prefilter(p) <- c(2, 3)
    checkEquals(prefilter(p), c(2, 3))
    p <- CentWavePredIsoParam(prefilter = c(3, 4))
    checkEquals(prefilter(p), c(3, 4))
    checkException(CentWavePredIsoParam(prefilter = 1:3))

    p <- new("CentWavePredIsoParam", mzCenterFun = "meanApex3")
    checkEquals(mzCenterFun(p), "meanApex3")
    mzCenterFun(p) <- "mean"
    checkEquals(mzCenterFun(p), "mean")
    p <- CentWavePredIsoParam(mzCenterFun = "mean")
    checkEquals(mzCenterFun(p), "mean")
    checkException(CentWavePredIsoParam(mzCenterFun = "median"))

    p <- new("CentWavePredIsoParam", integrate = 2L)
    checkEquals(integrate(p), 2L)
    integrate(p) <- 1L
    checkEquals(integrate(p), 1L)
    p <- CentWavePredIsoParam(integrate = 2L)
    checkEquals(integrate(p), 2L)
    checkException(CentWavePredIsoParam(integrate = 3L))

    p <- new("CentWavePredIsoParam", mzdiff = 1)
    checkEquals(mzdiff(p), 1)
    mzdiff(p) <- 2
    checkEquals(mzdiff(p), 2)
    p <- CentWavePredIsoParam(mzdiff = 3)
    checkEquals(mzdiff(p), 3)
    checkException(CentWavePredIsoParam(mzdiff = 2:3))

    p <- new("CentWavePredIsoParam", fitgauss = TRUE)
    checkEquals(fitgauss(p), TRUE)
    fitgauss(p) <- FALSE
    checkEquals(fitgauss(p), FALSE)
    p <- CentWavePredIsoParam(fitgauss = TRUE)
    checkEquals(fitgauss(p), TRUE)
    checkException(CentWavePredIsoParam(fitgauss = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", noise = 3)
    checkEquals(noise(p), 3)
    noise(p) <- 6
    checkEquals(noise(p), 6)
    p <- CentWavePredIsoParam(noise = 8)
    checkEquals(noise(p), 8)
    checkException(CentWavePredIsoParam(noise = c(3, 5)))

    p <- new("CentWavePredIsoParam", verboseColumns = TRUE)
    checkEquals(verboseColumns(p), TRUE)
    verboseColumns(p) <- FALSE
    checkEquals(verboseColumns(p), FALSE)
    p <- CentWavePredIsoParam(verboseColumns = TRUE)
    checkEquals(verboseColumns(p), TRUE)
    checkException(CentWavePredIsoParam(verboseColumns = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", firstBaselineCheck = TRUE)
    checkEquals(firstBaselineCheck(p), TRUE)
    firstBaselineCheck(p) <- FALSE
    checkEquals(firstBaselineCheck(p), FALSE)
    p <- CentWavePredIsoParam(firstBaselineCheck = FALSE)
    checkEquals(firstBaselineCheck(p), FALSE)
    checkException(CentWavePredIsoParam(firstBaselineCheck = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", snthreshIsoROIs = 5)
    checkEquals(snthreshIsoROIs(p), 5)
    snthreshIsoROIs(p) <- 6
    checkEquals(snthreshIsoROIs(p), 6)
    p <- CentWavePredIsoParam(snthreshIsoROIs = 7)
    checkEquals(snthreshIsoROIs(p), 7)
    checkException(CentWavePredIsoParam(snthreshIsoROIs = 1:3))

    p <- new("CentWavePredIsoParam", maxCharge = 5L)
    checkEquals(maxCharge(p), 5L)
    maxCharge(p) <- 6
    checkEquals(maxCharge(p), 6L)
    p <- CentWavePredIsoParam(maxCharge = 7)
    checkEquals(maxCharge(p), 7L)
    checkException(CentWavePredIsoParam(maxCharge = 1:3))

    p <- new("CentWavePredIsoParam", maxIso = 5L)
    checkEquals(maxIso(p), 5L)
    maxIso(p) <- 6
    checkEquals(maxIso(p), 6L)
    p <- CentWavePredIsoParam(maxIso = 7)
    checkEquals(maxIso(p), 7L)
    checkException(CentWavePredIsoParam(maxIso = 1:3))

    p <- new("CentWavePredIsoParam", mzIntervalExtension = FALSE)
    checkEquals(mzIntervalExtension(p), FALSE)
    mzIntervalExtension(p) <- TRUE
    checkEquals(mzIntervalExtension(p), TRUE)
    p <- CentWavePredIsoParam(mzIntervalExtension = FALSE)
    checkEquals(mzIntervalExtension(p), FALSE)
    checkException(CentWavePredIsoParam(mzIntervalExtension = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", polarity = "positive")
    checkEquals(polarity(p), "positive")
    polarity(p) <- "negative"
    checkEquals(polarity(p), "negative")
    p <- CentWavePredIsoParam(polarity = "negative")
    checkEquals(polarity(p), "negative")
    checkException(CentWavePredIsoParam(polarity = "bla"))

    ## Check the .param2list method:
    p <- new("CentWavePredIsoParam", snthresh = 123)
    L <- xcms:::.param2list(p)
    checkEquals(L$snthresh, 123)
}
