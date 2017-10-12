## Test Param and Subclasses.

############################################################
## CentWaveParam
test_CentWaveParam <- function() {
    ##checkTrue(FALSE)
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
    checkException(validObject(p))
    ##plist <- as.list(p)

    ## Check getter/setter methods:
    p <- new("CentWaveParam")
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
    firstBaselineCheck(p) <- TRUE
    checkEquals(firstBaselineCheck(p), TRUE)
    checkException(firstBaselineCheck(p) <- c(TRUE, TRUE))
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

test_PeakDensityParam <- function() {
    ## Check getter/setter methods:
    p <- new("PeakDensityParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    checkEquals(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    checkEquals(sampleGroups(p), 1:4)
    p <- PeakDensityParam(sampleGroups = c("a", "a", "b"))
    checkEquals(sampleGroups(p), c("a", "a", "b"))
    checkException(sampleGroups(p) <- NULL)
    checkException(sampleGroups(p) <- c(2, 2, NA))
    checkException(PeakDensityParam())
    checkException(PeakDensityParam(sampleGroups = c(1, 1, NA)))
    
    p <- new("PeakDensityParam", bw = 3)
    checkEquals(bw(p), 3)
    bw(p) <- 20
    checkEquals(bw(p), 20)
    p <- PeakDensityParam(bw = 33, sampleGroups = "a")
    checkEquals(bw(p), 33)
    checkException(PeakDensityParam(bw = -4, sampleGroups = 1))

    ## minFraction
    p <- new("PeakDensityParam", minFraction = 0.7)
    checkEquals(minFraction(p), 0.7)
    minFraction(p) <- 0.2
    checkEquals(minFraction(p), 0.2)
    p <- PeakDensityParam(minFraction = 0.4, sampleGroups = 3)
    checkEquals(minFraction(p), 0.4)
    checkException(PeakDensityParam(minFraction = -4, sampleGroups = 3))
    checkException(minFraction(p) <- c(0.3, 0.2, 0.4))
    checkException(minFraction(p) <- 2)

    ## minSamples
    p <- new("PeakDensityParam", minSamples = 3)
    checkEquals(minSamples(p), 3)
    minSamples(p) <- 20
    checkEquals(minSamples(p), 20)
    p <- PeakDensityParam(minSamples = 33, sampleGroups = rep(3, 100))
    checkEquals(minSamples(p), 33)
    checkException(PeakDensityParam(minSamples = -4, sampleGroups = 4))

    ## binSize
    p <- new("PeakDensityParam", binSize = 3)
    checkEquals(binSize(p), 3)
    binSize(p) <- 20
    checkEquals(binSize(p), 20)
    p <- PeakDensityParam(binSize = 0.3, sampleGroups = 4)
    checkEquals(binSize(p), 0.3)
    checkException(PeakDensityParam(binSize = -4, sampleGroups = 4))

    ## maxFeatures
    p <- new("PeakDensityParam", maxFeatures = 3)
    checkEquals(maxFeatures(p), 3)
    maxFeatures(p) <- 20
    checkEquals(maxFeatures(p), 20)
    p <- PeakDensityParam(maxFeatures = 33, sampleGroups = 4)
    checkEquals(maxFeatures(p), 33)
    checkException(PeakDensityParam(maxFeatures = -4, sampleGroups = 3))
}

test_MzClustParam <- function() {
    ## Check getter/setter methods:
    p <- new("MzClustParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    checkEquals(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    checkEquals(sampleGroups(p), 1:4)
    p <- MzClustParam(sampleGroups = c("a", "a", "b"))
    checkEquals(sampleGroups(p), c("a", "a", "b"))

    p <- new("MzClustParam", ppm = 3)
    checkEquals(ppm(p), 3)
    ppm(p) <- 20
    checkEquals(ppm(p), 20)
    p <- MzClustParam(ppm = 33)
    checkEquals(ppm(p), 33)
    checkException(MzClustParam(ppm = -4))

    p <- new("MzClustParam", absMz = 3)
    checkEquals(absMz(p), 3)
    absMz(p) <- 20
    checkEquals(absMz(p), 20)
    p <- MzClustParam(absMz = 33)
    checkEquals(absMz(p), 33)
    checkException(MzClustParam(absMz = -4))
    
    ## minFraction
    p <- new("MzClustParam", minFraction = 0.7)
    checkEquals(minFraction(p), 0.7)
    minFraction(p) <- 0.2
    checkEquals(minFraction(p), 0.2)
    p <- MzClustParam(minFraction = 0.4)
    checkEquals(minFraction(p), 0.4)
    checkException(MzClustParam(minFraction = -4))
    checkException(minFraction(p) <- c(0.3, 0.2, 0.4))
    checkException(minFraction(p) <- 2)

    ## minSamples
    p <- new("MzClustParam", minSamples = 3)
    checkEquals(minSamples(p), 3)
    minSamples(p) <- 20
    checkEquals(minSamples(p), 20)
    p <- MzClustParam(minSamples = 33)
    checkEquals(minSamples(p), 33)
    checkException(MzClustParam(minSamples = -4))
}

test_NearestPeaksParam <- function() {
    ## Check getter/setter methods:
    p <- new("NearestPeaksParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    checkEquals(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    checkEquals(sampleGroups(p), 1:4)
    p <- NearestPeaksParam(sampleGroups = c("a", "a", "b"))
    checkEquals(sampleGroups(p), c("a", "a", "b"))

    p <- new("NearestPeaksParam", mzVsRtBalance = 3)
    checkEquals(mzVsRtBalance(p), 3)
    mzVsRtBalance(p) <- 20
    checkEquals(mzVsRtBalance(p), 20)
    p <- NearestPeaksParam(mzVsRtBalance = 33)
    checkEquals(mzVsRtBalance(p), 33)
    checkException(NearestPeaksParam(mzVsRtBalance = -4))
    checkException(NearestPeaksParam(mzVsRtBalance = 1:4))

    p <- new("NearestPeaksParam", absMz = 3)
    checkEquals(absMz(p), 3)
    absMz(p) <- 20
    checkEquals(absMz(p), 20)
    p <- NearestPeaksParam(absMz = 33)
    checkEquals(absMz(p), 33)
    checkException(NearestPeaksParam(absMz = -4))
    checkException(NearestPeaksParam(absMz = 1:3))

    p <- new("NearestPeaksParam", absRt = 3)
    checkEquals(absRt(p), 3)
    absRt(p) <- 20
    checkEquals(absRt(p), 20)
    p <- NearestPeaksParam(absRt = 33)
    checkEquals(absRt(p), 33)
    checkException(NearestPeaksParam(absRt = -4))
    checkException(NearestPeaksParam(absRt = 1:3))
    
    p <- new("NearestPeaksParam", kNN = 3)
    checkEquals(kNN(p), 3)
    kNN(p) <- 20
    checkEquals(kNN(p), 20)
    p <- NearestPeaksParam(kNN = 33)
    checkEquals(kNN(p), 33)
    checkException(NearestPeaksParam(kNN = -4))
    checkException(NearestPeaksParam(kNN = 1:3))
}

test_PeakGroupsParam <- function() {
    ## Check getter/setter methods:
    p <- new("PeakGroupsParam", minFraction = 0.8)
    checkEquals(minFraction(p), 0.8)
    minFraction(p) <- 0.3
    checkEquals(minFraction(p), 0.3)
    p <- PeakGroupsParam(minFraction = 0.7)
    checkEquals(minFraction(p), 0.7)
    checkException(minFraction(p) <- c(2, 2))
    checkException(minFraction(p) <- -1)
    checkException(minFraction(p) <- 3)
    
    p <- new("PeakGroupsParam", extraPeaks = 2)
    checkEquals(extraPeaks(p), 2)
    extraPeaks(p) <- 0.3
    checkEquals(extraPeaks(p), 0.3)
    p <- PeakGroupsParam(extraPeaks = 7)
    checkEquals(extraPeaks(p), 7)
    checkException(extraPeaks(p) <- c(2, 2))
    checkException(extraPeaks(p) <- -1)

    p <- new("PeakGroupsParam", span = 0.5)
    checkEquals(span(p), 0.5)
    span(p) <- 0.3
    checkEquals(span(p), 0.3)
    p <- PeakGroupsParam(span = 7)
    checkEquals(span(p), 7)
    checkException(span(p) <- c(2, 2))
    checkException(span(p) <- -1)

    p <- new("PeakGroupsParam", smooth = "linear")
    checkEquals(smooth(p), "linear")
    smooth(p) <- "loess"
    checkEquals(smooth(p), "loess")
    p <- PeakGroupsParam(smooth = "linear")
    checkEquals(smooth(p), "linear")
    checkException(smooth(p) <- "other")
    checkException(smooth(p) <- c("linear", "loess"))

    p <- new("PeakGroupsParam", family = "symmetric")
    checkEquals(family(p), "symmetric")
    family(p) <- "gaussian"
    checkEquals(family(p), "gaussian")
    p <- PeakGroupsParam(family = "symmetric")
    checkEquals(family(p), "symmetric")
    checkException(family(p) <- "other")
    checkException(family(p) <- c("symmetric", "gaussian"))

    mt <- matrix(1:4, 1:4)
    p <- new("PeakGroupsParam", peakGroupsMatrix = mt)
    checkEquals(peakGroupsMatrix(p), mt)
    peakGroupsMatrix(p) <- mt + 2
    checkEquals(peakGroupsMatrix(p), mt + 2)
    p <- PeakGroupsParam(peakGroupsMatrix = mt)
    checkEquals(peakGroupsMatrix(p), mt)
}


test_ObiwarpParam <- function() {
    library(xcms)
    library(RUnit)
    ## Check getter/setter methods:
    p <- new("ObiwarpParam", binSize = 0.8)
    checkEquals(binSize(p), 0.8)
    binSize(p) <- 0.3
    checkEquals(binSize(p), 0.3)
    p <- ObiwarpParam(binSize = 0.7)
    checkEquals(binSize(p), 0.7)
    checkException(binSize(p) <- c(2, 2))
    checkException(binSize(p) <- -1)
    
    p <- new("ObiwarpParam", centerSample = 2L)
    checkEquals(centerSample(p), 2L)
    centerSample(p) <- 1
    checkEquals(centerSample(p), 1L)
    p <- ObiwarpParam(centerSample = 7)
    checkEquals(centerSample(p), 7)
    checkException(centerSample(p) <- c(2, 2))
    checkException(centerSample(p) <- -1)

    p <- new("ObiwarpParam", response = 3L)
    checkEquals(response(p), 3L)
    response(p) <- 5
    checkEquals(response(p), 5L)
    p <- ObiwarpParam(response = 7)
    checkEquals(response(p), 7)
    checkException(response(p) <- c(2, 2))
    checkException(response(p) <- -1)
    checkException(response(p) <- 200)

    p <- new("ObiwarpParam", distFun = "euc")
    checkEquals(distFun(p), "euc")
    checkEquals(gapInit(p), 0.9)
    checkEquals(gapExtend(p), 1.8)
    distFun(p) <- "cor"
    checkEquals(distFun(p), "cor")
    checkEquals(gapInit(p), 0.3)
    checkEquals(gapExtend(p), 2.4)
    distFun(p) <- "cov"
    checkEquals(distFun(p), "cov")
    checkEquals(gapInit(p), 0)
    checkEquals(gapExtend(p), 11.7)
    distFun(p) <- "prd"
    checkEquals(distFun(p), "prd")
    checkEquals(gapInit(p), 0)
    checkEquals(gapExtend(p), 7.8)
    p <- ObiwarpParam(distFun = "cov")
    checkEquals(distFun(p), "cov")
    checkException(distFun(p) <- c("a", "cov"))
    checkException(distFun(p) <- "other")

    p <- new("ObiwarpParam", gapInit = 4.2)
    checkEquals(gapInit(p), 4.2)
    gapInit(p) <- 5.2
    checkEquals(gapInit(p), 5.2)
    p <- ObiwarpParam(gapInit = 3.1)
    checkEquals(gapInit(p), 3.1)
    checkException(gapInit(p) <- c(2, 2))
    checkException(gapInit(p) <- -1)

    p <- new("ObiwarpParam", gapExtend = 4.2)
    checkEquals(gapExtend(p), 4.2)
    gapExtend(p) <- 5.2
    checkEquals(gapExtend(p), 5.2)
    p <- ObiwarpParam(gapExtend = 3.1)
    checkEquals(gapExtend(p), 3.1)
    checkException(gapExtend(p) <- c(2, 2))
    checkException(gapExtend(p) <- -1)

    p <- new("ObiwarpParam", factorDiag = 4.2)
    checkEquals(factorDiag(p), 4.2)
    factorDiag(p) <- 1.2
    checkEquals(factorDiag(p), 1.2)
    p <- ObiwarpParam(factorDiag = 3.1)
    checkEquals(factorDiag(p), 3.1)
    checkException(factorDiag(p) <- c(2, 2))
    checkException(factorDiag(p) <- -1)

    p <- new("ObiwarpParam", factorGap = 4.2)
    checkEquals(factorGap(p), 4.2)
    factorGap(p) <- 4.2
    checkEquals(factorGap(p), 4.2)
    p <- ObiwarpParam(factorGap = 3.1)
    checkEquals(factorGap(p), 3.1)
    checkException(factorGap(p) <- c(2, 2))
    checkException(factorGap(p) <- -1)

    p <- new("ObiwarpParam", localAlignment = TRUE)
    checkEquals(localAlignment(p), TRUE)
    localAlignment(p) <- FALSE
    checkEquals(localAlignment(p), FALSE)
    p <- ObiwarpParam(localAlignment = TRUE)
    checkEquals(localAlignment(p), TRUE)
    checkException(localAlignment(p) <- c(TRUE, FALSE))

    p <- new("ObiwarpParam", initPenalty = 4.2)
    checkEquals(initPenalty(p), 4.2)
    initPenalty(p) <- 2.2
    checkEquals(initPenalty(p), 2.2)
    p <- ObiwarpParam(initPenalty = 3.1)
    checkEquals(initPenalty(p), 3.1)
    checkException(factorGap(p) <- c(2, 2))
    checkException(factorGap(p) <- -1)
}

test_GenericParam <- function() {
    prm <- GenericParam(fun = "mean")
    checkEquals(prm@fun, "mean")
    ## Errors
    checkException(GenericParam(args = list(na.rm = TRUE)))
    checkException(GenericParam(fun = c("a", "b")))
}

test_FillChromPeaksParam <- function() {
    ## Check getter/setter methods:
    p <- new("FillChromPeaksParam", expandMz = 0.8)
    checkEquals(expandMz(p), 0.8)
    expandMz(p) <- 0.3
    checkEquals(expandMz(p), 0.3)
    p <- FillChromPeaksParam(expandMz = 0.7)
    checkEquals(expandMz(p), 0.7)
    checkException(expandMz(p) <- c(2, 2))
    checkException(expandMz(p) <- -2)

    p <- new("FillChromPeaksParam", expandRt = 0.8)
    checkEquals(expandRt(p), 0.8)
    expandRt(p) <- 0.3
    checkEquals(expandRt(p), 0.3)
    p <- FillChromPeaksParam(expandRt = 0.7)
    checkEquals(expandRt(p), 0.7)
    checkException(expandRt(p) <- c(2, 2))
    checkException(expandRt(p) <- -2)

    p <- new("FillChromPeaksParam", ppm = 8)
    checkEquals(ppm(p), 8)
    ppm(p) <- 3
    checkEquals(ppm(p), 3)
    p <- FillChromPeaksParam(ppm = 7)
    checkEquals(ppm(p), 7)
    checkException(ppm(p) <- c(2, 2))
    checkException(ppm(p) <- -2)
}


test_CalibrantMassParam <- function() {

    p <- new("CalibrantMassParam")
    checkTrue(validObject(p))
    p@method <- "other"
    checkException(validObject(p))
    mzs <- rnorm(200, mean = 500)
    p <- new("CalibrantMassParam")
    p@mz <- list(mzs)
    checkException(validObject(p))
    
    ## Constructor.
    p <- xcms:::CalibrantMassParam(mz = mzs, mzabs = 3, mzppm = 9,
                                   neighbors = 4, method = "shift")
    checkEquals(xcms:::.mz(p)[[1]], sort(mzs))
    checkEquals(xcms:::.mzabs(p), 3)
    checkEquals(xcms:::.mzppm(p), 9)
    checkEquals(xcms:::.neighbors(p), 4L)
    checkEquals(xcms:::.method(p), "shift")
    
}
