test_that("CentWaveParam works", {
    p <- new("CentWaveParam")
    checkDefaultValues <- function(x) {
        expect_equal(x@ppm, 25)
        expect_equal(x@peakwidth, c(20, 50))
        expect_equal(x@snthresh, 10)
        expect_equal(x@prefilter, c(3, 100))
        expect_equal(x@mzCenterFun, "wMean")
        expect_equal(x@integrate, 1)
        expect_equal(x@mzdiff, -0.001)
        expect_equal(x@fitgauss, FALSE)
        expect_equal(x@noise, 0)
        expect_equal(x@verboseColumns, FALSE)
        expect_equal(x@roiList, list())
        expect_equal(x@firstBaselineCheck, TRUE)
        expect_equal(x@roiScales, numeric())
    }
    checkDefaultValues(p)
    p <- new("CentWaveParam", fitgauss = TRUE)
    expect_equal(p@fitgauss, TRUE)
    p@fitgauss <- FALSE
    checkDefaultValues(p)
    p <- CentWaveParam()
    checkDefaultValues(p)
    p <- CentWaveParam(verboseColumns = TRUE, ppm = 30)
    expect_equal(p@verboseColumns, TRUE)
    expect_equal(p@ppm, 30)
    expect_equal(p@snthresh, 10)
    ## Check errors and validity
    expect_error(CentWaveParam(ppm = -1))
    expect_error(new("CentWaveParam", prefilter = c(3, 4, 5)))

    p@mzCenterFun <- "other"
    expect_error(validObject(p))
    ##plist <- as.list(p)

    ## Check getter/setter methods:
    p <- new("CentWaveParam")
    ppm(p) <- 1
    expect_equal(ppm(p), 1)
    expect_error(ppm(p) <- "a")
    peakwidth(p) <- c(5, 100)
    expect_equal(peakwidth(p), c(5, 100))
    expect_error(peakwidth(p) <- c(1, 2, 3))
    snthresh(p) <- 300
    expect_equal(snthresh(p), 300)
    prefilter(p) <- c(1, 1000)
    expect_equal(prefilter(p), c(1, 1000))
    expect_error(prefilter(p) <- 3)
    mzCenterFun(p) <- "meanApex3"
    expect_equal(mzCenterFun(p), "meanApex3")
    expect_error(mzCenterFun(p) <- "other")
    integrate(p) <- 2
    expect_equal(integrate(p), 2L)
    expect_error(integrate(p) <- 3L)
    mzdiff(p) <- 123.3
    expect_equal(mzdiff(p), 123.3)
    expect_error(mzdiff(p) <- 1:4)
    fitgauss(p) <- TRUE
    expect_equal(fitgauss(p), TRUE)
    expect_error(fitgauss(p) <- c(TRUE, TRUE))
    noise(p) <- 5
    expect_equal(noise(p), 5)
    expect_error(noise(p) <- 1:4)
    verboseColumns(p) <- TRUE
    expect_equal(verboseColumns(p), TRUE)
    expect_error(verboseColumns(p) <- c(TRUE, TRUE))
    firstBaselineCheck(p) <- TRUE
    expect_equal(firstBaselineCheck(p), TRUE)
    expect_error(firstBaselineCheck(p) <- c(TRUE, TRUE))
})

test_that("MatchedFilterParam works", {
    p <- new("MatchedFilterParam")
    checkDefaultValues <- function(x) {
        expect_equal(x@binSize, 0.1)
        expect_equal(x@impute, "none")
        expect_equal(x@baseValue, numeric())
        expect_equal(x@distance, numeric())
        expect_equal(x@fwhm, 30)
        expect_equal(x@sigma, 12.73994)
        expect_equal(x@max, 5)
        expect_equal(x@snthresh, 10)
        expect_equal(x@steps, 2)
        expect_equal(x@mzdiff, 0.6)
        expect_equal(x@index, FALSE)
    }
    checkDefaultValues(p)
    p <- new("MatchedFilterParam", binSize = 0.5)
    expect_equal(binSize(p), 0.5)
    ## binSize
    binSize(p) <- 0.6
    expect_equal(binSize(p), 0.6)
    p <- MatchedFilterParam(binSize = 0.6)
    expect_equal(binSize(p), 0.6)
    expect_error(binSize(p) <- 1:5)
    ## impute
    impute(p) <- "lin"
    expect_equal(impute(p), "lin")
    p <- MatchedFilterParam(impute = "lin")
    expect_equal(impute(p), "lin")
    expect_error(impute(p) <- "other")
    ## baseValue
    baseValue(p) <- 3
    expect_equal(baseValue(p), 3)
    p <- MatchedFilterParam(baseValue = 3)
    expect_equal(baseValue(p), 3)
    expect_error(baseValue(p) <- 1:3)
    ## distance
    distance(p) <- 4
    expect_equal(distance(p), 4)
    p <- MatchedFilterParam(distance = 4)
    expect_equal(distance(p), 4)
    expect_error(distance(p) <- 1:3)
    ## fwhm
    fwhm(p) <- 11
    expect_equal(fwhm(p), 11)
    p <- MatchedFilterParam(fwhm = 11)
    expect_equal(fwhm(p), 11)
    expect_error(fwhm(p) <- 1:3)
    ## sigma
    sigma(p) <- 12
    expect_equal(sigma(p), 12)
    p <- MatchedFilterParam(sigma = 12)
    expect_equal(sigma(p), 12)
    expect_error(sigma(p) <- 1:3)
    ## max
    max(p) <- 13
    expect_equal(max(p), 13)
    p <- MatchedFilterParam(max = 13)
    expect_equal(max(p), 13)
    expect_error(max(p) <- 1:3)
    ## snthresh
    snthresh(p) <- 14
    expect_equal(snthresh(p), 14)
    p <- MatchedFilterParam(snthresh = 14)
    expect_equal(snthresh(p), 14)
    expect_error(snthresh(p) <- 1:3)
    ## steps
    steps(p) <- 15
    expect_equal(steps(p), 15)
    p <- MatchedFilterParam(steps = 15)
    expect_equal(steps(p), 15)
    expect_error(steps(p) <- 1:3)
    ## mzdiff
    mzdiff(p) <- 16
    expect_equal(mzdiff(p), 16)
    p <- MatchedFilterParam(mzdiff = 16)
    expect_equal(mzdiff(p), 16)
    expect_error(mzdiff(p) <- 1:3)
    ## index
    index(p) <- TRUE
    expect_equal(index(p), TRUE)
    p <- MatchedFilterParam(index = TRUE)
    expect_equal(index(p), TRUE)
    expect_error(index(p) <- 1:3)
})

test_that("MassifquantParam works", {    
    ## Check getter/setter methods:
    p <- new("MassifquantParam")
    ppm(p) <- 1
    expect_equal(ppm(p), 1)
    p <- MassifquantParam(ppm = 1)
    expect_equal(ppm(p), 1)
    expect_error(ppm(p) <- 1:4)
    peakwidth(p) <- c(5, 100)
    expect_equal(peakwidth(p), c(5, 100))
    p <- MassifquantParam(peakwidth = c(5, 100))
    expect_equal(peakwidth(p), c(5, 100))
    expect_error(peakwidth(p) <- c(1, 2, 3))
    snthresh(p) <- 300
    expect_equal(snthresh(p), 300)
    p <- MassifquantParam(snthresh = 300)
    expect_equal(snthresh(p), 300)
    expect_error(snthresh(p) <- 2:4)
    prefilter(p) <- c(1, 1000)
    expect_equal(prefilter(p), c(1, 1000))
    p <- MassifquantParam(prefilter = c(1, 1000))
    expect_error(prefilter(p) <- 3)
    mzCenterFun(p) <- "meanApex3"
    expect_equal(mzCenterFun(p), "meanApex3")
    p <- MassifquantParam(mzCenterFun = "meanApex3")
    expect_equal(mzCenterFun(p), "meanApex3")
    expect_error(mzCenterFun(p) <- "other")
    integrate(p) <- 2
    expect_equal(integrate(p), 2L)
    p <- MassifquantParam(integrate = 2)
    expect_equal(integrate(p), 2L)
    expect_error(integrate(p) <- 3L)
    mzdiff(p) <- 123.3
    expect_equal(mzdiff(p), 123.3)
    p <- MassifquantParam(mzdiff = 123.3)
    expect_equal(mzdiff(p), 123.3)
    expect_error(mzdiff(p) <- 1:4)
    fitgauss(p) <- TRUE
    expect_equal(fitgauss(p), TRUE)
    p <- MassifquantParam(fitgauss = TRUE)
    expect_equal(fitgauss(p), TRUE)
    expect_error(fitgauss(p) <- c(TRUE, TRUE))
    noise(p) <- 5
    expect_equal(noise(p), 5)
    p <- MassifquantParam(noise = 5)
    expect_equal(noise(p), 5)
    expect_error(noise(p) <- 1:4)
    verboseColumns(p) <- TRUE
    expect_equal(verboseColumns(p), TRUE)
    p <- MassifquantParam(verboseColumns = TRUE)
    expect_equal(verboseColumns(p), TRUE)
    expect_error(verboseColumns(p) <- c(TRUE, TRUE))
    criticalValue(p) <- 34
    expect_equal(criticalValue(p), 34)
    p <- MassifquantParam(criticalValue = 34)
    expect_equal(criticalValue(p), 34)
    expect_error(criticalValue(p) <- as.numeric(3:5))
    consecMissedLimit(p) <- 4
    expect_equal(consecMissedLimit(p), 4L)
    p <- MassifquantParam(consecMissedLimit = 4)
    expect_equal(consecMissedLimit(p), 4L)
    expect_error(consecMissedLimit(p) <- 3:5)
    ## unions, checkBack, withWave
    unions(p) <- 0
    expect_equal(unions(p), 0L)
    p <- MassifquantParam(unions = 1)
    expect_equal(unions(p), 1L)
    expect_error(unions(p) <- 4)
    checkBack(p) <- 0
    expect_equal(checkBack(p), 0L)
    p <- MassifquantParam(checkBack = 1)
    expect_equal(checkBack(p), 1L)
    expect_error(checkBack(p) <- 4)
    withWave(p) <- TRUE
    expect_equal(withWave(p), TRUE)
    p <- MassifquantParam(withWave = TRUE)
    expect_equal(withWave(p), TRUE)
    expect_error(withWave(p) <- c(TRUE, FALSE))
})

test_that("MSWParam works", {
    ## Check getter/setter methods:
    p <- new("MSWParam", snthresh = 14)
    expect_equal(snthresh(p), 14)
    snthresh(p) <- 13
    expect_equal(snthresh(p), 13)
    p <- MSWParam(snthresh = 21)
    expect_equal(snthresh(p), 21)
    expect_error(MSWParam(snthresh = -4))

    p <- new("MSWParam", verboseColumns = TRUE)
    expect_equal(verboseColumns(p), TRUE)
    verboseColumns(p) <- FALSE
    expect_equal(verboseColumns(p), FALSE)
    p <- MSWParam(verboseColumns = TRUE)
    expect_equal(verboseColumns(p), TRUE)
    expect_error(MSWParam(verboseColumns = c(FALSE, TRUE)))

    p <- new("MSWParam", scales = 1:5)
    expect_equal(scales(p), 1:5)
    scales(p) <- 6:12
    expect_equal(scales(p), 6:12)
    p <- MSWParam(scales = 1:3)
    expect_equal(scales(p), 1:3)

    p <- new("MSWParam", nearbyPeak = FALSE)
    expect_equal(nearbyPeak(p), FALSE)
    nearbyPeak(p) <- TRUE
    expect_equal(nearbyPeak(p), TRUE)
    p <- MSWParam(nearbyPeak = FALSE)
    expect_equal(nearbyPeak(p), FALSE)
    expect_error(nearbyPeak(p) <- c(TRUE, FALSE))

    p <- new("MSWParam", peakScaleRange = 8)
    expect_equal(peakScaleRange(p), 8)
    peakScaleRange(p) <- 7
    expect_equal(peakScaleRange(p), 7)
    p <- MSWParam(peakScaleRange = 9)
    expect_equal(peakScaleRange(p), 9)
    expect_error(peakScaleRange(p) <- 4:6)

    p <- new("MSWParam", ampTh = 0.4)
    expect_equal(ampTh(p), 0.4)
    ampTh(p) <- 7
    expect_equal(ampTh(p), 7)
    p <- MSWParam(ampTh = 9)
    expect_equal(ampTh(p), 9)
    expect_error(ampTh(p) <- 4:6)

    p <- new("MSWParam", minNoiseLevel = 0.4)
    expect_equal(minNoiseLevel(p), 0.4)
    minNoiseLevel(p) <- 7
    expect_equal(minNoiseLevel(p), 7)
    p <- MSWParam(minNoiseLevel = 9)
    expect_equal(minNoiseLevel(p), 9)
    expect_error(minNoiseLevel(p) <- -4)

    p <- new("MSWParam", ridgeLength = 14)
    expect_equal(ridgeLength(p), 14)
    ridgeLength(p) <- 7
    expect_equal(ridgeLength(p), 7)
    p <- MSWParam(ridgeLength = 9)
    expect_equal(ridgeLength(p), 9)
    expect_error(ridgeLength(p) <- -4)

    p <- new("MSWParam", peakThr = 14)
    expect_equal(peakThr(p), 14)
    peakThr(p) <- 7
    expect_equal(peakThr(p), 7)
    p <- MSWParam(peakThr = 9)
    expect_equal(peakThr(p), 9)
    expect_error(peakThr(p) <- 3:5)

    p <- new("MSWParam", tuneIn = FALSE)
    expect_equal(tuneIn(p), FALSE)
    tuneIn(p) <- TRUE
    expect_equal(tuneIn(p), TRUE)
    p <- MSWParam(tuneIn = FALSE)
    expect_equal(tuneIn(p), FALSE)
    expect_error(tuneIn(p) <- c(TRUE, TRUE))

    p <- new("MSWParam", addParams = list(a = 3, b = 4))
    expect_equal(addParams(p), list(a = 3, b = 4))
    addParams(p) <- list(d = "a", e = 2:4)
    expect_equal(addParams(p), list(d = "a", e = 2:4))
    p <- MSWParam(otherp = 1:4, z = "a")
    expect_equal(addParams(p), list(otherp = 1:4, z = "a"))

    ## Check the .param2list method:
    p <- new("MSWParam", addParams = list(z = "z", bla = 1:4))
    L <- xcms:::.param2list(p)
    expect_equal(L$z, "z")
    expect_equal(L$bla, 1:4)
    expect_equal(L$snthresh, 3)
    p <- new("MSWParam")
    L <- xcms:::.param2list(p)
    expect_true(!any(names(L) == "addParams"))
    expect_equal(L$snthresh, 3)
})

test_that("CentWavePredIsoParam works", {
    ## Check getter/setter methods:
    p <- new("CentWavePredIsoParam", ppm = 14)
    expect_equal(ppm(p), 14)
    ppm(p) <- 13
    expect_equal(ppm(p), 13)
    p <- CentWavePredIsoParam(ppm = 21)
    expect_equal(ppm(p), 21)
    expect_error(CentWavePredIsoParam(ppm = -4))

    p <- new("CentWavePredIsoParam", peakwidth = c(1, 2))
    expect_equal(peakwidth(p), c(1, 2))
    peakwidth(p) <- c(2, 3)
    expect_equal(peakwidth(p), c(2, 3))
    p <- CentWavePredIsoParam(peakwidth = c(3, 4))
    expect_equal(peakwidth(p), c(3, 4))
    expect_error(CentWavePredIsoParam(peakwidth = 1:3))

    p <- new("CentWavePredIsoParam", snthresh = 5)
    expect_equal(snthresh(p), 5)
    snthresh(p) <- 6
    expect_equal(snthresh(p), 6)
    p <- CentWavePredIsoParam(snthresh = 7)
    expect_equal(snthresh(p), 7)
    expect_error(CentWavePredIsoParam(snthresh = 1:3))

    p <- new("CentWavePredIsoParam", prefilter = c(1, 2))
    expect_equal(prefilter(p), c(1, 2))
    prefilter(p) <- c(2, 3)
    expect_equal(prefilter(p), c(2, 3))
    p <- CentWavePredIsoParam(prefilter = c(3, 4))
    expect_equal(prefilter(p), c(3, 4))
    expect_error(CentWavePredIsoParam(prefilter = 1:3))

    p <- new("CentWavePredIsoParam", mzCenterFun = "meanApex3")
    expect_equal(mzCenterFun(p), "meanApex3")
    mzCenterFun(p) <- "mean"
    expect_equal(mzCenterFun(p), "mean")
    p <- CentWavePredIsoParam(mzCenterFun = "mean")
    expect_equal(mzCenterFun(p), "mean")
    expect_error(CentWavePredIsoParam(mzCenterFun = "median"))

    p <- new("CentWavePredIsoParam", integrate = 2L)
    expect_equal(integrate(p), 2L)
    integrate(p) <- 1L
    expect_equal(integrate(p), 1L)
    p <- CentWavePredIsoParam(integrate = 2L)
    expect_equal(integrate(p), 2L)
    expect_error(CentWavePredIsoParam(integrate = 3L))

    p <- new("CentWavePredIsoParam", mzdiff = 1)
    expect_equal(mzdiff(p), 1)
    mzdiff(p) <- 2
    expect_equal(mzdiff(p), 2)
    p <- CentWavePredIsoParam(mzdiff = 3)
    expect_equal(mzdiff(p), 3)
    expect_error(CentWavePredIsoParam(mzdiff = 2:3))

    p <- new("CentWavePredIsoParam", fitgauss = TRUE)
    expect_equal(fitgauss(p), TRUE)
    fitgauss(p) <- FALSE
    expect_equal(fitgauss(p), FALSE)
    p <- CentWavePredIsoParam(fitgauss = TRUE)
    expect_equal(fitgauss(p), TRUE)
    expect_error(CentWavePredIsoParam(fitgauss = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", noise = 3)
    expect_equal(noise(p), 3)
    noise(p) <- 6
    expect_equal(noise(p), 6)
    p <- CentWavePredIsoParam(noise = 8)
    expect_equal(noise(p), 8)
    expect_error(CentWavePredIsoParam(noise = c(3, 5)))

    p <- new("CentWavePredIsoParam", verboseColumns = TRUE)
    expect_equal(verboseColumns(p), TRUE)
    verboseColumns(p) <- FALSE
    expect_equal(verboseColumns(p), FALSE)
    p <- CentWavePredIsoParam(verboseColumns = TRUE)
    expect_equal(verboseColumns(p), TRUE)
    expect_error(CentWavePredIsoParam(verboseColumns = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", firstBaselineCheck = TRUE)
    expect_equal(firstBaselineCheck(p), TRUE)
    firstBaselineCheck(p) <- FALSE
    expect_equal(firstBaselineCheck(p), FALSE)
    p <- CentWavePredIsoParam(firstBaselineCheck = FALSE)
    expect_equal(firstBaselineCheck(p), FALSE)
    expect_error(CentWavePredIsoParam(firstBaselineCheck = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", snthreshIsoROIs = 5)
    expect_equal(snthreshIsoROIs(p), 5)
    snthreshIsoROIs(p) <- 6
    expect_equal(snthreshIsoROIs(p), 6)
    p <- CentWavePredIsoParam(snthreshIsoROIs = 7)
    expect_equal(snthreshIsoROIs(p), 7)
    expect_error(CentWavePredIsoParam(snthreshIsoROIs = 1:3))

    p <- new("CentWavePredIsoParam", maxCharge = 5L)
    expect_equal(maxCharge(p), 5L)
    maxCharge(p) <- 6
    expect_equal(maxCharge(p), 6L)
    p <- CentWavePredIsoParam(maxCharge = 7)
    expect_equal(maxCharge(p), 7L)
    expect_error(CentWavePredIsoParam(maxCharge = 1:3))

    p <- new("CentWavePredIsoParam", maxIso = 5L)
    expect_equal(maxIso(p), 5L)
    maxIso(p) <- 6
    expect_equal(maxIso(p), 6L)
    p <- CentWavePredIsoParam(maxIso = 7)
    expect_equal(maxIso(p), 7L)
    expect_error(CentWavePredIsoParam(maxIso = 1:3))

    p <- new("CentWavePredIsoParam", mzIntervalExtension = FALSE)
    expect_equal(mzIntervalExtension(p), FALSE)
    mzIntervalExtension(p) <- TRUE
    expect_equal(mzIntervalExtension(p), TRUE)
    p <- CentWavePredIsoParam(mzIntervalExtension = FALSE)
    expect_equal(mzIntervalExtension(p), FALSE)
    expect_error(CentWavePredIsoParam(mzIntervalExtension = c(TRUE, FALSE)))

    p <- new("CentWavePredIsoParam", polarity = "positive")
    expect_equal(polarity(p), "positive")
    polarity(p) <- "negative"
    expect_equal(polarity(p), "negative")
    p <- CentWavePredIsoParam(polarity = "negative")
    expect_equal(polarity(p), "negative")
    expect_error(CentWavePredIsoParam(polarity = "bla"))

    ## Check the .param2list method:
    p <- new("CentWavePredIsoParam", snthresh = 123)
    L <- xcms:::.param2list(p)
    expect_equal(L$snthresh, 123)
})

test_that("PeakDensityParam works", {    
    ## Check getter/setter methods:
    p <- new("PeakDensityParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    expect_equal(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    expect_equal(sampleGroups(p), 1:4)
    p <- PeakDensityParam(sampleGroups = c("a", "a", "b"))
    expect_equal(sampleGroups(p), c("a", "a", "b"))
    expect_error(sampleGroups(p) <- NULL)
    expect_error(sampleGroups(p) <- c(2, 2, NA))
    expect_error(PeakDensityParam())
    expect_error(PeakDensityParam(sampleGroups = c(1, 1, NA)))
    
    p <- new("PeakDensityParam", bw = 3)
    expect_equal(bw(p), 3)
    bw(p) <- 20
    expect_equal(bw(p), 20)
    p <- PeakDensityParam(bw = 33, sampleGroups = "a")
    expect_equal(bw(p), 33)
    expect_error(PeakDensityParam(bw = -4, sampleGroups = 1))

    ## minFraction
    p <- new("PeakDensityParam", minFraction = 0.7)
    expect_equal(minFraction(p), 0.7)
    minFraction(p) <- 0.2
    expect_equal(minFraction(p), 0.2)
    p <- PeakDensityParam(minFraction = 0.4, sampleGroups = 3)
    expect_equal(minFraction(p), 0.4)
    expect_error(PeakDensityParam(minFraction = -4, sampleGroups = 3))
    expect_error(minFraction(p) <- c(0.3, 0.2, 0.4))
    expect_error(minFraction(p) <- 2)

    ## minSamples
    p <- new("PeakDensityParam", minSamples = 3)
    expect_equal(minSamples(p), 3)
    minSamples(p) <- 20
    expect_equal(minSamples(p), 20)
    p <- PeakDensityParam(minSamples = 33, sampleGroups = rep(3, 100))
    expect_equal(minSamples(p), 33)
    expect_error(PeakDensityParam(minSamples = -4, sampleGroups = 4))

    ## binSize
    p <- new("PeakDensityParam", binSize = 3)
    expect_equal(binSize(p), 3)
    binSize(p) <- 20
    expect_equal(binSize(p), 20)
    p <- PeakDensityParam(binSize = 0.3, sampleGroups = 4)
    expect_equal(binSize(p), 0.3)
    expect_error(PeakDensityParam(binSize = -4, sampleGroups = 4))

    ## maxFeatures
    p <- new("PeakDensityParam", maxFeatures = 3)
    expect_equal(maxFeatures(p), 3)
    maxFeatures(p) <- 20
    expect_equal(maxFeatures(p), 20)
    p <- PeakDensityParam(maxFeatures = 33, sampleGroups = 4)
    expect_equal(maxFeatures(p), 33)
    expect_error(PeakDensityParam(maxFeatures = -4, sampleGroups = 3))
})

test_that("MzClustParam works", {
    ## Check getter/setter methods:
    p <- new("MzClustParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    expect_equal(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    expect_equal(sampleGroups(p), 1:4)
    p <- MzClustParam(sampleGroups = c("a", "a", "b"))
    expect_equal(sampleGroups(p), c("a", "a", "b"))

    p <- new("MzClustParam", ppm = 3)
    expect_equal(ppm(p), 3)
    ppm(p) <- 20
    expect_equal(ppm(p), 20)
    p <- MzClustParam(ppm = 33)
    expect_equal(ppm(p), 33)
    expect_error(MzClustParam(ppm = -4))

    p <- new("MzClustParam", absMz = 3)
    expect_equal(absMz(p), 3)
    absMz(p) <- 20
    expect_equal(absMz(p), 20)
    p <- MzClustParam(absMz = 33)
    expect_equal(absMz(p), 33)
    expect_error(MzClustParam(absMz = -4))
    
    ## minFraction
    p <- new("MzClustParam", minFraction = 0.7)
    expect_equal(minFraction(p), 0.7)
    minFraction(p) <- 0.2
    expect_equal(minFraction(p), 0.2)
    p <- MzClustParam(minFraction = 0.4)
    expect_equal(minFraction(p), 0.4)
    expect_error(MzClustParam(minFraction = -4))
    expect_error(minFraction(p) <- c(0.3, 0.2, 0.4))
    expect_error(minFraction(p) <- 2)

    ## minSamples
    p <- new("MzClustParam", minSamples = 3)
    expect_equal(minSamples(p), 3)
    minSamples(p) <- 20
    expect_equal(minSamples(p), 20)
    p <- MzClustParam(minSamples = 33)
    expect_equal(minSamples(p), 33)
    expect_error(MzClustParam(minSamples = -4))
})

test_that("NearestPeaksParam works", {
    ## Check getter/setter methods:
    p <- new("NearestPeaksParam", sampleGroups = c(1, 1, 1, 2, 2, 3, 4))
    expect_equal(sampleGroups(p), c(1, 1, 1, 2, 2, 3, 4))
    sampleGroups(p) <- 1:4
    expect_equal(sampleGroups(p), 1:4)
    p <- NearestPeaksParam(sampleGroups = c("a", "a", "b"))
    expect_equal(sampleGroups(p), c("a", "a", "b"))

    p <- new("NearestPeaksParam", mzVsRtBalance = 3)
    expect_equal(mzVsRtBalance(p), 3)
    mzVsRtBalance(p) <- 20
    expect_equal(mzVsRtBalance(p), 20)
    p <- NearestPeaksParam(mzVsRtBalance = 33)
    expect_equal(mzVsRtBalance(p), 33)
    expect_error(NearestPeaksParam(mzVsRtBalance = -4))
    expect_error(NearestPeaksParam(mzVsRtBalance = 1:4))

    p <- new("NearestPeaksParam", absMz = 3)
    expect_equal(absMz(p), 3)
    absMz(p) <- 20
    expect_equal(absMz(p), 20)
    p <- NearestPeaksParam(absMz = 33)
    expect_equal(absMz(p), 33)
    expect_error(NearestPeaksParam(absMz = -4))
    expect_error(NearestPeaksParam(absMz = 1:3))

    p <- new("NearestPeaksParam", absRt = 3)
    expect_equal(absRt(p), 3)
    absRt(p) <- 20
    expect_equal(absRt(p), 20)
    p <- NearestPeaksParam(absRt = 33)
    expect_equal(absRt(p), 33)
    expect_error(NearestPeaksParam(absRt = -4))
    expect_error(NearestPeaksParam(absRt = 1:3))
    
    p <- new("NearestPeaksParam", kNN = 3)
    expect_equal(kNN(p), 3)
    kNN(p) <- 20
    expect_equal(kNN(p), 20)
    p <- NearestPeaksParam(kNN = 33)
    expect_equal(kNN(p), 33)
    expect_error(NearestPeaksParam(kNN = -4))
    expect_error(NearestPeaksParam(kNN = 1:3))
})

test_that("PeakGroupsParam works", {
    ## Check getter/setter methods:
    p <- new("PeakGroupsParam", minFraction = 0.8)
    expect_equal(minFraction(p), 0.8)
    minFraction(p) <- 0.3
    expect_equal(minFraction(p), 0.3)
    p <- PeakGroupsParam(minFraction = 0.7)
    expect_equal(minFraction(p), 0.7)
    expect_error(minFraction(p) <- c(2, 2))
    expect_error(minFraction(p) <- -1)
    expect_error(minFraction(p) <- 3)
    
    p <- new("PeakGroupsParam", extraPeaks = 2)
    expect_equal(extraPeaks(p), 2)
    extraPeaks(p) <- 0.3
    expect_equal(extraPeaks(p), 0.3)
    p <- PeakGroupsParam(extraPeaks = 7)
    expect_equal(extraPeaks(p), 7)
    expect_error(extraPeaks(p) <- c(2, 2))
    expect_error(extraPeaks(p) <- -1)

    p <- new("PeakGroupsParam", span = 0.5)
    expect_equal(span(p), 0.5)
    span(p) <- 0.3
    expect_equal(span(p), 0.3)
    p <- PeakGroupsParam(span = 7)
    expect_equal(span(p), 7)
    expect_error(span(p) <- c(2, 2))
    expect_error(span(p) <- -1)

    p <- new("PeakGroupsParam", smooth = "linear")
    expect_equal(smooth(p), "linear")
    smooth(p) <- "loess"
    expect_equal(smooth(p), "loess")
    p <- PeakGroupsParam(smooth = "linear")
    expect_equal(smooth(p), "linear")
    expect_error(smooth(p) <- "other")
    expect_error(smooth(p) <- c("linear", "loess"))

    p <- new("PeakGroupsParam", family = "symmetric")
    expect_equal(family(p), "symmetric")
    family(p) <- "gaussian"
    expect_equal(family(p), "gaussian")
    p <- PeakGroupsParam(family = "symmetric")
    expect_equal(family(p), "symmetric")
    expect_error(family(p) <- "other")
    expect_error(family(p) <- c("symmetric", "gaussian"))

    mt <- matrix(1:4, 1:4)
    p <- new("PeakGroupsParam", peakGroupsMatrix = mt)
    expect_equal(peakGroupsMatrix(p), mt)
    peakGroupsMatrix(p) <- mt + 2
    expect_equal(peakGroupsMatrix(p), mt + 2)
    p <- PeakGroupsParam(peakGroupsMatrix = mt)
    expect_equal(peakGroupsMatrix(p), mt)
})

test_that("ObiwarpParam works", {
    ## Check getter/setter methods:
    p <- new("ObiwarpParam", binSize = 0.8)
    expect_equal(binSize(p), 0.8)
    binSize(p) <- 0.3
    expect_equal(binSize(p), 0.3)
    p <- ObiwarpParam(binSize = 0.7)
    expect_equal(binSize(p), 0.7)
    expect_error(binSize(p) <- c(2, 2))
    expect_error(binSize(p) <- -1)
    
    p <- new("ObiwarpParam", centerSample = 2L)
    expect_equal(centerSample(p), 2L)
    centerSample(p) <- 1
    expect_equal(centerSample(p), 1L)
    p <- ObiwarpParam(centerSample = 7)
    expect_equal(centerSample(p), 7)
    expect_error(centerSample(p) <- c(2, 2))
    expect_error(centerSample(p) <- -1)

    p <- new("ObiwarpParam", response = 3L)
    expect_equal(response(p), 3L)
    response(p) <- 5
    expect_equal(response(p), 5L)
    p <- ObiwarpParam(response = 7)
    expect_equal(response(p), 7)
    expect_error(response(p) <- c(2, 2))
    expect_error(response(p) <- -1)
    expect_error(response(p) <- 200)

    p <- new("ObiwarpParam", distFun = "euc")
    expect_equal(distFun(p), "euc")
    expect_equal(gapInit(p), 0.9)
    expect_equal(gapExtend(p), 1.8)
    distFun(p) <- "cor"
    expect_equal(distFun(p), "cor")
    expect_equal(gapInit(p), 0.3)
    expect_equal(gapExtend(p), 2.4)
    distFun(p) <- "cov"
    expect_equal(distFun(p), "cov")
    expect_equal(gapInit(p), 0)
    expect_equal(gapExtend(p), 11.7)
    distFun(p) <- "prd"
    expect_equal(distFun(p), "prd")
    expect_equal(gapInit(p), 0)
    expect_equal(gapExtend(p), 7.8)
    p <- ObiwarpParam(distFun = "cov")
    expect_equal(distFun(p), "cov")
    expect_error(distFun(p) <- c("a", "cov"))
    expect_error(distFun(p) <- "other")

    p <- new("ObiwarpParam", gapInit = 4.2)
    expect_equal(gapInit(p), 4.2)
    gapInit(p) <- 5.2
    expect_equal(gapInit(p), 5.2)
    p <- ObiwarpParam(gapInit = 3.1)
    expect_equal(gapInit(p), 3.1)
    expect_error(gapInit(p) <- c(2, 2))
    expect_error(gapInit(p) <- -1)

    p <- new("ObiwarpParam", gapExtend = 4.2)
    expect_equal(gapExtend(p), 4.2)
    gapExtend(p) <- 5.2
    expect_equal(gapExtend(p), 5.2)
    p <- ObiwarpParam(gapExtend = 3.1)
    expect_equal(gapExtend(p), 3.1)
    expect_error(gapExtend(p) <- c(2, 2))
    expect_error(gapExtend(p) <- -1)

    p <- new("ObiwarpParam", factorDiag = 4.2)
    expect_equal(factorDiag(p), 4.2)
    factorDiag(p) <- 1.2
    expect_equal(factorDiag(p), 1.2)
    p <- ObiwarpParam(factorDiag = 3.1)
    expect_equal(factorDiag(p), 3.1)
    expect_error(factorDiag(p) <- c(2, 2))
    expect_error(factorDiag(p) <- -1)

    p <- new("ObiwarpParam", factorGap = 4.2)
    expect_equal(factorGap(p), 4.2)
    factorGap(p) <- 4.2
    expect_equal(factorGap(p), 4.2)
    p <- ObiwarpParam(factorGap = 3.1)
    expect_equal(factorGap(p), 3.1)
    expect_error(factorGap(p) <- c(2, 2))
    expect_error(factorGap(p) <- -1)

    p <- new("ObiwarpParam", localAlignment = TRUE)
    expect_equal(localAlignment(p), TRUE)
    localAlignment(p) <- FALSE
    expect_equal(localAlignment(p), FALSE)
    p <- ObiwarpParam(localAlignment = TRUE)
    expect_equal(localAlignment(p), TRUE)
    expect_error(localAlignment(p) <- c(TRUE, FALSE))

    p <- new("ObiwarpParam", initPenalty = 4.2)
    expect_equal(initPenalty(p), 4.2)
    initPenalty(p) <- 2.2
    expect_equal(initPenalty(p), 2.2)
    p <- ObiwarpParam(initPenalty = 3.1)
    expect_equal(initPenalty(p), 3.1)
    expect_error(factorGap(p) <- c(2, 2))
    expect_error(factorGap(p) <- -1)
})

test_that("GenericParam works", {
    prm <- GenericParam(fun = "mean")
    expect_equal(prm@fun, "mean")
    ## Errors
    expect_error(GenericParam(args = list(na.rm = TRUE)))
    expect_error(GenericParam(fun = c("a", "b")))
})

test_that("FillChromPeaksParam works", {
    ## Check getter/setter methods:
    p <- new("FillChromPeaksParam", expandMz = 0.8)
    expect_equal(expandMz(p), 0.8)
    expandMz(p) <- 0.3
    expect_equal(expandMz(p), 0.3)
    p <- FillChromPeaksParam(expandMz = 0.7)
    expect_equal(expandMz(p), 0.7)
    expect_error(expandMz(p) <- c(2, 2))
    expect_error(expandMz(p) <- -2)

    p <- new("FillChromPeaksParam", expandRt = 0.8)
    expect_equal(expandRt(p), 0.8)
    expandRt(p) <- 0.3
    expect_equal(expandRt(p), 0.3)
    p <- FillChromPeaksParam(expandRt = 0.7)
    expect_equal(expandRt(p), 0.7)
    expect_error(expandRt(p) <- c(2, 2))
    expect_error(expandRt(p) <- -2)

    p <- new("FillChromPeaksParam", ppm = 8)
    expect_equal(ppm(p), 8)
    ppm(p) <- 3
    expect_equal(ppm(p), 3)
    p <- FillChromPeaksParam(ppm = 7)
    expect_equal(ppm(p), 7)
    expect_error(ppm(p) <- c(2, 2))
    expect_error(ppm(p) <- -2)
})

test_that("CalibrantMassParam works", {
    p <- new("CalibrantMassParam")
    expect_true(validObject(p))
    p@method <- "other"
    expect_error(validObject(p))
    mzs <- rnorm(200, mean = 500)
    p <- new("CalibrantMassParam")
    p@mz <- list(mzs)
    expect_error(validObject(p))
    
    ## Constructor.
    p <- CalibrantMassParam(mz = mzs, mzabs = 3, mzppm = 9,
                            neighbors = 4, method = "shift")
    expect_equal(.mz(p)[[1]], sort(mzs))
    expect_equal(.mzabs(p), 3)
    expect_equal(.mzppm(p), 9)
    expect_equal(.neighbors(p), 4L)
    expect_equal(.method(p), "shift")
})
