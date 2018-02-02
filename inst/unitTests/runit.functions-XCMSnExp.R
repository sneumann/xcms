## Unit tests for functions in functions-XCMSnExp.R

test_plotChromPeakDensity <- function() {
    mzr <- c(305.05, 305.15)
    plotChromPeakDensity(xod_x, mz = mzr)

    ## Use the full range.
    plotChromPeakDensity(xod_x)
    
    plotChromPeakDensity(xod_x, mz = c(0, 1))
    plotChromPeakDensity(xod_x, mz = c(300, 310), pch = 16, xlim = c(2500, 4000))
}

test_plotChromPeaks <- function() {
    ## Plot the full range.
    plotChromPeaks(xod_x)

    ## mz range
    plotChromPeaks(xod_x, ylim = c(453, 455))
    plotChromPeaks(xod_x, ylim = c(453.2, 453.201), xlim = c(2500, 3500))

    ## mzr <- c(453.2, 453.201)
    ## chrs <- chromatogram(xod_x, mz = mzr, rt = c(2500, 3500))
    ## plot(chrs)
    ## highlightChromPeaks(xod_x, mz = mzr)
}

test_plotChromPeakImage <- function() {
    plotChromPeakImage(xod_x, binSize = 30, log = FALSE)
    ## Check that it works if no peaks were found in one sample.
    tmp <- xod_x
    pks <- chromPeaks(tmp)
    pks <- pks[pks[, "sample"] != 1, ]
    chromPeaks(tmp) <- pks
    plotChromPeakImage(tmp, binSize = 30, log = FALSE)
    plotChromPeakImage(tmp, binSize = 20, log = FALSE)
    plotChromPeakImage(tmp, binSize = 10, log = FALSE, col = topo.colors(64))
}

test_applyAdjustedRtime <- function() {
    checkException(applyAdjustedRtime(faahko_od))
    checkEquals(applyAdjustedRtime(faahko_xod), faahko_xod)
    ## Now really replacing the stuff.
    tmp <- applyAdjustedRtime(xod_r)
    checkTrue(!hasAdjustedRtime(tmp))
    checkEquals(rtime(tmp), adjustedRtime(xod_r))
    checkTrue(length(processHistory(tmp)) == 1)
    tmp <- applyAdjustedRtime(xod_xgrg)
    checkTrue(!hasAdjustedRtime(tmp))
    checkEquals(rtime(tmp), adjustedRtime(xod_xgr))
    checkEquals(tmp, dropAdjustedRtime(tmp))
}

test_concatenate_OnDiskMSnExp <- function() {
    od1 <- readMSData(faahko_3_files[1], mode = "onDisk")
    od2 <- readMSData(faahko_3_files[2:3], mode = "onDisk")
    res <- xcms:::.concatenate_OnDiskMSnExp(od1, od2)
    checkEquals(fileNames(faahko_od), fileNames(res))
    checkEquals(experimentData(faahko_od), experimentData(res))
    checkEquals(fData(faahko_od), fData(res))
    checkEquals(pData(faahko_od), pData(res))
    ## checkEquals(spectra(faahko_od), spectra(res))
    ## Checking data manipulations
    od1 <- smooth(od1)
    od1 <- pickPeaks(od1)
    od2 <- smooth(od2)
    suppressWarnings(res <- xcms:::.concatenate_OnDiskMSnExp(od1, od2))
    checkTrue(length(res@spectraProcessingQueue) == 0)
    checkEquals(fileNames(faahko_od), fileNames(res))
    checkEquals(experimentData(faahko_od), experimentData(res))
    checkEquals(pData(faahko_od), pData(res))

    od2 <- pickPeaks(od2)
    res <- xcms:::.concatenate_OnDiskMSnExp(od1, od2)
    checkTrue(length(res@spectraProcessingQueue) == 2)
    checkEquals(fileNames(faahko_od), fileNames(res))
    checkEquals(experimentData(faahko_od), experimentData(res))
    checkEquals(pData(faahko_od), pData(res))

    checkEquals(xcms:::.concatenate_OnDiskMSnExp(od1), od1)
}

test_concatenate_XCMSnExp <- function() {
    od1 <- readMSData(faahko_3_files[1], mode = "onDisk")
    od2 <- readMSData(faahko_3_files[2], mode = "onDisk")
    od3 <- readMSData(faahko_3_files[3], mode = "onDisk")
    xod1 <- findChromPeaks(od1, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    xod2 <- findChromPeaks(od2, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    xod3 <- findChromPeaks(od3, param = CentWaveParam(noise = 10000,
                                                      snthresh = 40))
    res <- xcms:::.concatenate_XCMSnExp(xod1, xod2, xod3)
    checkEquals(pData(res), pData(faahko_xod))
    checkEquals(fData(res), fData(faahko_xod))
    checkEquals(chromPeaks(res), chromPeaks(faahko_xod))

    res2 <- c(xod1, xod2, xod3)
    checkEquals(pData(res), pData(res2))
    checkEquals(fData(res), fData(res2))
    checkEquals(chromPeaks(res), chromPeaks(res2))
}
