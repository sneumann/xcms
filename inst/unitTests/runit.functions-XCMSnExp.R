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
