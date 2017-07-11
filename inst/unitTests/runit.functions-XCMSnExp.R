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
