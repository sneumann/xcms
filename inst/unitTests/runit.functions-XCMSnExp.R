## Unit tests for functions in functions-XCMSnExp.R

test_plotChromPeakDensity <- function() {
    mzr <- c(305.05, 305.15)
    plotChromPeakDensity(xod_x, mz = mzr)

    ## Use the full range.
    plotChromPeakDensity(xod_x)
    
    plotChromPeakDensity(xod_x, mz = c(0, 1))
    plotChromPeakDensity(xod_x, mz = c(300, 310), pch = 16, xlim = c(2500, 4000))
}
