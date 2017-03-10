## Test functions to evaluate binning of values.
## library(faahKO)
## fs <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"))
## xr <- xcmsRaw(fs, profstep = 0)

############################################################
## Test profBin
##
## profBin matches with the implementation in R.
## Does however not match with binYonX, step = 0.2!
## profBin is now deprecated!
dontrun_test_profBin <- function() {
    xr <- deepCopy(faahko_xr_1)
    X <- xr@env$mz
    Y <- xr@env$intensity
    scanidx <- xr@scanindex
    ## Get the data from the first spectrum:
    X <- X[1:scanidx[2]]
    Y <- Y[1:scanidx[2]]
    step <- 0.1
    ## Calculate the number of bins we need (that's taken from findPeaks.matchedFilter):
    mass <- seq(floor(min(xr@env$mz)/step)*step,
                ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBin(X, Y, nBins, xstart = min(mass),
                           xend = max(mass))
    ## Do the same with the reference implementation
    resR <- xcms:::profBinR(X, Y, nBins = nBins, fromX = min(mass), toX = max(mass))
    checkEquals(resX, resR)

    ## The same with step 0.2:
    step <- 0.2
    mass <- seq(floor(min(xr@env$mz)/step)*step,
                ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBin(X, Y, nBins, xstart = min(mass),
                           xend = max(mass))
    resR <- xcms:::profBinR(X, Y, nBins = nBins, fromX = min(mass),
                            toX = max(mass))
    checkEquals(resX, resR)

    ## Using an odd step:
    step <- 0.13
    mass <- seq(floor(min(xr@env$mz)/step)*step, ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBin(X, Y, nBins, xstart = min(mass),
                           xend = max(mass))
    resR <- xcms:::profBinR(X, Y, nBins = nBins, fromX = min(mass),
                            toX = max(mass))
    checkEquals(resX, resR)

    ## Test iterative binning with profBin
    X <- 1:30
    Y <- 1:30
    step <- 0.2
    bufsize <- 10
    mass <- seq(floor(min(X)/step) * step, ceiling(max(X)/step) * step,
                by = step)
    bufidx <- integer(length(mass))
    idxrange <- c(1, bufsize)
    bufidx[idxrange[1]:idxrange[2]] <- 1:bufsize
    lookahead <- 1
    lookbehind <- 1
    res1 <- xcms:::profBin(X, Y, bufsize, mass[1], mass[bufsize])
    res1_dx <- xcms:::profBin_test(X, Y, bufsize, mass[1], mass[bufsize])
    i <- 10
    bufidx[i + lookahead] == 0
    ## re-read buffer.
    bufidx[idxrange[1]:idxrange[2]] <- 0
    idxrange <- c(max(1, i - lookbehind), min(bufsize + i -1 - lookbehind,
                                              length(mass)))
    bufidx[idxrange[1]:idxrange[2]] <- 1:(diff(idxrange) + 1)
    res2 <- xcms:::profBin(X, Y, diff(idxrange) + 1, mass[idxrange[1]],
                           mass[idxrange[2]])
    res2_dx <- xcms:::profBin_test(X, Y, diff(idxrange) + 1, mass[idxrange[1]],
                                   mass[idxrange[2]])
}


############################################################
## Test profBinLin
fails_test_profBinLin <- function() {

    library(xcms)
    library(RUnit)

    X <- 1:11
    ## If there are no bins with missing values we expect the
    ## results to be identical.
    res <- xcms:::profBinLin(X, X, 5L)
    ## THAT FAILS!!!
    ## checkEquals(res, xcms:::profBin(X, X, 5L))
    brks <- breaks_on_nBins(1, 11, 5L, TRUE)

    ## Test with different Ys:
    Y <- c(2, 1, 3, 4, 5, 6, 7, 8, 9, 10, 11)
    res <- xcms:::profBinLin(X, Y, 5L, 1, 11)
    Y <- c(2, 1, 6, 7, 8, 9, 10, 11)
    X <- c(1, 2, 6, 7, 8, 9, 10, 11)
    xcms:::profBinLin(X, Y, 5L, 1, 11)
    xcms:::profBinLinBase(X, Y, 5L, 1, 11)
    xcms:::profBin(X, Y, 5)
    ## ??????
    Y <- sort(Y)
    xcms:::profBin(X, Y, 5)
    xcms:::profBinLin(X, Y, 5L, 1, 11)
    xcms:::profBinLinBase(X, Y, 5L, 1, 11)

    ## That's to understand what binLin is doing:
    ## Empty values further up.
    Y <- c(1, 2, 3, 4, 8, 9, 10, 11)
    X <- c(1, 2, 3, 4, 8, 9, 10, 11)
    xcms:::profBin(X, Y, 5)
    xcms:::profBinLin(X, Y, 5)
    ## The 3rd bin gets populated with (4+9)/2 = 6.5
    ## Is it using the actual values or the bin values?
    Y <- c(1, 2, 3, 3, 8, 11, 10, 11)
    ## If it uses the closest value to the bin it should be (3+8) / 2 = 5.5
    ## if it's using the binned values it should be (3+11) / 2 = 7
    xcms:::profBinLin(X, Y, 5)
    ## = 7: so it's using the max bin value!
    ## To be sure:
    Y <- c(1, 2, 4, 3, 10, 20, 10, 11)
    X <- c(1, 2, 3, 4, 8, 9, 10, 11)
    xcms:::profBinLin(X, Y, 5)
    ## The 3rd bin value is 12, thus (4+20) / 2.

    ## Next example: 2 empty bins.
    Y <- c(1, 2, 3, 10, 11)
    X <- c(1, 2, 3, 10, 11)
    xcms:::profBinLin(X, Y, 5)
    ## 1, 3, 5.666, 8.333, 11
    ## it's x + (10 - 3) / 3


    ## Reporting this as issue: (issue #46)
    X <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    ## Binsize if we want to have 5 bins: (10-1)/ (5-1) = 2.25
    ## bins will thus be: (starting from 1, centered at 1)
    ## [-0.125,2.125), [2.125,4.375), [4.375,6.625), [6.625,8.875), [8.875,11.125)
    ## Binning values from above and selecting the max:
    ## bin1: 1, 2 -> 2
    ## bin2: 3, 4 -> 4
    ## bin3: 5, 6 -> 6
    ## bin4: 7, 8 -> 8
    ## bin5: 9, 10 -> 10
    ## profBin: bins as described above, select max value within bin
    xcms:::profBin(X, X, 5L)  ## OK
    ## profBinLin: same as profBin, but linearly interpolate values into which nothing was binned;
    ## all bins have values, thus we expect it to be the same
    xcms:::profBinLin(X, X, 5L)
    ## The value for the first bin is always wrong.
    ## profBinLinBase: same as profBinLin
    xcms:::profBinLinBase(X, X, 5L)  ## OK

    X <- sort(abs(rnorm(5000, mean = 500, sd = 300)))
    a <- xcms:::profBinLinBase(X, X, 200)
    b <- xcms:::profBinLin(X, X, 200)
    head(a)
    head(b)
    head(binYonX(X, X, nBins = 200, shiftByHalfBinSize = TRUE, method = "min")$y)

    ## We've got 1 missing bin.
    X <- c(1, 2, 7, 8, 9, 10)
    xcms:::profBinLin(X, X, 5L)
    xcms:::profBinLinBase(X, X, 5L)

}

############################################################
## Test profBinLinBase
## profBinLinBase is now deprecated!
dontrun_test_profBinLinBase <- function() {
    ## According to the help:
    ## o baselevel: half of minimal signal in the data.
    ## o basespace: 0.075. Linear interpolation is used for data points falling
    ##   within 2*basespace of each other. -> m/z length after which the
    ##   signal drops to baselevel.

    X <- 1:16
    Y <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
    nas <- is.na(Y)
    brks <- breaks_on_nBins(1, 16, nBins = 16, shiftByHalfBinSize = TRUE)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16)
    checkEquals(resX, binYonX(X, Y, nBins = 16, baseValue = 0.5,
                              shiftByHalfBinSize = TRUE)$y)
    ## Modify basespace such that 1 neighboring bin is considered.
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, param = list(basespace = 1))
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5 + 1.5/2, 2)
    checkEquals(resX, expect)
    ## Modify basespace such that 2 neighboring bins are considered.
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, param = list(basespace = 2))
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 1/4, 4 - 2/4, 4 - 3/4, 3,
                3 - 1/5, 3 - 2/5, 3 - 3/5, 3 - 4/5, 2)
    checkEquals(resX, expect)
    ## What happens if there are all NAs up to the end?
    Y[16] <- NA
    nas <- is.na(Y)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, xstart = 1, xend = 16,
                                  param = list(basespace = 1))
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(resX, expect)
    ## And if we start with NAs?
    Y[1] <- NA
    nas <- is.na(Y)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, xstart = 1, xend = 16,
                                  param = list(basespace = 1))
    expect <- c(0.5, 0.75, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(resX, expect)
    ## Only the first is an NA:
    Y[2] <- 2
    nas <- is.na(Y)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, xstart = 1, xend = 16,
                                  param = list(basespace = 1))
    expect <- c(0.5 + 1.5/2, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(resX, expect)

    ## Stretch of NAs from the beginning.
    Y[1:3] <- NA
    nas <- is.na(Y)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, xstart = 1, xend = 16,
                                  param = list(basespace = 1, baselevel = 0.5))
    expect <- c(0.5, 0.5, 0.5 + 1.5/2, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(resX, expect)

    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16, xstart = 1, xend = 16,
                                  param = list(basespace = 2, baselevel = 0.5))
    ## The same but considering two neighboring elements.
    expect <- c(0.5, 0.5 + 1.5/3, 0.5 + 2* 1.5/3, 2, 2+2/3, 2+2*2/3, 4, 4 - 1/4,
                4 - 2/4, 4 - 3/4, 3, 3 - 2.5/3, 3 - 2*2.5/3, 0.5, 0.5, 0.5)
    checkEquals(resX, expect)
}

############################################################
## Test NA handling in binYonX
test_binning_NA_handling <- function() {

    library(xcms)
    library(RUnit)

    ## NA values in y should be ignored internally.
    X <- 1:10
    Y <- X
    Y[4] <- NA

    res <- binYonX(X, Y, nBins = 5L)
    checkEquals(res$y, c(2, 3, 6, 8, 10))
    Y[3] <- NA
    res <- binYonX(X, Y, nBins = 5L, baseValue = 0)
    checkEquals(res$y, c(2, 0, 6, 8, 10))

    ## Initialize with NA
    res <- binYonX(X, Y, nBins = 5L, baseValue = NA)
    checkEquals(res$y, c(2, NA, 6, 8, 10))

    X <- 1:10
    res <- binYonX(X, X, binSize = 0.5, baseValue = 0)
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0,
                         8, 0, 9, 10))
    res <- binYonX(X, X, binSize = 0.5)
    checkEquals(res$y, c(1, NA, 2, NA, 3, NA, 4, NA, 5, NA, 6, NA, 7,
                         NA, 8, NA, 9, 10))
}


############################################################
## binYonX with method max.
test_binning_max <- function() {

    library(xcms)
    library(RUnit)

    X <- 1:10
    Y <- 1:10

    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, Y, nBins = 5L)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(breaks_on_nBins(min(X), max(X), 5L)))
    res <- binYonX(X, Y, nBins = 5L, returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2.05, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(breaks_on_nBins(1, 10, 5L)))
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10,
                   returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  Define fromIdx, toIdx; binning will only take place in that sub-set.
    res <- binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = 10)
    brks <- breaks_on_nBins(min(X[1]), max(X[10]), nBins = 5L)
    checkEquals(res$y, c(2.05, 4, 6, 8, 10))
    checkEquals(res$x, breakMidPoint(brks))
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  If we provide fromX and toX the result will be different, as the breaks
    ##  are calculated differently.
    ##  This calculates 5 bins on range(X)
    res <- binYonX(X, X, nBins = 5L, binFromX = min(X), binToX = max(X),
                   fromIdx = 1, toIdx = 10)
    brks <- breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkEquals(res$y, c(3, 5, 8, 10, NA)) ## because toIdx = 10 the last is 0
    checkEquals(res$x, breakMidPoint(brks))
    res <- binYonX(X, nBins = 5L, binFromX = min(X), binToX = max(X),
                   fromIdx = 1, toIdx = 10, returnIndex = TRUE)
    checkEquals(res$index, c(3, 5, 8, 10, NA))
    ##  Check the fromIdx and toIdx again:
    X <- rep(1:10, 3)
    Y <- 1:30
    brks <- breaks_on_nBins(1, 10, nBins = 5L)
    ##  In this case we have to set sortedX = TRUE
    res <- binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                   fromIdx = 11, toIdx = 20, sortedX = TRUE)
    checkEquals(res$y, c(12, 14, 16, 18, 20))
    checkEquals(res$x, breakMidPoint(brks))
    res <- binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                   fromIdx = 11, toIdx = 20, sortedX = TRUE,
                   returnIndex = TRUE)
    checkEquals(res$index, c(12, 14, 16, 18, 20))
    ##  re-order X
    X <- sample(X, size = length(X))
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(breaks_on_nBins(min(X), max(X), 5L)))

    ## Exceptions:
    X <- 1:10
    checkException(binYonX(X, X, nBins = -4))
    checkException(binYonX(X, X, nBins = 0))
    checkException(binYonX(X, X, nBins = 4, fromIdx = 0))
    checkException(binYonX(X, X, nBins = 4, toIdx = 0))
    checkException(binYonX(X, X, nBins = 4, fromIdx = 7, toIdx = 3))
    checkException(binYonX(X, X, nBins = 4, toIdx = 30))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, X, binSize = 0.5, baseValue = 0)
    brks <- breaks_on_binSize(1, 10, binSize = 0.5)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkIdentical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
                            0, 9, 10))
    checkIdentical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45)
    brks <- breaks_on_binSize(1, 10, binSize = 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.9, 3.9-5.35, 5.35-6.8, 6.8-8.25, 8.25-10
    checkIdentical(res$y, c(2, 3, 5, 6, 8, 10))
    checkIdentical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45, binFromX = 3, binToX = 9)
    brks <- breaks_on_binSize(3, 9, 1.45)
    checkIdentical(res$y, c(4, 5, 7, 9))
    checkIdentical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45, fromIdx = 3, toIdx = 10,
                   binFromX = min(X), binToX = max(X), baseValue = 0)
    brks <- breaks_on_binSize(1, 10, 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.90, 3.90-5.35, 5.35-6.80, 6.80-8.25, 8.25-9.70, 9.70-10.00
    checkIdentical(res$y, c(0, 3, 5, 6, 8, 10))
    checkIdentical(res$x, breakMidPoint(brks))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks)
    checkEquals(res$y, binYonX(X, X, nBins = 4L)$y)
    checkEquals(res$x, breakMidPoint(brks))
    ##  same breaks but sub-setting x
    res <- binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                   baseValue = 0)
    checkEquals(res$y, c(0, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, X, breaks = brks)
    checkEquals(res$y, c(4, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))

    ## Test on real data:
    xr <- deepCopy(faahko_xr_1)
    X <- xr@env$mz
    Y <- xr@env$intensity
    xRangeFull <- range(X)
    scanidx <- xr@scanindex
    ## Get the data from the first spectrum:
    X1 <- X[1:scanidx[2]]
    Y1 <- Y[1:scanidx[2]]

    ## ####
    ## Define the number of bins.
    step <- 0.1
    shift <- TRUE
    mass <- seq(floor(min(xRangeFull)/step)*step,
                ceiling(max(xRangeFull)/step)*step, by = step)
    nBins <- length(mass)
    resR <- xcms:::profBinR(X1, Y1, nBins = nBins, fromX = min(xRangeFull),
                            toX = max(xRangeFull), shiftByHalfBinSize = shift)
    res <- binYonX(X1, Y1, nBins = nBins, binFromX = min(xRangeFull),
                   binToX = max(xRangeFull), shiftByHalfBinSize = shift,
                   baseValue = 0)
    checkEquals(res$y, resR)

    ## Next
    step <- 0.2
    shift <- TRUE
    mass <- seq(floor(min(xRangeFull)/step)*step,
                ceiling(max(xRangeFull)/step)*step, by = step)
    nBins <- length(mass)
    resR <- xcms:::profBinR(X1, Y1, nBins = nBins, fromX = min(xRangeFull),
                            toX = max(xRangeFull), shiftByHalfBinSize = shift)
    res <- binYonX(X1, Y1, nBins = nBins, binFromX = min(xRangeFull),
                   binToX = max(xRangeFull), shiftByHalfBinSize = shift,
                   baseValue = 0)
    checkEquals(res$y, resR)
    shift <- FALSE
    mass <- seq(floor(min(xRangeFull)/step)*step,
                ceiling(max(xRangeFull)/step)*step, by = step)
    nBins <- length(mass)
    resR <- xcms:::profBinR(X1, Y1, nBins = nBins, fromX = min(xRangeFull),
                            toX = max(xRangeFull), shiftByHalfBinSize = shift)
    res <- binYonX(X1, Y1, nBins = nBins, binFromX = min(xRangeFull),
                   binToX = max(xRangeFull), shiftByHalfBinSize = shift,
                   baseValue = 0)
    checkEquals(res$y, resR)

    step <- 0.13
    shift <- TRUE
    mass <- seq(floor(min(xRangeFull)/step)*step,
                ceiling(max(xRangeFull)/step)*step, by = step)
    nBins <- length(mass)
    resR <- xcms:::profBinR(X1, Y1, nBins = nBins, fromX = min(xRangeFull),
                            toX = max(xRangeFull), shiftByHalfBinSize = shift)
    res <- binYonX(X1, Y1, nBins = nBins, binFromX = min(xRangeFull),
                   binToX = max(xRangeFull), shiftByHalfBinSize = shift,
                   baseValue = 0)
    checkEquals(res$y, resR)
}

## Test binning using min
test_binning_min <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "min")
    checkEquals(res$y, c(1, 3, 5, 7, 9))
    res <- binYonX(X, nBins = 5L, method = "min", returnIndex = TRUE)
    checkEquals(res$index, c(1, 3, 5, 7, 9))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10, method = "min")
    checkEquals(res$y, c(1, 3, 5, 7, 9))
    ##  Define fromIdx, toIdx
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10, method = "min",
                   binFromX = min(X), binToX = max(X),
                   baseValue = 0)
    brks <- breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkEquals(res$y, c(1, 4, 6, 9, 0))
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   method = "min", binFromX = min(X), binToX = max(X),
                   returnIndex = TRUE)
    checkEquals(res$index, c(1, 4, 6, 9, NA))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "min", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    breaks_on_binSize(min(X), max(X), binSize = 0.5)
    checkIdentical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    ##
    res <- binYonX(X, binSize = 1.45, method = "min")
    breaks_on_binSize(min(X), max(X), binSize = 1.45)
    checkIdentical(res$y, c(1, 3, 4, 6, 7, 9))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "min")
    checkEquals(res$y, c(1, 4, 6, 8))
    ##  same breaks but sub-setting x
    res <- binYonX(X, breaks = brks, fromIdx = 4, toIdx = 9, method = "min",
                   baseValue = 0)
    checkIdentical(res$y, c(0, 4, 6, 8))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "min")
    checkIdentical(res$y, c(3, 5, 6, 8))
}

## Test binning using sum
test_binning_sum <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "sum")
    breaks_on_nBins(1, 10, nBins = 5)
    checkEquals(res$y, c(3, 7, 11, 15, 19))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                   method = "sum")
    checkEquals(res$y, c(3.05, 7, 11, 15, 19))
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                   method = "sum", returnIndex = TRUE)
    checkEquals(length(res), 2)
    ##  Define fromIdx, toIdx
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   method = "sum", binFromX = min(X),
                   binToX = max(X), baseValue = 0)
    breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkIdentical(res$y, c(6.05, 9, 21, 19, 0))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "sum", baseValue = 0)
    breaks_on_binSize(min(X), max(X), binSize = 0.5)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkIdentical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    res <- binYonX(X, binSize = 1.45, method = "sum")
    breaks_on_binSize(min(X), max(X), binSize = 1.45)
    checkIdentical(res$y, c(3, 3, 9, 6, 15, 19))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "sum")
    checkIdentical(res$y, c(6, 9, 13, 27))
    ##  same breaks but sub-setting x
    res <- binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                   method = "sum", baseValue = 0)
    checkIdentical(res$y, c(0, 9, 13, 17))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "sum")
    checkIdentical(res$y, c(7, 5, 13, 17))
}

## Test binning using mean
test_binning_mean <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "mean")
    breaks_on_nBins(1, 10, nBins = 5)
    checkEquals(res$y, c(1.5, 3.5, 5.5, 7.5, 9.5))
    Y <- X
    Y[3] <- NA
    res <- binYonX(X, Y, nBins = 5L, method = "mean")
    checkEquals(res$y, c(1.5, 4, 5.5, 7.5, 9.5))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "mean", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "mean")
    checkEquals(res$y, c(2, 4.5, 6.5, 9))
}

test_breaks <- function() {

    library(xcms)
    library(RUnit)

    ## Test generation of breaks for binning.
    ## o nBins
    res <- breaks_on_nBins(1, 10, 4)
    checkEquals(res, seq(1, 10, length.out = 5))
    res <- breaks_on_nBins(2, 8, 20)
    checkEquals(res, seq(2, 8, length.out = 21))
    ##
    brksR <- seq(200, 600, length.out = 2002)
    brks <- breaks_on_nBins(200, 600, nBins = 2001)
    checkEquals(brks, brksR)
    ## Simulate shift by half bin size
    brksR <- seq((200 - 0.1), (600 + 0.1), length.out = 2002)
    brks <- breaks_on_nBins(200, 600, nBins = 2001,
                            shiftByHalfBinSize = TRUE)
    checkEquals(brks, brksR)

    ## o binSize
    res <- breaks_on_binSize(1, 10, 0.13)
    resR <- seq(1, 10, by = 0.13)
    checkEquals(res[-length(res)], resR[-length(resR)])
    checkEquals(res[length(res)], 10)
    ##   Will create one bin more.
    res <- breaks_on_binSize(1, 10, 0.51)
    resR <- seq(1, 10, by = 0.51)
    checkEquals(res[-length(res)], resR)
    checkEquals(res[length(res)], 10)
    ##
    brksR <- seq(200, 600, by = 0.2)
    brks <- breaks_on_binSize(200, 600, binSize = 0.2)
    checkEquals(brks, brksR)
    ## Simulate shift by half bin size
    brksR <- seq((200 - 0.1), (600 + 0.1), by = 0.2)
    brks <- breaks_on_binSize((200 - 0.1), (600 + 0.1),
                              binSize = 0.2)
    checkEquals(brks, brksR)
    ##
    ## Ultimate fix for issue #118 
    brksR <- seq((200 - 0.1), (600), by = 0.2)
    brks <- breaks_on_binSize((200 - 0.1), (600), binSize = 0.2)
    ## Compare them up to the last value, since in R that will be 600-01, while
    ## breaks_on_binSize will ensure that the upper limit (600) is still within
    ## the breaks
    cmn <- 1:(length(brksR) - 1)
    checkEquals(brks[cmn], brksR[cmn])
    ## Now, that below breaks on a windows build machine (issue #127)
    ## checkTrue(length(brks) > length(brksR))
    ## checkEquals(brks[-length(brks)], brksR)
    checkEquals(brks[length(brks)], 600)
}

############################################################
## Test imputation by linear interpolation
test_imputation_lin <- function() {

    library(xcms)
    library(RUnit)

    X <- 1:11
    brks <- breaks_on_nBins(1, 11, 5L)

    Y <- c(1, NA, NA, NA, 5, 6, NA, NA, 9, 10, 11)
    ## No imputation
    resVals <- binYonX(X, Y, nBins = 5L)$y
    ## Expect NA for bins 2 and 4
    checkIdentical(resVals, c(1, NA, 6, NA, 11))
    res <- imputeLinInterpol(resVals, method = "lin")
    checkIdentical(res, c(1, 3.5, 6, 8.5, 11))

    ## Empty bin at the start.
    Y <- c(NA, 4, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    checkIdentical(res, c(2, 4, 6, 8, 11))
    Y <- c(NA, 3, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    checkIdentical(res, c(1.5, 3, 6, 8, 11))

    ## Two empty bins at the start.
    Y <- c(NA, NA, 5, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    checkIdentical(res, c(5/3, 2*5/3, 5, 8, 11))

    ## Empty bin at the end.
    Y <- c(2, 4, 6, 9, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 9, 4.5))

    ## Two empty bins at the end.
    Y <- c(2, 4, 6, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 4, 2))

    ## 3 empty bins at the end.
    Y <- c(2, 4, NA, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 3, 2, 1))

    ## 2 consecutive empty bins.
    Y <- c(2, NA, NA, 8, 10)
    res <- imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 8, 10))

    ## Simulate the profBinLin: no interpolation at both ends:
    Y <- c(2, NA, NA, 8, 10)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    checkEquals(res, c(2, 4, 6, 8, 10))

    ## One NA at start:
    Y <- c(NA, 4, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    checkIdentical(res, c(0, 4, 6, 8, 11))

    ## Two NAs at the start.
    Y <- c(NA, NA, 5, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    checkIdentical(res, c(0, 0, 5, 8, 11))

    ## NA at the end.
    Y <- c(2, 4, 6, 9, NA)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    checkEquals(res, c(2, 4, 6, 9, 0))

    ## Two NAs at the end.
    Y <- c(2, 4, 6, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    checkEquals(res, c(2, 4, 6, 0, 0))
}

############################################################
## Test binning with linbase imputation
test_imputation_linbin <- function() {
    doPlot <- FALSE

    library(xcms)
    library(RUnit)

    ## Construct example:
    ## We're using the same test than we did for profBinLinBase.
    Y <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
    nas <- is.na(Y)
    res <- imputeLinInterpol(Y, method = "linbase",
                             baseValue = 0.5, distance = 0)
    dummy <- Y
    dummy[is.na(dummy)] <- 0.5
    checkEquals(res, dummy)

    ## Modify basespace such that 1 neighboring bin is considered.
    res <- imputeLinInterpol(Y, method = "linbase")
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5 + 1.5/2, 2)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## Modify basespace such that 2 neighboring bins are considered.
    res <- imputeLinInterpol(Y, method = "linbase", distance = 2L)
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 1/4, 4 - 2/4,
                4 - 3/4, 3, 3 - 1/5, 3 - 2/5, 3 - 3/5, 3 - 4/5, 2)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## What happens if there are all NAs up to the end?
    Y[16] <- NA
    res <- imputeLinInterpol(Y, method = "linbase", distance = 1)
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## And if we start with NAs?
    Y[1] <- NA
    res <- imputeLinInterpol(Y, method = "linbase", distance = 1)
    expect <- c(0.5, 0.75, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## Only the first is an NA:
    Y[2] <- 2
    res <- imputeLinInterpol(Y, method = "linbase", distance = 1)
    expect <- c(0.5 + 1.5/2, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## Stretch of NAs from the beginning.
    Y[1:3] <- NA
    res <- imputeLinInterpol(Y, method = "linbase", distance = 1,
                             baseValue = 0.5)
    expect <- c(0.5, 0.5, 0.5 + 1.5/2, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5, 0.5)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## The same but considering two neighboring elements.
    res <- imputeLinInterpol(Y, method = "linbase", distance = 2,
                             baseValue = 0.5)
    expect <- c(0.5, 0.5 + 1.5/3, 0.5 + 2* 1.5/3, 2, 2+2/3, 2+2*2/3, 4, 4 - 1/4,
                4 - 2/4, 4 - 3/4, 3, 3 - 2.5/3, 3 - 2*2.5/3, 0.5, 0.5, 0.5)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }
}


.test_inputs <- function() {

    library(xcms)

    xcms:::testIntegerInput(4)
    xcms:::testIntegerInput(c(4, 5))

    xcms:::testRealInput(4)
    xcms:::testRealInput(c(4.2231444, 5))

}

############################################################
## binYonX_multi
##
## Test binning in subsets.
test_binYonX_multi <- function() {

    library(xcms)
    library(RUnit)

    ## Simple test without actually needing subsets.
    X <- 1:11
    Y <- 1:11
    nBins = 5L

    X <- 1:30
    fromIdx <- c(1, 3, 14)
    toIdx <- c(8, 11, 28)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], binYonX(X, X, nBins = nBins,
                                       fromIdx = fromIdx[i],
                                       toIdx = toIdx[i]))
    }

    ## GC-torture; check if objects in C are save.
    gctorture(on = TRUE)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    gctorture(on = FALSE)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], binYonX(X, X, nBins = nBins,
                                       fromIdx = fromIdx[i],
                                       toIdx = toIdx[i]))
    }

    X <- 1:30
    fromIdx <- c(1, 3, 14)
    toIdx <- c(8, 11, 28)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], binYonX(X, X, nBins = nBins,
                                       fromIdx = fromIdx[i],
                                       toIdx = toIdx[i]))
    }
    ## Returning also the index
    resM <- binYonX(X, nBins = nBins, fromIdx = fromIdx, toIdx = toIdx,
                    returnIndex = TRUE)
    checkTrue(all(unlist(lapply(resM, function(z) {
        return(any(names(z) == "index"))
    })) == TRUE))

    ## With defining fromX and toX before. The breaks will thus
    ## be defined on these values.
    fromX <- 1
    toX <- 30
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx, binFromX = fromX, binToX = toX,
                    baseValue = 17)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], binYonX(X, X, nBins = nBins,
                                       fromIdx = fromIdx[i],
                                       toIdx = toIdx[i],
                                       binFromX = fromX,
                                       binToX = toX,
                                       baseValue = 17))
    }

    ## Error checks.
    checkException(binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = c(2, 3)))
}


############################################################
##
##        COMPARE WITH profBin METHODS
##
############################################################

############################################################
## Compare binYonX, max, with plain profBin
##
## We might drop this test once we drop the profBin method.
dontrun_test_compare_with_profBin <- function() {

    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    a <- xcms:::profBin(X, X, 1000)
    b <- binYonX(X, X, nBins = 1000, sortedX = TRUE,
                 shiftByHalfBinSize = TRUE,
                 baseValue = 0)
    checkIdentical(a, b$y)

    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBin(X, X, 1000),
    ##                binYonX(X, X, nBins = 1000, sortedX = TRUE,
    ##                               shiftByHalfBinSize = TRUE,
    ##                               baseValue = 0))

    ## To be on the save side, compare also with binSize:
    binS <- unique(diff(b$x))[1]
    c <- binYonX(X, binSize = binS, sortedX = TRUE,
                 shiftByHalfBinSize = TRUE,
                 baseValue = 0)
    checkIdentical(length(b$x), length(c$x))
    ## The bins will be the same except for the last one.
    checkEquals(b$x[-length(b$x)], c$x[-length(c$x)])
    checkEquals(b$y, c$y)

    ## NA handling
    X[20:30] <- NA
    checkException(a <- xcms:::profBin(X, X, 1000))
    checkException(b <- binYonX(X, X, nBins = 1000, sortedX = TRUE,
                                shiftByHalfBinSize = TRUE))

    xr <- deepCopy(faahko_xr_1)
    X <- xr@env$mz
    Y <- xr@env$intensity
    scanidx <- xr@scanindex
    ## Get the data from the first spectrum:
    X <- X[1:scanidx[2]]
    Y <- Y[1:scanidx[2]]
    step <- 0.1
    ## Calculate the number of bins we need (that's taken from findPeaks.matchedFilter):
    mass <- seq(floor(min(xr@env$mz)/step)*step, ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBin(X, Y, nBins, xstart = min(mass),
                           xend = max(mass))
    resB <- binYonX(X, Y, nBins = nBins, binFromX = min(mass),
                    binToX = max(mass), sortedX = TRUE,
                    shiftByHalfBinSize = TRUE,
                    baseValue = 0)
    checkIdentical(resX, resB$y)
    ## Do the same with the reference implementation
    resR <- xcms:::profBinR(X, Y, nBins = nBins, fromX = min(mass), toX = max(mass))
    checkEquals(resX, resR)
    ##
    step <- 0.2
    ## Calculate the number of bins we need (that's taken from findPeaks.matchedFilter):
    mass <- seq(floor(min(xr@env$mz)/step)*step, ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBin(X, Y, nBins, xstart = min(mass),
                           xend = max(mass))
    resB <- binYonX(X, Y, nBins = nBins, binFromX = min(mass),
                    binToX = max(mass), sortedX = TRUE,
                    shiftByHalfBinSize = TRUE,
                    baseValue = 0)
    checkIdentical(resX, resB$y)
    ## Do the same with the reference implementation
    resR <- xcms:::profBinR(X, Y, nBins = nBins, fromX = min(mass), toX = max(mass))
    checkEquals(resX, resR)
}


############################################################
## Compare with profBinM
##
## We might drop this test once we deprecate the profBinM method.
dontrun_test_compare_with_profBinM <- function() {

    ## profBinM: takes vectors of data and a second one specifying
    ## non-overlapping sub-sets.
    X <- 1:30
    fromX <- c(1, 11, 19)
    toX <- c(10, 18, 30)
    scIndex <- fromX - 1L
    nBins <- 5L

    resX <- xcms:::profBinM(X, X, scIndex, nBins)
    resM <- binYonX(X, nBins = nBins, binFromX = min(X), binToX = max(X),
                    fromIdx = fromX, toIdx = toX,
                    shiftByHalfBinSize = TRUE, baseValue = 0)
    M <- do.call(cbind, lapply(resM, function(x) {x$y}))
    checkIdentical(resX, M)

    ## A more realistic example.
    set.seed(18011977)
    X <- c(sort(abs(rnorm(3545, mean = 400, sd = 100))),
           sort(abs(rnorm(10000, mean = 500, sd = 200))),
           sort(abs(rnorm(7500, mean = 300, sd = 134))))
    fromX <- c(1, 3546, 13546)
    toX <- c(3545, 13545, 21045)
    nBins <- 2000L
    scIndex <- fromX - 1L
    fX <- min(X)
    tX <- max(X)
    ## Calculate
    resX <- xcms:::profBinM(X, X, scIndex, nBins, xstart = fX, xend = tX)
    resM <- binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                    fromIdx = fromX, toIdx = toX, baseValue = 0,
                    shiftByHalfBinSize = TRUE, sortedX = TRUE)
    M <- do.call(cbind, lapply(resM, function(x) {x$y}))
    checkIdentical(resX, M)

    xr <- deepCopy(faahko_xr_1)
    X <- xr@env$mz
    Y <- xr@env$intensity
    scanidx <- xr@scanindex
    step <- 0.1
    valsPerSpect <- diff(c(scanidx, length(xr@env$mz)))
    mass <- seq(floor(min(xr@env$mz)/step)*step,
                ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    ## profBinM
    resX <- xcms:::profBinM(X, Y, scanidx, nBins, xstart = mass[1],
                            xend = max(mass))
    toX <- cumsum(valsPerSpect)
    fromX <- c(1L, toX[-length(toX)] + 1L)
    ## binYonX
    resB <- binYonX(X, Y, nBins = nBins, binFromX = mass[1],
                    binToX = max(mass), fromIdx = fromX,
                    toIdx = toX, sortedX = TRUE,
                    shiftByHalfBinSize = TRUE,
                    baseValue = 0)
    resB <- do.call(cbind, lapply(resB, function(z) return(z$y)))
    checkIdentical(resB, resX)

    ## step 0.2.
    step <- 0.2
    mass <- seq(floor(min(xr@env$mz)/step)*step,
                ceiling(max(xr@env$mz)/step)*step, by = step)
    nBins <- length(mass)
    resX <- xcms:::profBinM(X, Y, scanidx, nBins, xstart = mass[1],
                            xend = max(mass))
    resB <- binYonX(X, Y, nBins = nBins, binFromX = mass[1],
                    binToX = max(mass), fromIdx = fromX,
                    toIdx = toX, sortedX = TRUE,
                    shiftByHalfBinSize = TRUE,
                    baseValue = 0)
    resB <- do.call(cbind, lapply(resB, function(z) return(z$y)))
    checkIdentical(resB, resX)
    binSize <- (max(mass) - min(mass)) / (nBins - 1)
    brks <- seq(min(mass) - binSize/2, max(mass) + binSize/2, by = binSize)
    resB <- binYonX(X, Y, breaks = brks, fromIdx = fromX,
                    toIdx = toX, sortedX = TRUE,
                    baseValue = 0)
    resB <- do.call(cbind, lapply(resB, function(z) return(z$y)))
    checkIdentical(resB, resX)
}

############################################################
## Compare with profBinLin
##
## We might drop this test once we deprecate the profBinLin method.
dontrun_test_compare_with_profBinLin <- function() {
    ## binLin does linear interpolation of values, in which nothing
    ## was binned... let's see.
    ## binLin does the imputation on the already binned values, this was also
    ## implemented in the impute = "lin" method.

    X <- 1:11
    brks <- breaks_on_nBins(1, 11, 5L, TRUE)

    ## One missing bin in the middle
    Y <- c(1, 2, 3, 4, NA, NA, NA, 8, 9, 10, 11)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5)
    res <- binYonX(X, Y, nBins = 5L,
                   shiftByHalfBinSize = TRUE)
    resM <- imputeLinInterpol(res$y, method = "lin")
    checkEquals(resX[-1], resM[-1])

    ## One missing bin at the end
    Y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, NA, NA)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5,
                              xstart = 1, xend = 11)
    res <- binYonX(X, Y, nBins = 5L,
                   shiftByHalfBinSize = TRUE)
    resM <- imputeLinInterpol(res$y, method = "lin")
    ## checkEquals(resX[-1], resM[-1])
    ## HEEEECK! profBinLin reports 0, which I think is wrong!
    ## -> reported in new_functionality.org

    ## Issue #49
    X <- 1:11
    set.seed(123)
    Y <- sort(rnorm(11, mean = 20, sd = 10))
    ## removing some of the values in the middle
    Y[5:9] <- NA
    nas <- is.na(Y)
    ## Do interpolation with profBinLin:
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5, xstart = min(X),
                              xend = max(X))
    resX
    ## Do the same with binYonX and imputeLinInterpol.
    res <- binYonX(X, Y, nBins = 5L, shiftByHalfBinSize = TRUE)
    resM <- imputeLinInterpol(res$y, method = "lin")
    resM
    ## First bin value is wrong (issue #46)

    ## Visualize
    plot(x = X, y = Y, pch = 16)
    ## Plot the breaks
    abline(v = breaks_on_nBins(min(X), max(X), 5L, TRUE), col = "grey")
    ## Result from profBinLin:
    points(x = res$x, y = resX, col = "blue", type = "b")
    ## Results from imputeLinInterpol
    points(x = res$x, y = resM, col = "green", type = "b",
           pch = 4, lty = 2)



    ## Remove values at the end
    set.seed(123)
    Y <- sort(rnorm(11, mean = 20, sd = 10))
    Y[9:11] <- NA
    nas <- is.na(Y)
    ## Do interpolation with profBinLin:
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5, xstart = min(X),
                              xend = max(X))
    resX
    ## No interpolation performed.
    ## Do the same with binYonX and imputeLinInterpol.
    res <- binYonX(X, Y, nBins = 5L, shiftByHalfBinSize = TRUE)
    resM <- imputeLinInterpol(res$y, method = "lin")
    resM

    ## Results form profBinLin are quite different.

    ## Visualize
    plot(x = X, y = Y, pch = 16, ylim = c(0, max(Y, na.rm = TRUE)))
    ## Plot the breaks
    abline(v = breaks_on_nBins(min(X), max(X), 5L, TRUE), col = "grey")
    ## Result from profBinLin:
    points(x = res$x, y = resX, col = "blue", type = "b")
    ## Results from imputeLinInterpol
    points(x = res$x, y = resM, col = "green", type = "b",
           pch = 4, lty = 2)
    ## END ISSUE


    ## Two missing in the middle.
    Y <- c(1, 2, 3, 4, NA, NA, NA, NA, NA, 10, 11)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5,
                              xstart = 1, xend = 11)
    res <- binYonX(X, Y, nBins = 5L,
                   shiftByHalfBinSize = TRUE)
    resM <- imputeLinInterpol(res$y, method = "lin")
    checkIdentical(resX[-1], resM[-1])

    ## Larger Test:
    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    a <- xcms:::profBinLin(X, X, 1000)
    a2 <- xcms:::profBin(X, X, 1000)
    b <- binYonX(X, X, nBins = 1000, sortedX = TRUE,
                 shiftByHalfBinSize = TRUE)
    bres <- imputeLinInterpol(b$y, method = "lin")
    checkIdentical(a, bres)

    ## Benchmark
    ## myProfBinLin <- function(x, nBins) {
    ##     res <- binYonX(x, nBins = nBins, sortedX = TRUE,
    ##                           shiftByHalfBinSize = TRUE)
    ##     return(imputeLinInterpol(res$y, method = "lin"))
    ## }
    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBinLin(X, X, 1000),
    ##                myProfBinLin(X, 1000))
}

############################################################
## Compare with profBinLinM
## profBinLinM is now deprecated!
dontrun_test_compare_with_profBinLinM <- function() {

    ## A more realistic example.
    set.seed(18011977)
    X <- c(sort(abs(rnorm(3545, mean = 400, sd = 100))),
           sort(abs(rnorm(10000, mean = 500, sd = 200))),
           sort(abs(rnorm(7500, mean = 300, sd = 134))))
    fromX <- c(1, 3546, 13546)
    toX <- c(3545, 13545, 21045)
    nBins <- 1000L
    scIndex <- fromX - 1L
    fX <- min(X)
    tX <- max(X)
    ## Calculate
    resX <- xcms:::profBinLinM(X, X, scIndex, nBins, xstart = fX, xend = tX)
    resM <- binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                    fromIdx = fromX, toIdx = toX,
                    shiftByHalfBinSize = TRUE, sortedX = TRUE)
    ## imputation:
    vals <- lapply(resM, function(z) {
        return(imputeLinInterpol(z$y, method = "lin"))
    })
    M <- do.call(cbind, vals)
    ## Can not compare that sucker!!!
    doplot <- FALSE
    if (doplot) {
        ## Evaluate with a plot.
        resM2 <- binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                         fromIdx = fromX, toIdx = toX,
                         shiftByHalfBinSize = TRUE, sortedX = TRUE)
        ## Get the first non-0 element in the first sub-set
        idx <- which(resM2[[1]]$y > 0)[1]
        ##plot(resM[[1]]$x[1:50], resM[[1]]$y[1:50])
        plot(resM[[1]]$x[1:50], vals[[1]][1:50])
        points(resM2[[1]]$x[1:50], resM2[[1]]$y[1:50], type = "h", col = "red")

        ## Plot the results:
        plot(x = 1:nrow(resX), resX[, 1], cex = 0.3)
        ## The binned data without interpolation:
        points(x = 1:nrow(resX), resM[[1]]$y, col = "red", cex = 0.5, pch = 4)
        ## Interpolated:
        points(x = 1:nrow(resX), M[, 1], col = "green", type = "l")
    }

    ## Benchmark.
    ## myProfBinLinM <- function(x, ...) {
    ##     res <- binYonX(x, ...)
    ##     return(lapply(res, function(z) {
    ##         return(imputeLinInterpol(z$y, method = "lin"))
    ##     }))
    ## }
    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBinLinM(X, X, scIndex, nBins, xstart = fX, xend = tX),
    ##                myProfBinLinM(X, nBins = nBins, binFromX = fX, binToX = tX,
    ##                               fromIdx = fromX, toIdx = toX,
    ##                               shiftByHalfBinSize = TRUE, sortedX = TRUE))
    ## We're faster!
}



## The point is that, for some settings, we get different results in the matchedFilter
## call depending on whether the full profile matrix is generated once, or if it's
## created in chunks (buffers). Actually, this can be nailed down to differences in
## the binning results. To explain these is a little bit more complicated, but is
## most likely due to numeric representation (or rather uncertainty). Thus the bin size
## or the breaks defining the bins can be slightly different for the `profBinM` call.
## This does not seem to be present for the binYonX implementation.
## This has been reported as an issue.
dontrun_tackle_binning_differences <- function() {

    library(xcms)
    library(RUnit)

    xr <- deepCopy(faahko_xr_1)
    mz <- xr@env$mz
    int <- xr@env$intensity
    scantime <- xr@scantime
    scanindex <- xr@scanindex

    ## The code below is extracted from the corresponding functions in xcms.
    step <- 0.2
    ## Create the full buffer in one go:
    mrange <- range(mz)
    mass <- seq(floor(mrange[1]/step)*step, ceiling(mrange[2]/step)*step,
                by = step)
    buf <- xcms:::profBinM(mz, int, scanindex, length(mass),
                           min(mass), max(mass), TRUE)
    ## Just generate a subset.
    buf_subset <- xcms:::profBinM(mz, int, scanindex, 100L, mass[1],
                                  mass[100], TRUE)
    checkEquals(buf[1:100, ], buf_subset)
    checkEquals(buf[1:100, 423], buf_subset[, 423])  ## FAILS!

    ## Using binYonX:
    toIdx <- cumsum(diff(c(scanindex, length(mz))))
    fromIdx <- c(1L, toIdx[-length(toIdx)] + 1L)
    ## Full matrix
    myBuf <- binYonX(mz, int, fromIdx = fromIdx, toIdx = toIdx,
                     binFromX = min(mass), binToX = max(mass),
                     binSize = 0.2, baseValue = 0,
                     shiftByHalfBinSize = TRUE, sortedX = TRUE)
    myBuf <- do.call(cbind, lapply(myBuf, function(z) z$y))
    checkEquals(myBuf, buf)  ## OK
    checkEquals(myBuf[1:100, ], buf_subset)  ## FAILS!
    ## subset
    myBuf_subset <- binYonX(mz, int, fromIdx = fromIdx, toIdx = toIdx,
                            binFromX = min(mass), binToX = mass[100],
                            binSize = 0.2, baseValue = 0,
                            shiftByHalfBinSize = TRUE, sortedX = TRUE)
    myBuf_subset <- do.call(cbind, lapply(myBuf_subset, function(z) z$y))
    checkEquals(myBuf_subset, myBuf[1:100, ])  ## OK
    ## Second subset
    myBuf_subset <- binYonX(mz, int, fromIdx = fromIdx, toIdx = toIdx,
                            binFromX = mass[101], binToX = mass[200],
                            binSize = 0.2, baseValue = 0,
                            shiftByHalfBinSize = TRUE, sortedX = TRUE)
    myBuf_subset <- do.call(cbind, lapply(myBuf_subset, function(z) z$y))
    checkEquals(myBuf_subset, myBuf[101:200, ])  ## OK
    ## Conclusion: binYonX is "stable".
}


############################################################
## Some performance tests...
## sortSplitTest <- function(){
##     out <- double(100)
##     out[sort(unique(idx))] <- unlist(lapply(split(ints, idx), max), use.names = FALSE)
##     return(out)
## }
## sortSplitTest2 <- function(){
##     out <- double(100)
##     uidx <- unique(idx)
##     out[uidx] <- unlist(lapply(split(ints, factor(idx, levels = uidx)), max), use.names = FALSE)
##     return(out)
## }
## First outperforms second!!!

.benchmark_binning <- function() {
    ## Compare new binning functions with profBin.
    library(xcms)
    library(microbenchmark)

    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 300)))
    a <- xcms:::profBin(X, X, 1000)
    b <- binYonX(X, X, nBins = 1000, sortedX = TRUE)

    microbenchmark(xcms:::profBin(X, X, 5),
                   binYonX(X, X, nBins = 1000, sortedX = TRUE))
    ## So, we're 2 times faster.
    microbenchmark(binYonX(X, X, binSize = 1.65, sortedX = TRUE),
                   binYonX(X, X, nBins = 1000, sortedX = TRUE))
    ## binLin:
    microbenchmark(xcms:::profBinLin(X, X, 1000),
                   binYonX(X, X, nBins = 1000, sortedX = TRUE,
                           impute = "lin"))
    ## Also here we're 2 times faster.
}
