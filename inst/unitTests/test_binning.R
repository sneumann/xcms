## Test functions to evaluate binning of values.


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
    brks <- xcms:::breaks_on_nBins(1, 11, 5L, TRUE)

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
    head(xcms:::binYonX(X, X, nBins = 200, shiftByHalfBinSize = TRUE, method = "min")$y)

    ## We've got 1 missing bin.
    X <- c(1, 2, 7, 8, 9, 10)
    xcms:::profBinLin(X, X, 5L)
    xcms:::profBinLinBase(X, X, 5L)

}

############################################################
## Test profBinLinBase
test_profBinLinBase <- function() {
    ## According to the help:
    ## o baselevel: half of minimal signal in the data.
    ## o basespace: 0.075. Linear interpolation is used for data points falling
    ##   within 2*basespace of each other. -> m/z length after which the
    ##   signal drops to baselevel.

    library(xcms)
    library(RUnit)

    X <- 1:16
    Y <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
    nas <- is.na(Y)
    brks <- xcms:::breaks_on_nBins(1, 16, nBins = 16, shiftByHalfBinSize = TRUE)
    resX <- xcms:::profBinLinBase(X[!nas], Y[!nas], 16)
    checkEquals(resX, xcms:::binYonX(X, Y, nBins = 16, baseValue = 0.5,
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

    res <- xcms:::binYonX(X, Y, nBins = 5L)
    checkEquals(res$y, c(2, 3, 6, 8, 10))
    Y[3] <- NA
    res <- xcms:::binYonX(X, Y, nBins = 5L, baseValue = 0)
    checkEquals(res$y, c(2, 0, 6, 8, 10))

    ## Initialize with NA
    res <- xcms:::binYonX(X, Y, nBins = 5L, baseValue = NA)
    checkEquals(res$y, c(2, NA, 6, 8, 10))

    X <- 1:10
    res <- xcms:::binYonX(X, X, binSize = 0.5, baseValue = 0)
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0,
                         8, 0, 9, 10))
    res <- xcms:::binYonX(X, X, binSize = 0.5)
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
    res <- xcms:::binYonX(X, Y, nBins = 5L)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(min(X), max(X), 5L)))
    res <- xcms:::binYonX(X, Y, nBins = 5L, returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- xcms:::binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2.05, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(1, 10, 5L)))
    res <- xcms:::binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10,
                          returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  Define fromIdx, toIdx; binning will only take place in that sub-set.
    res <- xcms:::binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = 10)
    brks <- xcms:::breaks_on_nBins(min(X[1]), max(X[10]), nBins = 5L)
    checkEquals(res$y, c(2.05, 4, 6, 8, 10))
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                          returnIndex = TRUE)
    checkEquals(res$index, c(2, 4, 6, 8, 10))
    ##  If we provide fromX and toX the result will be different, as the breaks
    ##  are calculated differently.
    ##  This calculates 5 bins on range(X)
    res <- xcms:::binYonX(X, X, nBins = 5L, binFromX = min(X), binToX = max(X),
                          fromIdx = 1, toIdx = 10)
    brks <- xcms:::breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkEquals(res$y, c(3, 5, 8, 10, NA)) ## because toIdx = 10 the last is 0
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, nBins = 5L, binFromX = min(X), binToX = max(X),
                          fromIdx = 1, toIdx = 10, returnIndex = TRUE)
    checkEquals(res$index, c(3, 5, 8, 10, NA))
    ##  Check the fromIdx and toIdx again:
    X <- rep(1:10, 3)
    Y <- 1:30
    brks <- xcms:::breaks_on_nBins(1, 10, nBins = 5L)
    ##  In this case we have to set sortedX = TRUE
    res <- xcms:::binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                              fromIdx = 11, toIdx = 20, sortedX = TRUE)
    checkEquals(res$y, c(12, 14, 16, 18, 20))
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                          fromIdx = 11, toIdx = 20, sortedX = TRUE,
                          returnIndex = TRUE)
    checkEquals(res$index, c(12, 14, 16, 18, 20))
    ##  re-order X
    X <- sample(X, size = length(X))
    res <- xcms:::binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(min(X), max(X), 5L)))

    ## Exceptions:
    X <- 1:10
    checkException(xcms:::binYonX(X, X, nBins = -4))
    checkException(xcms:::binYonX(X, X, nBins = 0))
    checkException(xcms:::binYonX(X, X, nBins = 4, fromIdx = 0))
    checkException(xcms:::binYonX(X, X, nBins = 4, toIdx = 0))
    checkException(xcms:::binYonX(X, X, nBins = 4, fromIdx = 7, toIdx = 3))
    checkException(xcms:::binYonX(X, X, nBins = 4, toIdx = 30))

    ## o binSize
    X <- 1:10
    res <- xcms:::binYonX(X, X, binSize = 0.5, baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, binSize = 0.5)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, X, binSize = 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.9, 3.9-5.35, 5.35-6.8, 6.8-8.25, 8.25-9.7, 9.7-10
    checkEquals(res$y, c(2, 3, 5, 6, 8, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, binSize = 1.45)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, X, binSize = 1.45, binFromX = 3, binToX = 9)
    checkEquals(res$y, c(4, 5, 7, 8, 9))
    brks <- xcms:::breaks_on_binSize(3, 9, 1.45)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX(X, X, binSize = 1.45, fromIdx = 3, toIdx = 10,
                          binFromX = min(X), binToX = max(X), baseValue = 0)
    ## bins:
    ## 1-2.45, 2.45-3.90, 3.90-5.35, 5.35-6.80, 6.80-8.25, 8.25-9.70, 9.70-10.00
    checkEquals(res$y, c(0, 3, 5, 6, 8, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, 1.45)
    checkEquals(res$x, breakMidPoint(brks))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks)
    checkEquals(res$y, xcms:::binYonX(X, X, nBins = 4L)$y)
    checkEquals(res$x, breakMidPoint(brks))
    ##  same breaks but sub-setting x
    res <- xcms:::binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                          baseValue = 0)
    checkEquals(res$y, c(0, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- xcms:::binYonX(X, X, breaks = brks)
    checkEquals(res$y, c(4, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))
}

## Test binning using min
test_binning_min <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- xcms:::binYonX(X, nBins = 5L, method = "min")
    checkEquals(res$y, c(1, 3, 5, 7, 9))
    res <- xcms:::binYonX(X, nBins = 5L, method = "min", returnIndex = TRUE)
    checkEquals(res$index, c(1, 3, 5, 7, 9))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- xcms:::binYonX(X, nBins = 5L, binFromX = 1, binToX = 10, method = "min")
    checkEquals(res$y, c(1, 3, 5, 7, 9))
    ##  Define fromIdx, toIdx
    res <- xcms:::binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10, method = "min",
                          binFromX = min(X), binToX = max(X),
                          baseValue = 0)
    brks <- xcms:::breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkEquals(res$y, c(1, 4, 6, 9, 0))
    res <- xcms:::binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                          method = "min", binFromX = min(X), binToX = max(X),
                          returnIndex = TRUE)
    checkEquals(res$index, c(1, 4, 6, 9, NA))

    ## o binSize
    X <- 1:10
    res <- xcms:::binYonX(X, binSize = 0.5, method = "min", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    res <- xcms:::binYonX(X, binSize = 1.45, method = "min")
    xcms:::breaks_on_binSize(min(X), max(X), binSize = 1.45)
    checkEquals(res$y, c(1, 3, 4, 6, 7, 9, 10))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks, method = "min")
    checkEquals(res$y, c(1, 4, 6, 8))
    ##  same breaks but sub-setting x
    res <- xcms:::binYonX(X, breaks = brks, fromIdx = 4, toIdx = 9, method = "min",
                          baseValue = 0)
    checkEquals(res$y, c(0, 4, 6, 8))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks, method = "min")
    checkEquals(res$y, c(3, 5, 6, 8))
}

## Test binning using sum
test_binning_sum <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- xcms:::binYonX(X, nBins = 5L, method = "sum")
    xcms:::breaks_on_nBins(1, 10, nBins = 5)
    checkEquals(res$y, c(3, 7, 11, 15, 19))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- xcms:::binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                          method = "sum")
    checkEquals(res$y, c(3.05, 7, 11, 15, 19))
    res <- xcms:::binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                          method = "sum", returnIndex = TRUE)
    checkEquals(length(res), 2)
    ##  Define fromIdx, toIdx
    res <- xcms:::binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                          method = "sum", binFromX = min(X),
                          binToX = max(X), baseValue = 0)
    xcms:::breaks_on_nBins(min(X), max(X), nBins = 5L)
    checkEquals(res$y, c(6.05, 9, 21, 19, 0))

    ## o binSize
    X <- 1:10
    res <- xcms:::binYonX(X, binSize = 0.5, method = "sum", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    res <- xcms:::binYonX(X, binSize = 1.45, method = "sum")
    xcms:::breaks_on_binSize(min(X), max(X), binSize = 1.45)
    checkEquals(res$y, c(3, 3, 9, 6, 15, 9, 10))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks, method = "sum")
    checkEquals(res$y, c(6, 9, 13, 27))
    ##  same breaks but sub-setting x
    res <- xcms:::binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                          method = "sum", baseValue = 0)
    checkEquals(res$y, c(0, 9, 13, 17))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks, method = "sum")
    checkEquals(res$y, c(7, 5, 13, 17))
}

## Test binning using mean
test_binning_mean <- function() {

    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- xcms:::binYonX(X, nBins = 5L, method = "mean")
    xcms:::breaks_on_nBins(1, 10, nBins = 5)
    checkEquals(res$y, c(1.5, 3.5, 5.5, 7.5, 9.5))
    Y <- X
    Y[3] <- NA
    res <- xcms:::binYonX(X, Y, nBins = 5L, method = "mean")
    checkEquals(res$y, c(1.5, 4, 5.5, 7.5, 9.5))

    ## o binSize
    X <- 1:10
    res <- xcms:::binYonX(X, binSize = 0.5, method = "mean", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- xcms:::binYonX(X, breaks = brks, method = "mean")
    checkEquals(res$y, c(2, 4.5, 6.5, 9))
}

test_breaks <- function() {

    library(xcms)
    library(RUnit)

    ## Test generation of breaks for binning.
    ## o nBins
    res <- xcms:::breaks_on_nBins(1, 10, 4)
    checkIdentical(res, seq(1, 10, length.out = 5))
    res <- xcms:::breaks_on_nBins(2, 8, 20)
    checkEquals(res, seq(2, 8, length.out = 21))

    ## o binSize
    res <- xcms:::breaks_on_binSize(1, 10, 0.13)
    checkEquals(res, c(seq(1, 10, by = 0.13), 10))

}

############################################################
## Test imputation by linear interpolation
test_imputation_lin <- function() {

    library(xcms)
    library(RUnit)

    X <- 1:11
    brks <- xcms:::breaks_on_nBins(1, 11, 5L)

    Y <- c(1, NA, NA, NA, 5, 6, NA, NA, 9, 10, 11)
    ## No imputation
    resVals <- xcms:::binYonX(X, Y, nBins = 5L)$y
    ## Expect NA for bins 2 and 4
    checkEquals(resVals, c(1, NA, 6, NA, 11))
    res <- xcms:::imputeLinInterpol(resVals, method = "lin")
    checkEquals(res, c(1, 3.5, 6, 8.5, 11))

    ## Empty bin at the start.
    Y <- c(NA, 4, 6, 8, 11)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 8, 11))
    Y <- c(NA, 3, 6, 8, 11)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(1.5, 3, 6, 8, 11))

    ## X <- 1:11
    ## Y <- c(NA, NA, 3, 4, 5, 6, 7, 8, 9, 10, 11)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(2, 4, 6, 8, 11))
    ## Y <- c(NA, NA, 3, NA, 5, 6, 7, 8, 9, 10, 11)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(1.5, 3, 6, 8, 11))

    ## Two empty bins at the start.
    Y <- c(NA, NA, 5, 8, 11)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(5/3, 2*5/3, 5, 8, 11))
    ## X <- 1:11
    ## Y <- c(NA, NA, NA, NA, 5, NA, 7, 8, 9, 10, 11)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "no")
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(5/3, 2*5/3, 5, 8, 11))
    ## Y <- c(NA, NA, NA, NA, 6, NA, 7, 8, 9, 10, 11)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(2, 4, 6, 8, 11))

    ## Empty bin at the end.
    Y <- c(2, 4, 6, 9, NA)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 9, 4.5))
    ## Y <- c(1, 2, 3, 4, 5, 6, 9, 7.2, NA, NA, NA)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "no")
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(2, 4, 6, 9, 4.5))

    ## Two empty bins at the end.
    Y <- c(2, 4, 6, NA, NA)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 4, 2))
    ## Y <- c(1, 2, 3, 4, 5, 6, NA, NA, NA, NA, NA)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "no")
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(2, 4, 6, 4, 2))

    ## 2 consecutive empty bins.
    Y <- c(2, NA, NA, 8, 10)
    res <- xcms:::imputeLinInterpol(Y, method = "lin")
    checkEquals(res, c(2, 4, 6, 8, 10))
    ## Y <- c(2, 1, NA, NA, NA, NA, 7, 8, 9, 10, NA)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(2, 4, 6, 8, 10))
    ## Y <- c(NA, 1, NA, NA, NA, NA, 9, 8, 9, 10, NA)
    ## res <- xcms:::binYonX(X, Y, nBins = 5L, impute = "lin")
    ## checkEquals(res$y, c(1, 1 + 8/3, 1 + 2*8/3, 9, 10))
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase",
                                    baseValue = 0.5, distance = 0)
    dummy <- Y
    dummy[is.na(dummy)] <- 0.5
    checkEquals(res, dummy)

    ## Modify basespace such that 1 neighboring bin is considered.
    res <- xcms:::imputeLinInterpol(Y, method = "linbase")
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5 + 1.5/2, 2)
    checkEquals(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## Modify basespace such that 2 neighboring bins are considered.
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 2L)
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 1)
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 1)
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 1)
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 1,
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
    res <- xcms:::imputeLinInterpol(Y, method = "linbase", distance = 2,
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
    resM <- xcms:::binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                           toIdx = toIdx)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], xcms:::binYonX(X, X, nBins = nBins,
                                              fromIdx = fromIdx[i],
                                              toIdx = toIdx[i]))
    }

    ## GC-torture; check if objects in C are save.
    gctorture(on = TRUE)
    resM <- xcms:::binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                           toIdx = toIdx)
    gctorture(on = FALSE)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], xcms:::binYonX(X, X, nBins = nBins,
                                              fromIdx = fromIdx[i],
                                              toIdx = toIdx[i]))
    }

    X <- 1:30
    fromIdx <- c(1, 3, 14)
    toIdx <- c(8, 11, 28)
    resM <- xcms:::binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                                 toIdx = toIdx)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], xcms:::binYonX(X, X, nBins = nBins,
                                              fromIdx = fromIdx[i],
                                              toIdx = toIdx[i]))
    }
    ## Returning also the index
    resM <- xcms:::binYonX(X, nBins = nBins, fromIdx = fromIdx, toIdx = toIdx,
                           returnIndex = TRUE)
    checkTrue(all(unlist(lapply(resM, function(z) {
        return(any(names(z) == "index"))
    })) == TRUE))

    ## With defining fromX and toX before. The breaks will thus
    ## be defined on these values.
    fromX <- 1
    toX <- 30
    resM <- xcms:::binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                           toIdx = toIdx, binFromX = fromX, binToX = toX,
                           baseValue = 17)
    checkEquals(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        checkEquals(resM[[i]], xcms:::binYonX(X, X, nBins = nBins,
                                              fromIdx = fromIdx[i],
                                              toIdx = toIdx[i],
                                              binFromX = fromX,
                                              binToX = toX,
                                              baseValue = 17))
    }

    ## Error checks.
    checkException(xcms:::binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = c(2, 3)))
}


############################################################
##
##        COMPARE WITH profBin METHODS
##
############################################################

############################################################
## Compare binYonX, max, with plain profBin
test_compare_with_profBin <- function() {

    library(xcms)
    library(RUnit)

    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    a <- xcms:::profBin(X, X, 1000)
    b <- xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE,
                        shiftByHalfBinSize = TRUE,
                        baseValue = 0)
    checkEquals(a, b$y)

    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBin(X, X, 1000),
    ##                xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE,
    ##                               shiftByHalfBinSize = TRUE,
    ##                               baseValue = 0))

    ## To be on the save side, compare also with binSize:
    binS <- unique(diff(b$x))[1]
    c <- xcms:::binYonX(X, binSize = binS, sortedX = TRUE,
                        shiftByHalfBinSize = TRUE,
                        baseValue = 0)
    checkIdentical(length(b$x), length(c$x))
    ## The bins will be the same except for the last one.
    checkEquals(b$x[-length(b$x)], c$x[-length(c$x)])
    checkEquals(b$y, c$y)

    ## NA handling
    X[20:30] <- NA
    checkException(a <- xcms:::profBin(X, X, 1000))
    checkException(b <- xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE,
                                       shiftByHalfBinSize = TRUE))
    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    Y <- X
    Y[20:30] <- NA
    ## Doesn't work! a <- xcms:::profBin(X, Y, 1000)
}


############################################################
## Compare with profBinM
test_compare_with_profBinM <- function() {

    ## profBinM: takes vectors of data and a second one specifying
    ## non-overlapping sub-sets.
    X <- 1:30
    fromX <- c(1, 11, 19)
    toX <- c(10, 18, 30)
    scIndex <- fromX - 1L
    nBins <- 5L

    resX <- xcms:::profBinM(X, X, scIndex, nBins)
    resM <- xcms:::binYonX(X, nBins = nBins, binFromX = min(X), binToX = max(X),
                           fromIdx = fromX, toIdx = toX,
                           shiftByHalfBinSize = TRUE)
    M <- do.call(cbind, lapply(resM, function(x) {x$y}))
    checkEquals(resX, M)

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
    resM <- xcms:::binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                           fromIdx = fromX, toIdx = toX,
                           shiftByHalfBinSize = TRUE, sortedX = TRUE)
    M <- do.call(cbind, lapply(resM, function(x) {x$y}))
    checkEquals(resX, M)

    ## Benchmark.
    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBinM(X, X, scIndex, nBins, xstart = fX, xend = tX),
    ##                xcms:::binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
    ##                               fromIdx = fromX, toIdx = toX,
    ##                               shiftByHalfBinSize = TRUE, sortedX = TRUE))
    ## we're slower speed.
}



############################################################
## Compare with profBinLin
test_compare_with_profBinLin <- function() {
    ## binLin does linear interpolation of values, in which nothing
    ## was binned... let's see.
    ## binLin does the imputation on the already binned values, this was also
    ## implemented in the impute = "lin" method.

    library(xcms)
    library(RUnit)

    X <- 1:11
    brks <- xcms:::breaks_on_nBins(1, 11, 5L, TRUE)

    ## One missing bin in the middle
    Y <- c(1, 2, 3, 4, NA, NA, NA, 8, 9, 10, 11)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5)
    res <- xcms:::binYonX(X, Y, nBins = 5L,
                          shiftByHalfBinSize = TRUE)
    resM <- xcms:::imputeEmptyBins(res$y, method = "lin")
    checkEquals(resX[-1], resM[-1])

    ## One missing bin at the end
    Y <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, NA, NA)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5,
                              xstart = 1, xend = 11)
    res <- xcms:::binYonX(X, Y, nBins = 5L,
                          shiftByHalfBinSize = TRUE)
    resM <- xcms:::imputeEmptyBins(res$y, method = "lin")
    ## checkEquals(resX[-1], res$y[-1])
    ## HEEEECK! profBinLin reports 0, which I think is wrong!

    ## Two missing in the middle.
    Y <- c(1, 2, 3, 4, NA, NA, NA, NA, NA, 10, 11)
    nas <- is.na(Y)
    resX <- xcms:::profBinLin(X[!nas], Y[!nas], 5,
                              xstart = 1, xend = 11)
    res <- xcms:::binYonX(X, Y, nBins = 5L,
                          shiftByHalfBinSize = TRUE)
    resM <- xcms:::imputeEmptyBins(res$y, method = "lin")
    checkEquals(resX[-1], resM[-1])

    ## Larger Test:
    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    a <- xcms:::profBinLin(X, X, 1000)
    a2 <- xcms:::profBin(X, X, 1000)
    b <- xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE,
                        shiftByHalfBinSize = TRUE)
    bres <- xcms:::imputeEmptyBins(b$y, method = "lin")
    checkEquals(a, bres)
    ## AMAZING!!!

    ## Benchmark
    ## myProfBinLin <- function(x, nBins) {
    ##     res <- xcms:::binYonX(x, nBins = nBins, sortedX = TRUE,
    ##                           shiftByHalfBinSize = TRUE)
    ##     return(xcms:::imputeEmptyBins(res$y, method = "lin"))
    ## }
    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBinLin(X, X, 1000),
    ##                myProfBinLin(X, 1000))
}

############################################################
## Compare with profBinLinM
test_compare_with_profBinLinM <- function() {

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
    resX <- xcms:::profBinLinM(X, X, scIndex, nBins, xstart = fX, xend = tX)
    resM <- xcms:::binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                           fromIdx = fromX, toIdx = toX,
                           shiftByHalfBinSize = TRUE, sortedX = TRUE)
    ## imputation:
    vals <- lapply(resM, function(z) {
        return(xcms:::imputeEmptyBins(z$y, method = "lin"))
    })
    M <- do.call(cbind, vals)
    ## Can not compare that sucker!!!
    doplot <- FALSE
    if (doplot) {
        ## Evaluate with a plot.
        resM2 <- xcms:::binYonX(X, nBins = nBins, binFromX = fX, binToX = tX,
                                fromIdx = fromX, toIdx = toX,
                                shiftByHalfBinSize = TRUE, sortedX = TRUE)
        ## Get the first non-0 element in the first sub-set
        idx <- which(resM2[[1]]$y > 0)[1]
        ##plot(resM[[1]]$x[1:50], resM[[1]]$y[1:50])
        plot(resM[[1]]$x[1:50], vals[[1]][1:50])
        points(resM2[[1]]$x[1:50], resM2[[1]]$y[1:50], type = "h", col = "red")
    }
    ## checkEquals(resX[, 2:3], M[, 2:3])

    ## Benchmark.
    ## myProfBinLinM <- function(x, ...) {
    ##     res <- xcms:::binYonX(x, ...)
    ##     return(lapply(res, function(z) {
    ##         return(xcms:::imputeEmptyBins(z$y, method = "lin"))
    ##     }))
    ## }
    ## library(microbenchmark)
    ## microbenchmark(xcms:::profBinLinM(X, X, scIndex, nBins, xstart = fX, xend = tX),
    ##                myProfBinLinM(X, nBins = nBins, binFromX = fX, binToX = tX,
    ##                               fromIdx = fromX, toIdx = toX,
    ##                               shiftByHalfBinSize = TRUE, sortedX = TRUE))
    ## We're faster!
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
    b <- xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE)
    microbenchmark(xcms:::profBin(X, X, 5),
                   xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE))
    ## So, we're 2 times faster.
    microbenchmark(xcms:::binYonX(X, X, binSize = 1.65, sortedX = TRUE),
                   xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE))
    ## binLin:
    microbenchmark(xcms:::profBinLin(X, X, 1000),
                   xcms:::binYonX(X, X, nBins = 1000, sortedX = TRUE,
                                  impute = "lin"))
    ## Also here we're 2 times faster.
}
