test_that("binYonX NA handling works", {
    ## NA values in y should be ignored internally.
    X <- 1:10
    Y <- X
    Y[4] <- NA
    res <- binYonX(X, Y, nBins = 5L)
    expect_equal(res$y, c(2, 3, 6, 8, 10))
    Y[3] <- NA
    res <- binYonX(X, Y, nBins = 5L, baseValue = 0)
    expect_equal(res$y, c(2, 0, 6, 8, 10))
    ## Initialize with NA
    res <- binYonX(X, Y, nBins = 5L, baseValue = NA)
    expect_equal(res$y, c(2, NA, 6, 8, 10))
    X <- 1:10
    res <- binYonX(X, X, binSize = 0.5, baseValue = 0)
    expect_equal(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0,
                          8, 0, 9, 10))
    res <- binYonX(X, X, binSize = 0.5)
    expect_equal(res$y, c(1, NA, 2, NA, 3, NA, 4, NA, 5, NA, 6, NA, 7,
                          NA, 8, NA, 9, 10))
})

test_that("binYonX max works", {
    X <- 1:10
    Y <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }
    
    ## o nBins
    res <- binYonX(X, Y, nBins = 5L)
    expect_equal(res$y, c(2, 4, 6, 8, 10))
    expect_equal(res$x,
                 breakMidPoint(breaks_on_nBins(min(X), max(X), 5L)))
    res <- binYonX(X, Y, nBins = 5L, returnIndex = TRUE)
    expect_equal(res$index, c(2, 4, 6, 8, 10))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    expect_equal(res$y, c(2.05, 4, 6, 8, 10))
    expect_equal(res$x,
                 breakMidPoint(breaks_on_nBins(1, 10, 5L)))
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10,
                   returnIndex = TRUE)
    expect_equal(res$index, c(2, 4, 6, 8, 10))
    ##  Define fromIdx, toIdx; binning will only take place in that sub-set.
    res <- binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = 10)
    brks <- breaks_on_nBins(min(X[1]), max(X[10]), nBins = 5L)
    expect_equal(res$y, c(2.05, 4, 6, 8, 10))
    expect_equal(res$x, breakMidPoint(brks))
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   returnIndex = TRUE)
    expect_equal(res$index, c(2, 4, 6, 8, 10))
    ##  If we provide fromX and toX the result will be different, as the breaks
    ##  are calculated differently.
    ##  This calculates 5 bins on range(X)
    res <- binYonX(X, X, nBins = 5L, binFromX = min(X), binToX = max(X),
                   fromIdx = 1, toIdx = 10)
    brks <- breaks_on_nBins(min(X), max(X), nBins = 5L)
    expect_equal(res$y, c(3, 5, 8, 10, NA)) ## because toIdx = 10 the last is 0
    expect_equal(res$x, breakMidPoint(brks))
    res <- binYonX(X, nBins = 5L, binFromX = min(X), binToX = max(X),
                   fromIdx = 1, toIdx = 10, returnIndex = TRUE)
    expect_equal(res$index, c(3, 5, 8, 10, NA))
    ##  Check the fromIdx and toIdx again:
    X <- rep(1:10, 3)
    Y <- 1:30
    brks <- breaks_on_nBins(1, 10, nBins = 5L)
    ##  In this case we have to set sortedX = TRUE
    res <- binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                   fromIdx = 11, toIdx = 20, sortedX = TRUE)
    expect_equal(res$y, c(12, 14, 16, 18, 20))
    expect_equal(res$x, breakMidPoint(brks))
    res <- binYonX(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                   fromIdx = 11, toIdx = 20, sortedX = TRUE,
                   returnIndex = TRUE)
    expect_equal(res$index, c(12, 14, 16, 18, 20))
    ##  re-order X
    X <- sample(X, size = length(X))
    res <- binYonX(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    expect_equal(res$y, c(2, 4, 6, 8, 10))
    expect_equal(res$x,
                 breakMidPoint(breaks_on_nBins(min(X), max(X), 5L)))

    ## Exceptions:
    X <- 1:10
    expect_error(binYonX(X, X, nBins = -4))
    expect_error(binYonX(X, X, nBins = 0))
    expect_error(binYonX(X, X, nBins = 4, fromIdx = 0))
    expect_error(binYonX(X, X, nBins = 4, toIdx = 0))
    expect_error(binYonX(X, X, nBins = 4, fromIdx = 7, toIdx = 3))
    expect_error(binYonX(X, X, nBins = 4, toIdx = 30))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, X, binSize = 0.5, baseValue = 0)
    brks <- breaks_on_binSize(1, 10, binSize = 0.5)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    expect_identical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8,
                              0, 9, 10))
    expect_identical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45)
    brks <- breaks_on_binSize(1, 10, binSize = 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.9, 3.9-5.35, 5.35-6.8, 6.8-8.25, 8.25-10
    expect_identical(res$y, c(2, 3, 5, 6, 8, 10))
    expect_identical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45, binFromX = 3, binToX = 9)
    brks <- breaks_on_binSize(3, 9, 1.45)
    expect_identical(res$y, c(4, 5, 7, 9))
    expect_identical(res$x, breakMidPoint(brks))
    ##
    res <- binYonX(X, X, binSize = 1.45, fromIdx = 3, toIdx = 10,
                   binFromX = min(X), binToX = max(X), baseValue = 0)
    brks <- breaks_on_binSize(1, 10, 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.90, 3.90-5.35, 5.35-6.80, 6.80-8.25, 8.25-9.70, 9.70-10.00
    expect_identical(res$y, c(0, 3, 5, 6, 8, 10))
    expect_identical(res$x, breakMidPoint(brks))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks)
    expect_equal(res$y, binYonX(X, X, nBins = 4L)$y)
    expect_equal(res$x, breakMidPoint(brks))
    ##  same breaks but sub-setting x
    res <- binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                   baseValue = 0)
    expect_equal(res$y, c(0, 5, 7, 9))
    expect_equal(res$x, breakMidPoint(brks))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, X, breaks = brks)
    expect_equal(res$y, c(4, 5, 7, 9))
    expect_equal(res$x, breakMidPoint(brks))

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
    expect_equal(res$y, resR)

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
    expect_equal(res$y, resR)
    shift <- FALSE
    mass <- seq(floor(min(xRangeFull)/step)*step,
                ceiling(max(xRangeFull)/step)*step, by = step)
    nBins <- length(mass)
    resR <- xcms:::profBinR(X1, Y1, nBins = nBins, fromX = min(xRangeFull),
                            toX = max(xRangeFull), shiftByHalfBinSize = shift)
    res <- binYonX(X1, Y1, nBins = nBins, binFromX = min(xRangeFull),
                   binToX = max(xRangeFull), shiftByHalfBinSize = shift,
                   baseValue = 0)
    expect_equal(res$y, resR)

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
    expect_equal(res$y, resR)
})

## Test binning using min
test_that("binYonX min works", {
    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "min")
    expect_equal(res$y, c(1, 3, 5, 7, 9))
    res <- binYonX(X, nBins = 5L, method = "min", returnIndex = TRUE)
    expect_equal(res$index, c(1, 3, 5, 7, 9))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10, method = "min")
    expect_equal(res$y, c(1, 3, 5, 7, 9))
    ##  Define fromIdx, toIdx
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10, method = "min",
                   binFromX = min(X), binToX = max(X),
                   baseValue = 0)
    brks <- breaks_on_nBins(min(X), max(X), nBins = 5L)
    expect_equal(res$y, c(1, 4, 6, 9, 0))
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   method = "min", binFromX = min(X), binToX = max(X),
                   returnIndex = TRUE)
    expect_equal(res$index, c(1, 4, 6, 9, NA))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "min", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    breaks_on_binSize(min(X), max(X), binSize = 0.5)
    expect_identical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    ##
    res <- binYonX(X, binSize = 1.45, method = "min")
    breaks_on_binSize(min(X), max(X), binSize = 1.45)
    expect_identical(res$y, c(1, 3, 4, 6, 7, 9))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "min")
    expect_equal(res$y, c(1, 4, 6, 8))
    ##  same breaks but sub-setting x
    res <- binYonX(X, breaks = brks, fromIdx = 4, toIdx = 9, method = "min",
                   baseValue = 0)
    expect_identical(res$y, c(0, 4, 6, 8))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "min")
    expect_identical(res$y, c(3, 5, 6, 8))
})

## Test binning using sum
test_that("binYonX sum works", {
    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "sum")
    breaks_on_nBins(1, 10, nBins = 5)
    expect_equal(res$y, c(3, 7, 11, 15, 19))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                   method = "sum")
    expect_equal(res$y, c(3.05, 7, 11, 15, 19))
    res <- binYonX(X, nBins = 5L, binFromX = 1, binToX = 10,
                   method = "sum", returnIndex = TRUE)
    expect_equal(length(res), 2)
    ##  Define fromIdx, toIdx
    res <- binYonX(X, nBins = 5L, fromIdx = 1, toIdx = 10,
                   method = "sum", binFromX = min(X),
                   binToX = max(X), baseValue = 0)
    breaks_on_nBins(min(X), max(X), nBins = 5L)
    expect_identical(res$y, c(6.05, 9, 21, 19, 0))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "sum", baseValue = 0)
    breaks_on_binSize(min(X), max(X), binSize = 0.5)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    expect_identical(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    res <- binYonX(X, binSize = 1.45, method = "sum")
    breaks_on_binSize(min(X), max(X), binSize = 1.45)
    expect_identical(res$y, c(3, 3, 9, 6, 15, 19))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "sum")
    expect_identical(res$y, c(6, 9, 13, 27))
    ##  same breaks but sub-setting x
    res <- binYonX(X, X, breaks = brks, fromIdx = 4, toIdx = 9,
                   method = "sum", baseValue = 0)
    expect_identical(res$y, c(0, 9, 13, 17))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "sum")
    expect_identical(res$y, c(7, 5, 13, 17))
})


## Test binning using mean
test_that("binYonX mean works", {
    X <- 1:10
    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- binYonX(X, nBins = 5L, method = "mean")
    breaks_on_nBins(1, 10, nBins = 5)
    expect_equal(res$y, c(1.5, 3.5, 5.5, 7.5, 9.5))
    Y <- X
    Y[3] <- NA
    res <- binYonX(X, Y, nBins = 5L, method = "mean")
    expect_equal(res$y, c(1.5, 4, 5.5, 7.5, 9.5))

    ## o binSize
    X <- 1:10
    res <- binYonX(X, binSize = 0.5, method = "mean", baseValue = 0)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    expect_equal(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- binYonX(X, breaks = brks, method = "mean")
    expect_equal(res$y, c(2, 4.5, 6.5, 9))
})

test_that("breaks defining functions work", {
    ## Test generation of breaks for binning.
    ## o nBins
    res <- breaks_on_nBins(1, 10, 4)
    expect_equal(res, seq(1, 10, length.out = 5))
    res <- breaks_on_nBins(2, 8, 20)
    expect_equal(res, seq(2, 8, length.out = 21))
    ##
    brksR <- seq(200, 600, length.out = 2002)
    brks <- breaks_on_nBins(200, 600, nBins = 2001)
    expect_equal(brks, brksR)
    ## Simulate shift by half bin size
    brksR <- seq((200 - 0.1), (600 + 0.1), length.out = 2002)
    brks <- breaks_on_nBins(200, 600, nBins = 2001,
                            shiftByHalfBinSize = TRUE)
    expect_equal(brks, brksR)

    ## o binSize
    res <- breaks_on_binSize(1, 10, 0.13)
    resR <- seq(1, 10, by = 0.13)
    expect_equal(res[-length(res)], resR[-length(resR)])
    expect_equal(res[length(res)], 10)
    ##   Will create one bin more.
    res <- breaks_on_binSize(1, 10, 0.51)
    resR <- seq(1, 10, by = 0.51)
    expect_equal(res[-length(res)], resR)
    expect_equal(res[length(res)], 10)
    ##
    brksR <- seq(200, 600, by = 0.2)
    brks <- breaks_on_binSize(200, 600, binSize = 0.2)
    expect_equal(brks, brksR)
    ## Simulate shift by half bin size
    brksR <- seq((200 - 0.1), (600 + 0.1), by = 0.2)
    brks <- breaks_on_binSize((200 - 0.1), (600 + 0.1),
                              binSize = 0.2)
    expect_equal(brks, brksR)
    ##
    ## Ultimate fix for issue #118 
    brksR <- seq((200 - 0.1), (600), by = 0.2)
    brks <- breaks_on_binSize((200 - 0.1), (600), binSize = 0.2)
    ## Compare them up to the last value, since in R that will be 600-01, while
    ## breaks_on_binSize will ensure that the upper limit (600) is still within
    ## the breaks
    cmn <- 1:(length(brksR) - 1)
    expect_equal(brks[cmn], brksR[cmn])
    ## Now, that below breaks on a windows build machine (issue #127)
    ## expect_true(length(brks) > length(brksR))
    ## expect_equal(brks[-length(brks)], brksR)
    expect_equal(brks[length(brks)], 600)
})

test_that("binYonX with imputation_lin works",  {
    X <- 1:11
    brks <- breaks_on_nBins(1, 11, 5L)
    Y <- c(1, NA, NA, NA, 5, 6, NA, NA, 9, 10, 11)
    ## No imputation
    resVals <- binYonX(X, Y, nBins = 5L)$y
    ## Expect NA for bins 2 and 4
    expect_identical(resVals, c(1, NA, 6, NA, 11))
    res <- imputeLinInterpol(resVals, method = "lin")
    expect_identical(res, c(1, 3.5, 6, 8.5, 11))

    ## Empty bin at the start.
    Y <- c(NA, 4, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_identical(res, c(2, 4, 6, 8, 11))
    Y <- c(NA, 3, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_identical(res, c(1.5, 3, 6, 8, 11))

    ## Two empty bins at the start.
    Y <- c(NA, NA, 5, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_identical(res, c(5/3, 2*5/3, 5, 8, 11))

    ## Empty bin at the end.
    Y <- c(2, 4, 6, 9, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_equal(res, c(2, 4, 6, 9, 4.5))

    ## Two empty bins at the end.
    Y <- c(2, 4, 6, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_equal(res, c(2, 4, 6, 4, 2))

    ## 3 empty bins at the end.
    Y <- c(2, 4, NA, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_equal(res, c(2, 4, 3, 2, 1))

    ## 2 consecutive empty bins.
    Y <- c(2, NA, NA, 8, 10)
    res <- imputeLinInterpol(Y, method = "lin")
    expect_equal(res, c(2, 4, 6, 8, 10))

    ## Simulate the profBinLin: no interpolation at both ends:
    Y <- c(2, NA, NA, 8, 10)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    expect_equal(res, c(2, 4, 6, 8, 10))

    ## One NA at start:
    Y <- c(NA, 4, 6, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    expect_identical(res, c(0, 4, 6, 8, 11))

    ## Two NAs at the start.
    Y <- c(NA, NA, 5, 8, 11)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    expect_identical(res, c(0, 0, 5, 8, 11))

    ## NA at the end.
    Y <- c(2, 4, 6, 9, NA)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    expect_equal(res, c(2, 4, 6, 9, 0))

    ## Two NAs at the end.
    Y <- c(2, 4, 6, NA, NA)
    res <- imputeLinInterpol(Y, method = "lin", noInterpolAtEnds = TRUE)
    expect_equal(res, c(2, 4, 6, 0, 0))
})

test_that("binYonX with imputation_linbin works", {
    doPlot <- FALSE
    ## Construct example:
    ## We're using the same test than we did for profBinLinBase.
    Y <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
    nas <- is.na(Y)
    res <- imputeLinInterpol(Y, method = "linbase",
                             baseValue = 0.5, distance = 0)
    dummy <- Y
    dummy[is.na(dummy)] <- 0.5
    expect_equal(res, dummy)

    ## Modify basespace such that 1 neighboring bin is considered.
    res <- imputeLinInterpol(Y, method = "linbase")
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 3.5/2,
                0.5, 0.5 + 2.5/2, 3, 3 - 2.5/2, 0.5, 0.5, 0.5 + 1.5/2, 2)
    expect_equal(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }

    ## Modify basespace such that 2 neighboring bins are considered.
    res <- imputeLinInterpol(Y, method = "linbase", distance = 2L)
    expect <- c(3, 2, 1, 2, 2+2/3, 2+2*2/3, 4, 4 - 1/4, 4 - 2/4,
                4 - 3/4, 3, 3 - 1/5, 3 - 2/5, 3 - 3/5, 3 - 4/5, 2)
    expect_equal(res, expect)
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
    expect_equal(res, expect)
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
    expect_equal(res, expect)
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
    expect_equal(res, expect)
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
    expect_equal(res, expect)
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
    expect_equal(res, expect)
    if (doPlot) {
        plot(x = 1:length(Y), Y, ylim = c(0, max(Y, na.rm = TRUE)), pch = 16)
        points(x = 1:length(Y), y = expect, col = "blue", type = "b")
        points(x = 1:length(Y), y = res, col = "red", type = "l", lty = 2)
    }
})

## test_that("testIntegerInput and testDoubleInput works", {

##     xcms:::testIntegerInput(4)
##     xcms:::testIntegerInput(c(4, 5))

##     xcms:::testRealInput(4)
##     xcms:::testRealInput(c(4.2231444, 5))

## })

test_that("binYonX on subsets works", {
    ## Simple test without actually needing subsets.
    X <- 1:11
    Y <- 1:11
    nBins = 5L

    X <- 1:30
    fromIdx <- c(1, 3, 14)
    toIdx <- c(8, 11, 28)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    expect_equal(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        expect_equal(resM[[i]], binYonX(X, X, nBins = nBins,
                                        fromIdx = fromIdx[i],
                                        toIdx = toIdx[i]))
    }

    ## GC-torture; check if objects in C are save.
    gctorture(on = TRUE)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    gctorture(on = FALSE)
    for (i in 1:length(fromIdx)) {
        expect_equal(resM[[i]], binYonX(X, X, nBins = nBins,
                                        fromIdx = fromIdx[i],
                                        toIdx = toIdx[i]))
    }

    X <- 1:30
    fromIdx <- c(1, 3, 14)
    toIdx <- c(8, 11, 28)
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx)
    expect_equal(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        expect_equal(resM[[i]], binYonX(X, X, nBins = nBins,
                                        fromIdx = fromIdx[i],
                                        toIdx = toIdx[i]))
    }
    ## Returning also the index
    resM <- binYonX(X, nBins = nBins, fromIdx = fromIdx, toIdx = toIdx,
                    returnIndex = TRUE)
    expect_true(all(unlist(lapply(resM, function(z) {
        return(any(names(z) == "index"))
    })) == TRUE))

    ## With defining fromX and toX before. The breaks will thus
    ## be defined on these values.
    fromX <- 1
    toX <- 30
    resM <- binYonX(X, X, nBins = nBins, fromIdx = fromIdx,
                    toIdx = toIdx, binFromX = fromX, binToX = toX,
                    baseValue = 17)
    expect_equal(length(resM), 3)
    for (i in 1:length(fromIdx)) {
        expect_equal(resM[[i]], binYonX(X, X, nBins = nBins,
                                        fromIdx = fromIdx[i],
                                        toIdx = toIdx[i],
                                        binFromX = fromX,
                                        binToX = toX,
                                        baseValue = 17))
    }
    ## Error checks.
    expect_error(binYonX(X, X, nBins = 5L, fromIdx = 1, toIdx = c(2, 3)))
})
