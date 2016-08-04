## Test functions to evaluate binning of values.

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

test_compare_with_profBin <- function() {

    library(xcms)
    library(RUnit)

    set.seed(18011977)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
    a <- xcms:::profBin(X, X, 1000)
    b <- xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE,
                            shiftByHalfBinSize = TRUE)
    checkEquals(a, b$y)

    ## To be on the save side, compare also with binSize:
    binS <- unique(diff(b$x))[1]
    c <- xcms:::binYonX_max(X, binSize = binS, sortedX = TRUE,
                            shiftByHalfBinSize = TRUE)
    checkIdentical(length(b$x), length(c$x))
    ## The bins will be the same except for the last one.
    checkEquals(b$x[-length(b$x)], c$x[-length(c$x)])
    checkEquals(b$y, c$y)
}

test_binning <- function() {

    library(xcms)
    library(RUnit)

    X <- 1:10
    Y <- 1:10

    breakMidPoint <- function(x) {
        return((x[-1L] + x[-length(x)])/2)
    }

    ## o nBins
    res <- xcms:::binYonX_max(X, Y, nBins = 5L)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(min(X), max(X), 5L)))
    ##  Defining fromX, toX
    X <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5,
           11, 11.124, 11.125, 11.126, 12, 13)
    res <- xcms:::binYonX_max(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2.05, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(1, 10, 5L)))
    ##  Define fromIdx, toIdx
    res <- xcms:::binYonX_max(X, X, nBins = 5L, fromIdx = 1, toIdx = 10)
    brks <- xcms:::breaks_on_nBins(min(X), max(X), nBins = 5L)
    ##  This calculates 5 bins on range(X)
    checkEquals(res$y, c(3, 5, 8, 10, 0)) ## because toIdx = 10 the last is 0
    checkEquals(res$x, breakMidPoint(brks))
    ##  Check the fromIdx and toIdx again:
    X <- rep(1:10, 3)
    Y <- 1:30
    brks <- xcms:::breaks_on_nBins(1, 10, nBins = 5L)
    ##  In this case we have to set sortedX = TRUE
    res <- xcms:::binYonX_max(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
                              fromIdx = 11, toIdx = 20, sortedX = TRUE)
    checkEquals(res$y, c(12, 14, 16, 18, 20))
    checkEquals(res$x, breakMidPoint(brks))
    ##  re-order X
    X <- sample(X, size = length(X))
    res <- xcms:::binYonX_max(X, X, nBins = 5L, binFromX = 1, binToX = 10)
    checkEquals(res$y, c(2, 4, 6, 8, 10))
    checkEquals(res$x,
                breakMidPoint(xcms:::breaks_on_nBins(min(X), max(X), 5L)))

    ## Exceptions:
    X <- 1:10
    checkException(xcms:::binYonX_max(X, X, nBins = -4))
    checkException(xcms:::binYonX_max(X, X, nBins = 0))
    checkException(xcms:::binYonX_max(X, X, nBins = 4, fromIdx = 0))
    checkException(xcms:::binYonX_max(X, X, nBins = 4, toIdx = 0))
    checkException(xcms:::binYonX_max(X, X, nBins = 4, fromIdx = 7, toIdx = 3))
    checkException(xcms:::binYonX_max(X, X, nBins = 4, toIdx = 30))

    ## o binSize
    X <- 1:10
    res <- xcms:::binYonX_max(X, X, binSize = 0.5)
    ##  We expect: every other bin is 0. The last value has to be 10, the first 1.
    ## 1-1.5, 1.5-2, 2-2.5, 2.5-3, 3-3.5, 3.5-4, 4-4.5, 4.5-5, 5-5.5, 5.5-6,
    ## 6-6.5, 6.5-7, 7-7.5, 7.5-8, 8-8.5, 8.5-9, 9-9.5. 9.5-10.
    checkEquals(res$y, c(1, 0, 2, 0, 3, 0, 4, 0, 5, 0, 6, 0, 7, 0, 8, 0, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, binSize = 0.5)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX_max(X, X, binSize = 1.45)
    ## bins:
    ## 1-2.45, 2.45-3.9, 3.9-5.35, 5.35-6.8, 6.8-8.25, 8.25-9.7, 9.7-10
    checkEquals(res$y, c(2, 3, 5, 6, 8, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, binSize = 1.45)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX_max(X, X, binSize = 1.45, binFromX = 3, binToX = 9)
    checkEquals(res$y, c(4, 5, 7, 8, 9))
    brks <- xcms:::breaks_on_binSize(3, 9, 1.45)
    checkEquals(res$x, breakMidPoint(brks))
    res <- xcms:::binYonX_max(X, X, binSize = 1.45, fromIdx = 3, toIdx = 10)
    ## bins:
    ## 1-2.45, 2.45-3.90, 3.90-5.35, 5.35-6.80, 6.80-8.25, 8.25-9.70, 9.70-10.00
    checkEquals(res$y, c(0, 3, 5, 6, 8, 9, 10))
    brks <- xcms:::breaks_on_binSize(1, 10, 1.45)
    checkEquals(res$x, breakMidPoint(brks))

    ## o breaks
    X <- 1:10
    brks <- seq(1, 10, length.out = 5)
    res <- xcms:::binYonX_max(X, breaks = brks)
    checkEquals(res$y, xcms:::binYonX_max(X, X, nBins = 4L)$y)
    checkEquals(res$x, breakMidPoint(brks))
    ##  same breaks but sub-setting x
    res <- xcms:::binYonX_max(X, X, breaks = brks, fromIdx = 4, toIdx = 9)
    checkEquals(res$y, c(0, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))
    ##  breaks spanning a sub-set of x
    brks <- seq(3, 9, length.out = 5)
    res <- xcms:::binYonX_max(X, X, breaks = brks)
    checkEquals(res$y, c(4, 5, 7, 9))
    checkEquals(res$x, breakMidPoint(brks))
}


.test_inputs <- function() {

    library(xcms)

    xcms:::testIntegerInput(4)
    xcms:::testIntegerInput(c(4, 5))

    xcms:::testRealInput(4)
    xcms:::testRealInput(c(4.2231444, 5))

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
    b <- xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE)
    microbenchmark(xcms:::profBin(X, X, 5),
                   xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE))
    ## So, we're 3 times faster.
    microbenchmark(xcms:::binYonX_max(X, X, binSize = 1.65, sortedX = TRUE),
                   xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE))
}
