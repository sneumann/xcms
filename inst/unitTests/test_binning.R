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

## test_compare_with_profBin <- function() {

##     library(xcms)
##     library(RUnit)

##     set.seed(18011977)
##     X <- sort(abs(rnorm(500000, mean = 500, sd = 50)))
##     a <- xcms:::profBin(X, X, 1000)
##     b <- xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE,
##                             shiftHalfBinSize = TRUE)
## }

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

############################################################
## Keep that for now but delete eventually!
## Some old tests to evaluate binning using profBin function
.notrun_yet <- function() {

    library(xcms)

    set.seed(18011977)
    set.seed(123)
    mz <- sort(abs(rnorm(50000, mean = 200, sd = 300)))
    ints <- rnorm(50000, mean = 400, sd = 4000)
    ints[ints < 0] <- 0

    ## !!! values are supposed to be sorted!!!

    ## Dot a for loop in R...

    x <- 1:10
    nBins <- 5

    ## That's how I would do the stuff in c
    myProfBinR <- function(x, y, nBins, xstart = min(x), xend = max(x),
                           shiftBinCenter = 0) {
        dx <- (xend - xstart) / (nBins)
        ## breaks are actually identical to seq(1, 10, length.out = 6)
        breaks <- numeric(nBins + 1)
        currentBr <- xstart + shiftBinCenter * dx
        for(i in 1:length(breaks)){
            breaks[i] <- currentBr
            currentBr <- currentBr + dx
        }
        Res <- numeric(nBins)
        numX <- length(x)
        startx <- 0
        ## Determine the index from which we start...
        for (i in 1:length(x)) {
            if (x[i] >= xstart) {
                startx <- i
                break
            }
        }
        ## cat("starting in x from position ", startx, "\n")
        ## Eventually die if startx is still 0
        ## Loop through result, within each, loop through the
        ## x vector.
        for (i in 1:length(Res)) {
            upperBreak <- breaks[i+1]
            ## cat("i ", i, " upperBreak ", upperBreak)
            ## if (i == length(Res)){
            ##     while (x[startI] <= upperBreak) {
            ##         if (Res[i] < x[startI])
            ##             Res[i] <- x[startI]
            ##         startI <- startI + 1
            ##     }
            ## } else {
                while (startx <= numX & (x[startx] < upperBreak | x[startx] == upperBreak)) {
                    if (Res[i] < x[startx])
                        Res[i] <- x[startx]
                    ## cat(" startx", startx, " x[startx] ", x[startx], ";")
                    startx <- startx + 1
                }
            ## cat("\n")
            ## }
        }
        return(Res)
    }

    ## Other problem:
    ## have 1:10, want to bin them into 5 bins:
    ## 1,2 3,4 5,6 7,8 9,10
    xcms:::profBin(1:10, 1:10, num = 5)
    myProfBinR(1:10, 1:10, nBin = 5, shiftBinCenter = -0.5)
    (10 - 1) / (5-1)

    vals <- c(1, 2.2, 3, 4, 5, 6, 7, 8, 9, 10)
    xcms:::profBin(vals, vals, num = 5)
    myProfBinR(vals, vals, nBin = 5, shiftBinCenter = -0.5)

    vals <- c(1, 2.25, 3, 4, 5, 6, 7, 8, 9, 10)
    xcms:::profBin(vals, vals, num = 5)
    myProfBinR(vals, vals, nBin = 5, shiftBinCenter = -0.5)

    vals <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10)
    xcms:::profBin(vals, vals, num = 5)
    myProfBinR(vals, vals, nBin = 5, shiftBinCenter = -0.5)

    vals <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5, 11, 12, 13)
    xcms:::profBin(vals, vals, num = 5, xend = 10)
    myProfBinR(vals, vals, nBin = 5, shiftBinCenter = -0.5, xend = 10)
    ## is 11 not too large?
    (10 - 1) / (5-1)
    ## now; it's supposed to go from -0.125 to 11.125
    vals <- c(1, 2.05, 3, 4, 5, 6, 7, 8, 9, 10, 10.5, 11, 11.124, 11.125, 11.126, 12, 13)
    xcms:::profBin(vals, vals, num = 5, xend = 10)
    myProfBinR(vals, vals, nBin = 5, shiftBinCenter = -0.5, xend = 10)



    ## Other stuff below!

    ## Do the profBin in R:
    ProfBinR <- function(x, y, nBins, xstart = min(x), xend = max(x)) {
        dx <- (xend - xstart) - (nBins - 1)
        out <- numeric(nBins)
        startx <- xstart - dx
        endx <- xend + dx
        ## get the index from which we start...
        starti <- 0
        for (i in 1:length(x)) {
            if (x[i] >= startx) {
                starti <- i
                break
            }
        }
        cat("starting with i ", starti, "\n")
        for (i in starti:length(x)) {
            if (x[i] >= endx)
                break
            ## This doesn't work in R!
            outi <- floor((x[i] - xstart) / dx + 0.5)
            cat("i: ", i, " outi ", outi, " out[outi] ", out[outi], " x[i] ", x[i], "\n")
            if (outi >= 0 & outi < nBins) {
                if (out[outi] < y[i])
                    out[outi] <- y[i]
            }
        }
        return(out)
    }


    ## R binning. something wrong here!
    ## use x to define bins and bin y on them.
    ## Additional arguments: binSize, nBins and breaks
    binYonX_R <- function(x, y, binSize, nBins, breaks, FUN = max) {
        if (missing(breaks)) {
            ## binSize has precedence over nBins
            if (!missing(binSize)) {
                breaks <- seq(min(x), max(x), by = binSize)
            } else {
                if (!missing(nBins)) {
                    ## Add +1 to get the correct number of bins!
                    breaks <- seq(min(x), max(x), length.out = (nBins + 1))
                } else {
                    stop("Either 'binSize', 'nBins' or 'breaks' has to",
                         " be specified!")
                }
            }
        }

        ## Shift breaks by diff / 2
        breaks <- breaks - diff(breaks[1:2])/2
        nVals <- length(x)
        nbr <- length(breaks)
        ## Assign x to bins. all.inside = TRUE ensures all values are within.
        ## left.open: if true all intervals are open at left and closed at right.
        idx <- findInterval(x, breaks, all.inside = TRUE)

        out <- double(nbr - 1)
        ## Have to sort the index! split does also sorting.
        out[sort(unique(idx))] <- unlist(lapply(split(y, idx), FUN),
                                         use.names = FALSE)
        return(out)
        return(cbind(bin = (breaks[-nbr] + breaks[-1L])/2, binned = out))
    }

    ## Speed test.
    ##microbenchmark(xcms:::profBin(mz, ints, num = 100), binYonX_R(mz, ints, nBins = 100))

    brks <- seq(min(mz), max(mz), length.out = 101)

    Test <- tapply(ints, cut(mz, 100), max)

    inR <- binYonX_R(x = mz, y = ints, nBins = 100)

    ## Now let's see if I can do the same with xcms:
    inX <- xcms:::profBin(mz, ints, num = 100)  ## ?????
    Res <- cbind(inR, inX)

    ## NOT the same; eventually due to dx being -1?
    inX2 <- xcms:::profBin(mz, ints, num = 101)
    Res <- cbind(Res, inX2[1:100])

    inMy <- xcms:::binYonX_max(mz, ints, nBin = 100)
    Res <- cbind(Res, inMy)

    ## Why for heaven's sake are the numbers different????
    cbrs <- seq(min(mz), max(mz), length.out = 101)

    ## What is xcms doing???
    max(ints[mz > brs[1] & mz <= brs[2]])
    max(ints[mz >= brs[1] & mz < brs[2]])
    ## the same

    max(ints[mz > brs[2] & mz <= brs[3]])
    max(ints[mz >= brs[2] & mz < brs[3]])
    ## also the same!

    ## Very slow loop.
    for(i in 1:(length(brs)-1)){
        maxval <- ints[mz >= brs[i] & mz < brs[i+1]]
        if (all(is.na(maxval)))
            next
        expect_identical(unname(inR[i, 2]),
                         unname(max(maxval, na.rm=TRUE))
                         )
    }
    for(i in 1:(length(brs)-1)){
        maxval <- ints[mz > brs[i] & mz <= brs[i+1]]
        if (all(is.na(maxval)))
            next
        expect_identical(unname(inX[i]),
                         unname(max(maxval, na.rm=TRUE))
                         )
    }

    ## Trying to figure out what xcms is doing...
    dx <- (max(mz) - min(mz))/(100-1)
    ## This dx is however different than the diff of the breaks above!
    unique(diff(brs))

    inX2 <- xcms:::profBin(mz, ints, num = 101)
    length(inX2)
    nrow(inR)

    ## cbind them
    Res <- cbind(Res, inX2[1:100])
    head(Res)
    ## WHY???
    dx <- (max(mz) - min(mz))/(100)

    ## My implementation:
    inMy <- xcms:::binYonX_max(mz, ints, nBin = 100)


}


.benchmark_binning <- function() {
    ## Compare new binning functions with profBin.
    library(xcms)
    library(microbenchmark)
    X <- sort(abs(rnorm(500000, mean = 500, sd = 300)))
    a <- xcms:::profBin(X, X, 1000)
    b <- xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE)
    microbenchmark(xcms:::profBin(X, X, 5),
                   xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE))
    ## So, we're 3 times faster.
    microbenchmark(xcms:::binYonX_max(X, X, binSize = 1.65, sortedX = TRUE),
                   xcms:::binYonX_max(X, X, nBins = 1000, sortedX = TRUE))
}
