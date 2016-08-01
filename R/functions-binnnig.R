## Functions to bin vectors.


## @description This function bins \code{x} into \code{nBin} equal-spaced bins
## and identifies the maximal \code{y}-value among the subset of \code{y} values
## corresponding to (having the same index as) values in \code{x} falling into
## each bin.
## @details The function first defines the bins/intervals on \code{x}
## (similar to \code{brs <- seq(min(x), max(x), length.out = (nBin + 1))}) and
## determines then for each bin the maximal \code{y} value (e.g.
## \code{max(y[x >= brs[1] & x < brs[2]])} for the first bin).
## @note \code{x} is supposed to be ordered increasingly. If it is not,
## \code{x} will be ordered internally (along with \code{y}).
## @title Get the maximal y value for bins on x
## @param x Numeric vector on which the bins should be defined. \code{x} is
## supposed to be ordered.
## @param y Numeric vector of values from which the max value per bin should
## be returned.
## @param nBin The number of bins in which \code{x} should be binned.
## @param fromX Optionally calculate the binning for a subset of \code{x} starting
## at \code{fromX}.
## @param toX Optionally calculate the binning for a subset of \code{x} up
## to \code{toX}.
## @param param
## @param sortedX Logical indicating whether the values in \code{x} are ordered.
## @return
## @author Johannes Rainer
binYonX_max <- function(x, y, nBins, fromX = min(x), toX = max(x),
                        param = list(), sortedX = !is.unsorted(x)) {
    ## get function name from function: deparse(quote(max))
    if (!sortedX) {
        o <- order(x, method = "radix")
        x <- x[o]
        y <- y[o]
    }
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    res <- .C("_binYonX_max",
              x,
              y,
              as.integer(length(x)),
              as.double(fromX),
              as.double(toX),
              as.integer(nBins),
              out = double(nBins),
              PACKAGE = "xcms")$out
    return(res)
}

binTest <- function(x, y, nBins, binSize, fromX = min(x), toX = max(x),
                    sortedX = !is.unsorted(x)) {
    if(missing(x) | missing(y))
        stop("Arguments 'x' and 'y' are mandatory!")
    if(missing(nBins) & missing(binSize))
        stop("One of 'nBins' or 'binSize' has to be defined!")
    if (!sortedX) {
        o <- order(x, method = "radix")
        x <- x[o]
        y <- y[o]
    }
    ## Check that input values have the correct data types.
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    if (!is.double(fromX)) fromX <- as.double(fromX)
    if (!is.double(toX)) toX <- as.double(toX)
    ## Cal the C function.
    if(!missing(nBins)) {
        if (!is.integer(nBins)) nBins <- as.integer(nBins)
        return(.Call("binXonY_nBins_max", x, y, nBins, fromX, toX,
                     PACKAGE = "xcms"))
    }
    if(!missing(binSize)) {
        if (!is.double(binSize)) binSize <- as.double(binSize)
        return(.Call("binXonY_binSize_max", x, y, binSize, fromX, toX,
                     PACKAGE = "xcms"))
        }
}



testIntegerInput <- function(x) {
    if (!is.integer(x)) x <- as.integer(x)
    res <- .Call("test_integer", x, PACKAGE = "xcms")
    return(res)
}

testRealInput <- function(x) {
    if (!is.double(x)) x <- as.double(x)
    res <- .Call("test_real", x, PACKAGE = "xcms")
    return(res)
}


## binYonX_max_with_breaks <- function(x, y, breaks, sortedX = !is.unsorted(x)) {
##     if (!sortedX) {
##         o <- order(x, method = "radix")
##         x <- x[o]
##         y <- y[o]
##     }
##     if (!is.double(x)) x <- as.double(x)
##     if (!is.double(y)) y <- as.double(y)
## }

