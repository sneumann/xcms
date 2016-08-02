## Functions to bin vectors.


## @description This functions takes two same-sized numeric vectors \code{x}
## and \code{y}, bins/cuts \code{x} into bins (either a pre-defined number
## of equal-sized bins or bins of a pre-defined size) and identifies for
## each bin the maximal value in \code{y} for the \code{x} values falling
## within each bin. \code{x} is expected to be incrementally sorted and, if
## not, it will be internally sorted (in which case also \code{y} will be
## ordered according to the order of \code{x}).
##
## @details The breaks for the bins can be defined with the input arguments
##  \code{nBins} or \code{binSize}. Pre-defined breaks can be passed with the
## argument \code{breaks}. If \code{nBins} are defined, the breaks are
## calculated similar to \code{brs <- seq(min(x), max(x), length.out = (nBins + 1))}).
## For each bin then the maximal \code{y} value
## (e.g. \code{max(y[x >= brs[1] & x < brs[2]])} for the first bin) is
## determined. See examples for more details.
## @note The function ensures that all values in the range defined by the breaks
## (including \code{binFromX} and \code{binToX}) are considered in the binning. This
## means that for all bins except the last one the condition is that values in
## \code{x} have to be \code{>= xlower} and \code{< xupper} (with \code{xlower}
## and \code{xupper} being the lower and upper boundary, respectively). For the
## last bin the condition is \code{x >= xlower & x <= xupper}.
## @title Get maximum values in y for bins defined on x
## @param x Numeric vector to be used for binning.
## @param y Numeric vector (same length than \code{x}) from which the maximum
## values for each bin should be defined. If not provided, \code{x} will
## be used.
## @param breaks Numeric vector defining the breaks for the bins, i.e. the
## lower and upper values for each bin. See examples below.
## @param nBins Integer of length one defining the number of desired bins.
## @param binSize Numeric of length one defining the desired bin size.
## @param binFromX Numeric of length one allowing to manually specify the range
## to be used for binning. This will affect only the calculation of the breaks
## for the bins (i.e. if \code{nBins} or \code{binSize} is provided).
## @param binToX Same as \code{binFromX}, but defining the maximum x-value to be
## used for binning.
## @param fromIdx Integer of length 1 allowing to define the sub-set of x
## that will be used for binning.
## @param toIdx Same as \code{toIdx}, but defining the maximum index in x to
## be used for binning.
## @param sortedX Whether \code{x} is sorted.
## @return Returns a numeric vector with the maximal y-values for each bin
## of x-values.
## @author Johannes Rainer
## @examples
## ########
## ## Simple example illustrating the breaks and the binning.
## ##
## ## Define breaks for 5 bins:
## brks <- seq(2, 12, length.out = 6)
## ## The first bin is then [2,4), the second [4,6) and so on.
## brks
## ## Get the max value falling within each bin.
## xcms:::binYonX_max(x = 1:16, y = 1:16, breaks = brks)
## ## Thus, the largest value in x = 1:16 falling into the bin [2,4) (i.e. being
## ## >= 2 and < 4) is 3, the largest one falling into [4,6) is 5 and so on.
## ## Note however the function ensures that binFromX and binToX (in this example 1 and
## ## 12) fall within a bin, i.e. 12 is considered for the last bin.
##
## #######
## ## Bin only values within a sub-set of x
## ##
## ## This example illustrates how the fromIdx and toIdx parameters can be used.
## ## x defines 3 times the sequence form 1 to 10, while y is the sequence from
## ## 1 to 30. In this very simple example x is supposed to represent M/Z values
## ## from 3 consecutive scans and y the intensities measured for each M/Z in
## ## each scan. We want to get the maximum intensities for M/Z value bins only
## ## for the second scan, and thus we use fromIdx = 11 and toIdx = 20. The breaks
## ## for the bins are defined with the nBins, binFromX and binToX.
## X <- rep(1:10, 3)
## Y <- 1:30
## ## Bin the M/Z values in the second scan into 5 bins and get the maximum
## ## intensity for each bin. Note that we have to specify sortedX = TRUE as
## ## the x and y vectors would be sorted otherwise.
## xcms:::binYonX_max(X, Y, nBins = 5L, binFromX = 1, binToX = 10,
##                    sortedX = TRUE, fromIdx = 11, toIdx = 20)
binYonX_max <- function(x, y = x, breaks, nBins, binSize, binFromX = min(x),
                        binToX = max(x), fromIdx = 1L, toIdx = length(x),
                        sortedX = !is.unsorted(x)) {
    if(missing(x) | missing(y))
        stop("Arguments 'x' and 'y' are mandatory!")
    if(missing(nBins) & missing(binSize) & missing(breaks))
        stop("One of 'breaks', 'nBins' or 'binSize' has to be defined!")
    if (!sortedX) {
        message("'x' is not sorted, will sort 'x' and 'y'.")
        o <- order(x, method = "radix")
        x <- x[o]
        y <- y[o]
    }
    ## Check that input values have the correct data types.
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    if (!is.double(binFromX)) binFromX <- as.double(binFromX)
    if (!is.double(binToX)) binToX <- as.double(binToX)
    if (!is.integer(fromIdx)) as.integer(fromIdx)
    if (!is.integer(toIdx)) as.integer(toIdx)
    ## Call the C function(s). breaks has precedence over nBins over binSize.
    if (!missing(breaks)) {
        if (!is.double(breaks)) breaks <- as.double(nBins)
        return(.Call("binXonY_breaks_max", x, y, breaks, force(fromIdx - 1L),
                     force(toIdx - 1L), PACKAGE = "xcms"))
    }
    if (!missing(nBins)) {
        if (!is.integer(nBins)) nBins <- as.integer(nBins)
        return(.Call("binXonY_nBins_max", x, y, nBins, binFromX, binToX,
                     force(fromIdx - 1L), force(toIdx - 1L), PACKAGE = "xcms"))
    }
    if (!missing(binSize)) {
        if (!is.double(binSize)) binSize <- as.double(binSize)
        return(.Call("binXonY_binSize_max", x, y, binSize, binFromX, binToX,
                     force(fromIdx - 1L), force(toIdx - 1L), PACKAGE = "xcms"))
    }
}

## @description Calculate breaks for same-sized bins for data values
## from \code{fromX} to \code{toX}.
##
## @details This generates bins such as a call to
## \code{seq(fromX, toX, length.out = nBins)} would. The first and second element
## in the result vector thus defines the lower and upper boundary for the first
## bin, the second and third value for the second bin and so on.
## @title Generate breaks for binning
## @param fromX Numeric of length 1 specifying the lowest value for the bins.
## @param toX Numeric of length 1 specifying the largest value for the bins.
## @param nBins Integer of length 1 defining the number of bins.
## @return A numeric vector of length \code{nBins + 1} defining the lower and
## upper bounds of the bins.
## @author Johannes Rainer
breaks_on_nBins <- function(fromX, toX, nBins) {
    if(missing(fromX) | missing(toX) | missing(nBins))
        stop("'fromX', 'toX' and 'nBins' are required!")
    if (!is.double(fromX)) fromX <- as.double(fromX)
    if (!is.double(toX)) toX <- as.double(toX)
    if (!is.integer(nBins)) nBins <- as.integer(nBins)
    return(.Call("breaks_on_nBins", fromX, toX, nBins, PACKAGE = "xcms"))
}

## @description Defines breaks for \code{binSize} sized bins for values ranging
## from \code{fromX} to \code{toX}.
##
## @details This generates breaks for bins such that
## \code{seq(fromX, toX, by = binSize)} would, with the exception that one last break
## is added to ensure that also \code{toX} would be included (see examples). This
## function can be up to twice as fast as the \code{seq} call.
## @title Generate breaks for binning using a defined bin size.
## @param fromX Numeric of length 1 specifying the lowest value for the bins.
## @param toX Numeric of length 1 specifying the largest value for the bins.
## @param binSize Numeric of length 1 defining the size of a bin.
## @return A numeric vector defining the lower and upper bounds of the bins.
## @author Johannes Rainer
##' @examples
##'
breaks_on_binSize <- function(fromX, toX, binSize) {
    if(missing(fromX) | missing(toX) | missing(binSize))
        stop("'fromX', 'toX' and 'binSize' are required!")
    if (!is.double(fromX)) fromX <- as.double(fromX)
    if (!is.double(toX)) toX <- as.double(toX)
    if (!is.integer(binSize)) binSize <- as.double(binSize)
    return(.Call("breaks_on_binSize", fromX, toX, binSize, PACKAGE = "xcms"))
}



############################################################
## Some other tests.
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

