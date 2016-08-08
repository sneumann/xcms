## Functions to bin vectors.


## @description This functions takes two same-sized numeric vectors \code{x}
## and \code{y}, bins/cuts \code{x} into bins (either a pre-defined number
## of equal-sized bins or bins of a pre-defined size) and aggregates values
## in \code{y} corresponding to \code{x} values falling within each bin. By
## default (i.e. \code{method = "max"}) the maximal \code{y} value for the
## corresponding \code{x} values is identified. \code{x} is expected to be
## incrementally sorted and, if not, it will be internally sorted (in which
## case also \code{y} will be ordered according to the order of \code{x}).
##
## @details The breaks for the bins can be defined with the input arguments
##  \code{nBins} or \code{binSize}. Pre-defined breaks can be passed with the
## argument \code{breaks}. If \code{nBins} are defined, the breaks are
## calculated similar to \code{brs <- seq(min(x), max(x), length.out = (nBins + 1))}).
## For each bin then the maximal \code{y} value
## (e.g. \code{max(y[x >= brs[1] & x < brs[2]])} for the first bin) is
## determined. See examples for more details.
##
## Imputation of missing bin values: by default (argument \code{impute = "no"})
## the value of \code{missingValue} is reported for a bin if none of the values
## in \code{x} fall into the bin. Alternatively, a value for such bins can be
## calculated by specifying an interpolation method with argument \code{impute}.
## \code{impute = "lin"} uses simple linear imputation to derive a value for a
## bin with missing values from its neighboring bins (for which a value was
## defined). This method is equivalent to the \code{\link{profBinLin}} method.
## For more details see examples below.
## @note The function ensures that all values in the range defined by the breaks
## (including \code{binFromX} and \code{binToX}) are considered in the binning. This
## means that for all bins except the last one values in \code{x} have to be
## \code{>= xlower} and \code{< xupper} (with \code{xlower}
## and \code{xupper} being the lower and upper boundary, respectively). For the
## last bin the condition is \code{x >= xlower & x <= xupper}.
## Note also that if \code{shiftByHalfBinSize} is \code{TRUE} the range of values
## that is used for binning is expanded by \code{binSize} (i.e. the lower boundary
## will be \code{fromX - binSize/2}, the upper \code{toX + binSize/2}). Setting
## this argument to \code{TRUE} resembles the binning that is/was used in
## \code{xcms} \code{\link{profBin}} and similar methods.
##
## \code{NA} handling: by default the function ignores \code{NA} values in
## \code{y} (thus inherently assumes \code{na.rm = TRUE}). No \code{NA} values are
## allowed in \code{x}.
##
## All methods to impute missing bin values are based on the already aggregated
## data, not the individual values in \code{y}.
##
## @title Aggregate values in y for bins defined on x
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
## @method A character string specifying the method that should be used to
## aggregate values in \code{y}. Allowed are \code{"max"}, \code{"min"},
## \code{"sum"} and \code{"mean"} to identify the maximal or minimal value or
## to sum all values within a bin or calculate their mean value.
## @param sortedX Whether \code{x} is sorted.
## @param shiftByHalfBinSize Logical specifying whether the bins should be
## shifted by half the bin size to the left. Thus, the first bin will have its
## center at \code{fromX} and its lower and upper boundary are
##  \code{fromX - binSize/2} and \code{fromX + binSize/2}. This argument is
## ignored if \code{breaks} are provided.
## @missingValue A numeric of length one specifying the value that should be
## used as \emph{default} value for a bin. Defaults to \code{0}, thus if no value
## in \code{x} falls within a bin, \code{0} is reported for this bin.
## @impute Allows to specify a method to be used for imputation of missing
## bin-values. Allowed values are \code{"no"} (no imputation), \code{"lin"}
## (linear interpolation). See details section for more information.
## @return Returns a list of length 2, the first element (named \code{"x"})
## contains the bin mid-points, the second element (named \code{"y"}) the
## aggregated values from input vector \code{y} within each bin.
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
##
## #######
## ## Use linear interpolation to impute missing bin values.
## ##
## ## The input vector y contains NAs for x-values falling within some bins. By
## ## default the value of argument missingValue is reported for these.
## X <- 1:11
## Y <- c(1, NA, NA, NA, 5, 6, NA, NA, 9, 10, 11)
## ## The breaks for the bins are:
## xcms:::breaks_on_nBins(1, 11, 5L)
## ## Thus bin 2 and 4 will get a value of 0.
## binYonX(X, Y, nBins = 5L)$y
## ## If we use impute = "lin" a value is inferred for these bins based on their
## ## neighboring bins, i.e. the value for the second bin will be (6 + 1)/2 and
## ## the value for the 4th bin will be (6 + 11)/2.
## binYonX(X, Y, nBins = 5L, impute = "lin")$y
binYonX <- function(x, y, breaks, nBins, binSize, binFromX = min(x),
                    binToX = max(x), fromIdx = 1L, toIdx = length(x),
                    method = "max",
                    sortedX = !is.unsorted(x), shiftByHalfBinSize = FALSE,
                    impute = "no",
                    missingValue = 0) {
    if (!missing(x) & missing(y))
        y <- x
    if (missing(x) | missing(y))
        stop("Arguments 'x' and 'y' are mandatory!")
    if (missing(nBins) & missing(binSize) & missing(breaks))
        stop("One of 'breaks', 'nBins' or 'binSize' has to be defined!")
    if (!sortedX) {
        message("'x' is not sorted, will sort 'x' and 'y'.")
        o <- order(x, method = "radix")
        x <- x[o]
        y <- y[o]
    }
    ## For now we don't allow NAs in x
    if (anyNA(x))
        stop("No 'NA' values are allowed in 'x'!")
    ## Check that input values have the correct data types.
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    if (!is.double(binFromX)) binFromX <- as.double(binFromX)
    if (!is.double(binToX)) binToX <- as.double(binToX)
    if (!is.integer(fromIdx)) as.integer(fromIdx)
    if (!is.integer(toIdx)) as.integer(toIdx)
    ## Call the C function(s). breaks has precedence over nBins over binSize.
    shiftIt <- 0L
    if (!missing(breaks)) {
        if (shiftByHalfBinSize)
            warning("Argument 'shiftByHalfBinSize' is ignored if 'breaks'",
                    " are provided.")
        if (!is.double(breaks)) breaks <- as.double(nBins)
        nBins <- NA_integer_
        binSize <- as.double(NA)
    } else {
        if (!missing(nBins)) {
            if (shiftByHalfBinSize)
                shiftIt <- 1L
            breaks <- as.double(NA)
            if (!is.integer(nBins)) nBins <- as.integer(nBins)
            binSize <- as.double(NA)
        } else{
            if (shiftByHalfBinSize)
                shiftIt <- 1L
            breaks <- as.double(NA)
            nBins <- NA_integer_
            if (!is.double(binSize)) binSize <- as.double(binSize)
        }
    }
    return(.Call("binYonX", x, y, breaks, nBins, binSize,
                 binFromX, binToX, force(fromIdx - 1L), force(toIdx - 1L),
                 shiftIt, as.double(missingValue),
                 as.integer(.aggregateMethod2int(method)),
                 as.integer(.imputeMethod2int(impute)),
                 PACKAGE = "xcms"))

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
breaks_on_nBins <- function(fromX, toX, nBins, shiftByHalfBinSize = FALSE) {
    if(missing(fromX) | missing(toX) | missing(nBins))
        stop("'fromX', 'toX' and 'nBins' are required!")
    if (!is.double(fromX)) fromX <- as.double(fromX)
    if (!is.double(toX)) toX <- as.double(toX)
    if (!is.integer(nBins)) nBins <- as.integer(nBins)
    shiftIt <- 0L
    if (shiftByHalfBinSize)
        shiftIt <- 1L
    return(.Call("breaks_on_nBins", fromX, toX, nBins, shiftIt,
                 PACKAGE = "xcms"))
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
## @examples
##
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

############################################################
## .aggregateMethod2int
## Converts the method character string to an integer that can
## be used by the downstream C-functions to switch method.
.aggregateMethod2int <- function(method = "max") {
    method <- match.arg(method, names(.aggregateMethods))
    return(.aggregateMethods[method])
}
.aggregateMethods <- c(1, 2, 3, 4)
names(.aggregateMethods) <- c("max", "min", "sum", "mean")


############################################################
## .imputeMethod2int
## Converts the imputation method into an integer that can be
## passed down to the C-function.
.imputeMethod2int <- function(method = "no") {
    method <- match.arg(method, names(.imputeMethods))
    return(.imputeMethods[method])
}
.imputeMethods <- c(0, 1)
names(.imputeMethods) <- c("no", "lin")
