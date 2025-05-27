## Functions to bin vectors.

#' @title Aggregate values in y for bins defined on x
#'
#' @description This functions takes two same-sized numeric vectors `x`
#'     and `y`, bins/cuts `x` into bins (either a pre-defined number
#'     of equal-sized bins or bins of a pre-defined size) and aggregates values
#'     in `y` corresponding to `x` values falling within each bin. By
#'     default (i.e. `method = "max"`) the maximal `y` value for the
#'     corresponding `x` values is identified. `x` is expected to be
#'     incrementally sorted and, if not, it will be internally sorted (in which
#'     case also `y` will be ordered according to the order of `x`).
#'
#' @details The breaks defining the boundary of each bin can be either passed
#'     directly to the function with the argument `breaks`, or are
#'     calculated on the data based on arguments `nBins` or `binSize`
#'     along with `fromIdx`, `toIdx` and optionally `binFromX`
#'     and `binToX`.
#'     Arguments `fromIdx` and `toIdx` allow to specify subset(s) of
#'     the input vector `x` on which bins should be calculated. The
#'     default the full `x` vector is considered. Also, if not specified
#'     otherwise with arguments `binFromX` and `binToX`, the range
#'     of the bins within each of the sub-sets will be from `x[fromIdx]`
#'     to `x[toIdx]`. Arguments `binFromX` and `binToX` allow to
#'     overwrite this by manually defining the a range on which the breaks
#'     should be calculated. See examples below for more details.
#'
#'     Calculation of breaks: for `nBins` the breaks correspond to
#'     `seq(min(x[fromIdx])), max(x[fromIdx], length.out = (nBins + 1))`.
#'     For `binSize` the breaks correspond to
#'     `seq(min(x[fromIdx]), max(x[toIdx]), by = binSize)` with the
#'     exception that the last break value is forced to be equal to
#'     `max(x[toIdx])`. This ensures that all values from the specified
#'     range are covered by the breaks defining the bins. The last bin could
#'     however in some instances be slightly larger than `binSize`. See
#'     [breaks_on_binSize()] and [breaks_on_nBins()] for
#'     more details.
#'
#' @note The function ensures that all values within the range used to define
#'     the breaks are considered in the binning (and assigned to a bin). This
#'     means that for all bins except the last one values in `x` have to be
#'     `>= xlower` and `< xupper` (with `xlower`
#'     and `xupper` being the lower and upper boundary, respectively). For
#'     the last bin the condition is `x >= xlower & x <= xupper`.
#'     Note also that if `shiftByHalfBinSize` is `TRUE` the range of
#'     values that is used for binning is expanded by `binSize` (i.e. the
#'     lower boundary will be `fromX - binSize/2`, the upper
#'     `toX + binSize/2`). Setting this argument to `TRUE` resembles
#'     the binning that is/was used in `profBin` function from
#'     *xcms* < 1.51.
#'
#'     `NA` handling: by default the function ignores `NA` values in
#'     `y` (thus inherently assumes `na.rm = TRUE`). No `NA`
#'     values are allowed in `x`.
#'
#' @param x Numeric vector to be used for binning.
#'
#' @param y Numeric vector (same length than `x`) from which the maximum
#'     values for each bin should be defined. If not provided, `x` will
#'     be used.
#'
#' @param breaks Numeric vector defining the breaks for the bins, i.e. the
#'     lower and upper values for each bin. See examples below.
#'
#' @param nBins `integer(1)` defining the number of desired bins.
#'
#' @param binSize `numeric(1)` defining the desired bin size.
#'
#' @param binFromX Optional `numeric(1)` allowing to manually specify
#'     the range of x-values to be used for binning.
#'     This will affect only the calculation of the breaks for the bins
#'     (i.e. if `nBins` or `binSize` is provided).
#'     If not provided the minimal value in the sub-set
#'     `fromIdx`-`toIdx` in input vector `x` will be used.
#'
#' @param binToX Same as `binFromX`, but defining the maximum x-value to be
#'     used for binning.
#'
#' @param fromIdx Integer vector defining the start position of one or multiple
#'     sub-sets of input vector `x` that should be used for binning.
#'
#' @param toIdx Same as `toIdx`, but defining the maximum index (or indices)
#'     in x to be used for binning.
#'
#' @param method A character string specifying the method that should be used to
#'     aggregate values in `y`. Allowed are `"max"`, `"min"`,
#'     `"sum"` and `"mean"` to identify the maximal or minimal value
#'     or to sum all values within a bin or calculate their mean value.
#'
#' @param sortedX Whether `x` is sorted.
#'
#' @param shiftByHalfBinSize Logical specifying whether the bins should be
#'     shifted by half the bin size to the left. Thus, the first bin will have
#'     its center at `fromX` and its lower and upper boundary are
#'     `fromX - binSize/2` and `fromX + binSize/2`. This argument is
#'     ignored if `breaks` are provided.
#'
#' @param baseValue The base value for empty bins (i.e. bins into which either
#'     no values in `x` did fall, or to which only `NA` values in
#'     `y` were assigned). By default (i.e. if not specified), `NA`
#'     is assigned to such bins.
#'
#' @param returnIndex Logical indicating whether the index of the max (if
#'     `method = "max"`) or min (if `method = "min"`) value within
#'     each bin in input vector `x` should also be reported. For methods
#'     other than `"max"` or `"min"` this argument is ignored.
#'
#' @param returnX `logical` allowing to avoid returning `$x`, i.e. the
#'     mid-points of the bins. `returnX = FALSE` might be useful in cases
#'     where `breaks` are pre-defined as it considerably reduces the memory
#'     demand.
#'
#' @return Returns a list of length 2, the first element (named `"x"`)
#'     contains the bin mid-points, the second element (named `"y"`) the
#'     aggregated values from input vector `y` within each bin. For
#'     `returnIndex = TRUE` the list contains an additional element
#'     `"index"` with the index of the max or min (depending on whether
#'     `method = "max"` or `method = "min"`) value within each bin in
#'     input vector `x`.
#'
#' @author Johannes Rainer
#'
#' @seealso [imputeLinInterpol()]
#'
#' @md
#'
#' @examples
#' ########
#' ## Simple example illustrating the breaks and the binning.
#' ##
#' ## Define breaks for 5 bins:
#' brks <- seq(2, 12, length.out = 6)
#' ## The first bin is then [2,4), the second [4,6) and so on.
#' brks
#' ## Get the max value falling within each bin.
#' binYonX(x = 1:16, y = 1:16, breaks = brks)
#' ## Thus, the largest value in x = 1:16 falling into the bin [2,4) (i.e. being
#' ## >= 2 and < 4) is 3, the largest one falling into [4,6) is 5 and so on.
#' ## Note however the function ensures that the minimal and maximal x-value
#' ## (in this example 1 and 12) fall within a bin, i.e. 12 is considered for
#' ## the last bin.
#'
#' #######
#' ## Performing the binning ons sub-set of x
#' ##
#' X <- 1:16
#' ## Bin X from element 4 to 10 into 5 bins.
#' X[4:10]
#' binYonX(X, X, nBins = 5L, fromIdx = 4, toIdx = 10)
#' ## This defines breaks for 5 bins on the values from 4 to 10 and bins
#' ## the values into these 5 bins. Alternatively, we could manually specify
#' ## the range for the binning, i.e. the minimal and maximal value for the
#' ## breaks:
#' binYonX(X, X, nBins = 5L, fromIdx = 4, toIdx = 10, binFromX = 1, binToX = 16)
#' ## In this case the breaks for 5 bins were defined from a value 1 to 16 and
#' ## the values 4 to 10 were binned based on these breaks.
#'
#' #######
#' ## Bin values within a sub-set of x, second example
#' ##
#' ## This example illustrates how the fromIdx and toIdx parameters can be used.
#' ## x defines 3 times the sequence form 1 to 10, while y is the sequence from
#' ## 1 to 30. In this very simple example x is supposed to represent M/Z values
#' ## from 3 consecutive scans and y the intensities measured for each M/Z in
#' ## each scan. We want to get the maximum intensities for M/Z value bins only
#' ## for the second scan, and thus we use fromIdx = 11 and toIdx = 20. The breaks
#' ## for the bins are defined with the nBins, binFromX and binToX.
#' X <- rep(1:10, 3)
#' Y <- 1:30
#' ## Bin the M/Z values in the second scan into 5 bins and get the maximum
#' ## intensity for each bin. Note that we have to specify sortedX = TRUE as
#' ## the x and y vectors would be sorted otherwise.
#' binYonX(X, Y, nBins = 5L, sortedX = TRUE, fromIdx = 11, toIdx = 20)
#'
#' #######
#' ## Bin in overlapping sub-sets of X
#' ##
#' ## In this example we define overlapping sub-sets of X and perform the binning
#' ## within these.
#' X <- 1:30
#' ## Define the start and end indices of the sub-sets.
#' fIdx <- c(2, 8, 21)
#' tIdx <- c(10, 25, 30)
#' binYonX(X, nBins = 5L, fromIdx = fIdx, toIdx = tIdx)
#' ## The same, but pre-defining also the desired range of the bins.
#' binYonX(X, nBins = 5L, fromIdx = fIdx, toIdx = tIdx, binFromX = 4, binToX = 28)
#' ## The same bins are thus used for each sub-set.
binYonX <- function(x, y, breaks, nBins, binSize, binFromX,
                    binToX, fromIdx = 1L, toIdx = length(x),
                    method = "max", baseValue,
                    sortedX = !is.unsorted(x),
                    shiftByHalfBinSize = FALSE,
                    returnIndex = FALSE, returnX = TRUE) {
    if (!missing(x) & missing(y))
        y <- x
    if (missing(x) | missing(y))
        stop("Arguments 'x' and 'y' are mandatory!")
    if (missing(nBins) & missing(binSize) & missing(breaks))
        stop("One of 'breaks', 'nBins' or 'binSize' has to be defined!")
    if (!sortedX) {
        message("'x' is not sorted, will sort 'x' and 'y'.")
        ## Sort method; see issue #180 for MSnbase
        ## Note: order method = "radix" is considerably faster - but there is no
        ## method argument for older R versions.
        o <- order(x)
        x <- x[o]
        y <- y[o]
    }
    ## Check fromIdx and toIdx
    if (any(fromIdx < 1) | any(toIdx > length(x)))
        stop("'fromIdx' and 'toIdx' have to be within 1 and lenght(x)!")
    if (length(toIdx) != length(fromIdx))
        stop("'fromIdx' and 'toIdx' have to have the same length!")
    if (missing(binFromX))
        binFromX <- as.double(NA)
    if (missing(binToX))
        binToX <- as.double(NA)
    ## For now we don't allow NAs in x
    if (anyNA(x))
        stop("No 'NA' values are allowed in 'x'!")
    ## Check that input values have the correct data types.
    if (!is.double(x)) x <- as.double(x)
    if (!is.double(y)) y <- as.double(y)
    if (!is.double(binFromX)) binFromX <- as.double(binFromX)
    if (!is.double(binToX)) binToX <- as.double(binToX)
    if (!is.integer(fromIdx)) fromIdx <- as.integer(fromIdx)
    if (!is.integer(toIdx)) toIdx <- as.integer(toIdx)
    ## breaks has precedence over nBins over binSize.
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
            breaks <- as.double(NA)
            if (!is.integer(nBins)) nBins <- as.integer(nBins)
            binSize <- as.double(NA)
        } else{
            breaks <- as.double(NA)
            nBins <- NA_integer_
            if (!is.double(binSize)) binSize <- as.double(binSize)
        }
    }
    if (shiftByHalfBinSize)
        shiftIt <- 1L
    ## Define default value for baseValue
    if (missing(baseValue)) {
        baseValue = as.double(NA)
    } else {
        if (!is.double(baseValue)) baseValue <- as.double(baseValue)
    }

    getIndex <- 0L
    if (returnIndex)
        getIndex <- 1L
    getX <- 0L
    if (returnX)
        getX <- 1L
    if (length(toIdx) > 1) {
        .Call("binYonX_multi", x, y, breaks, nBins, binSize,
              binFromX, binToX, force(fromIdx - 1L), force(toIdx - 1L),
              shiftIt,
              as.integer(.aggregateMethod2int(method)),
              baseValue,
              getIndex,
              getX,
              PACKAGE = "xcms")
    } else {
        .Call("binYonX", x, y, breaks, nBins, binSize,
              binFromX, binToX, force(fromIdx - 1L), force(toIdx - 1L),
              shiftIt,
              as.integer(.aggregateMethod2int(method)),
              baseValue,
              getIndex,
              getX,
              PACKAGE = "xcms")
    }
}


############################################################
## imputeLinInterpol
##
#' @title Impute values for empty elements in a vector using linear interpolation
#'
#' @description This function provides missing value imputation based on linear
#'     interpolation and resembles some of the functionality of the
#'     `profBinLin()` and `profBinLinBase()` functions deprecated from
#'     version 1.51 on.
#'
#' @details Values for NAs in input vector `x` can be imputed using methods
#'     `"lin"` and `"linbase"`:
#'
#'     `impute = "lin"` uses simple linear imputation to derive a value
#'     for an empty element in input vector `x` from its neighboring
#'     non-empty elements. This method is equivalent to the linear
#'     interpolation in the `profBinLin` method. Whether interpolation is
#'     performed if missing values are present at the beginning and end of
#'     `x` can be set with argument `noInterpolAtEnds`. By default
#'     interpolation is also performed at the ends interpolating from `0`
#'     at the beginning and towards `0` at the end. For
#'     `noInterpolAtEnds = TRUE` no interpolation is performed at both
#'     ends replacing the missing values at the beginning and/or the end of
#'     `x` with `0`.
#'
#'     `impute = "linbase"` uses linear interpolation to impute values for
#'     empty elements within a user-definable proximity to non-empty elements
#'     and setting the element's value to the `baseValue` otherwise. The
#'     default for the `baseValue` is half of the smallest value in
#'     `x` (`NA`s being removed). Whether linear interpolation based
#'     imputation is performed for a missing value depends on the
#'     `distance` argument. Interpolation is only performed if one of the
#'     next `distance` closest neighbors to the current empty element has
#'     a value other than `NA`. No interpolation takes place for
#'     `distance = 0`, while `distance = 1` means that the value for
#'     an empty element is interpolated from directly adjacent non-empty
#'     elements while, if the next neighbors of the current empty element are
#'     also `NA`, it's vale is set to `baseValue`.
#'     This corresponds to the linear interpolation performed by the
#'     `profBinLinBase` method. For more details see examples below.
#'
#' @param x A numeric vector with eventual missing (`NA`) values.
#'
#' @param baseValue The base value to which empty elements should be set. This
#'     is only considered for `method = "linbase"` and corresponds to the
#'     `profBinLinBase()`'s `baselevel` argument.
#'
#' @param method One of `"none"`, `"lin"` or `"linbase"`.
#'
#' @param distance For `method = "linbase"`: number of non-empty
#'     neighboring element of an empty element that should be considered for
#'     linear interpolation. See details section for more information.
#'
#' @param noInterpolAtEnds For `method = "lin"`: Logical indicating
#'     whether linear interpolation should also be performed at the ends of the
#'     data vector (i.e. if missing values are present at the beginning or the
#'     end of the vector).
#'
#' @return A numeric vector with empty values imputed based on the selected
#'     `method`.
#'
#' @md
#'
#' @author Johannes Rainer
#'
#' @examples
#' #######
#' ## Impute missing values by linearly interpolating from neighboring
#' ## non-empty elements
#' x <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
#' imputeLinInterpol(x, method = "lin")
#' ## visualize the interpolation:
#' plot(x = 1:length(x), y = x)
#' points(x = 1:length(x), y = imputeLinInterpol(x, method = "lin"), type = "l", col = "grey")
#'
#' ## If the first or last elements are NA, interpolation is performed from 0
#' ## to the first non-empty element.
#' x <- c(NA, 2, 1, 4, NA)
#' imputeLinInterpol(x, method = "lin")
#' ## visualize the interpolation:
#' plot(x = 1:length(x), y = x)
#' points(x = 1:length(x), y = imputeLinInterpol(x, method = "lin"), type = "l", col = "grey")
#'
#' ## If noInterpolAtEnds is TRUE no interpolation is performed at both ends
#' imputeLinInterpol(x, method = "lin", noInterpolAtEnds = TRUE)
#'
#' ######
#' ## method = "linbase"
#' ## "linbase" performs imputation by interpolation for empty elements based on
#' ## 'distance' adjacent non-empty elements, setting all remaining empty elements
#' ## to the baseValue
#' x <- c(3, NA, 1, 2, NA, NA, 4, NA, NA, NA, 3, NA, NA, NA, NA, 2)
#' ## Setting distance = 0 skips imputation by linear interpolation
#' imputeLinInterpol(x, method = "linbase", distance = 0)
#'
#' ## With distance = 1 for all empty elements next to a non-empty element the value
#' ## is imputed by linear interpolation.
#' xInt <- imputeLinInterpol(x, method = "linbase", distance = 1L)
#' xInt
#'
#' plot(x = 1:length(x), y = x, ylim = c(0, max(x, na.rm = TRUE)))
#' points(x = 1:length(x), y = xInt, type = "l", col = "grey")
#'
#' ## Setting distance = 2L would cause that for all empty elements for which the
#' ## distance to the next non-empty element is <= 2 the value is imputed by
#' ## linear interpolation:
#' xInt <- imputeLinInterpol(x, method = "linbase", distance = 2L)
#' xInt
#'
#' plot(x = 1:length(x), y = x, ylim = c(0, max(x, na.rm = TRUE)))
#' points(x = 1:length(x), y = xInt, type = "l", col = "grey")
imputeLinInterpol <- function(x, baseValue, method = "lin", distance = 1L,
                              noInterpolAtEnds = FALSE) {
    method <- match.arg(method, c("none", "lin", "linbase")) ## interDist? distance = 2
    if (method == "none") {
        return(x)
    }
    if (!is.double(x)) x <- as.double(x)
    if (method == "lin") {
        noInter <- 0L
        if (noInterpolAtEnds)
            noInter <- 1L
        return(.Call("impute_with_linear_interpolation", x, noInter,
                     PACKAGE = "xcms"))
    }
    if (method == "linbase") {
        if (missing(baseValue))
            baseValue <- min(x, na.rm = TRUE) / 2
        if (!is.double(baseValue)) baseValue <- as.double(baseValue)
        if (!is.integer(distance)) distance <- as.integer(distance)
        return(.Call("impute_with_linear_interpolation_base", x, baseValue,
                     distance, PACKAGE = "xcms"))
    }
}

#' @description Calculate breaks for same-sized bins for data values
#' from `fromX` to `toX`.
#'
#' @details This generates bins such as a call to
#' `seq(fromX, toX, length.out = nBins)` would. The first and second element
#' in the result vector thus defines the lower and upper boundary for the first
#' bin, the second and third value for the second bin and so on.
#'
#' @title Generate breaks for binning
#'
#' @param fromX `numeric(1)` specifying the lowest value for the bins.
#'
#' @param toX `numeric(1)` specifying the largest value for the bins.
#'
#' @param nBins `numeric(1)` defining the number of bins.
#'
#' @param shiftByHalfBinSize Logical indicating whether the bins should be
#'     shifted left by half bin size. This results centered bins, i.e. the
#'     first bin being centered at `fromX` and the last around `toX`.
#'
#' @return A numeric vector of length `nBins + 1` defining the lower and
#'     upper bounds of the bins.
#'
#' @author Johannes Rainer
#'
#' @family functions to define bins
#'
#' @seealso [binYonX()] for a binning function.
#'
#' @md
#'
#' @examples
#' ## Create breaks to bin values from 3 to 20 into 20 bins
#' breaks_on_nBins(3, 20, nBins = 20)
#' ## The same call but using shiftByHalfBinSize
#' breaks_on_nBins(3, 20, nBins = 20, shiftByHalfBinSize = TRUE)
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

#' @description
#'
#' Defines breaks for `binSize` sized bins for values ranging
#' from `fromX` to `toX`.
#'
#' @details
#'
#' This function creates breaks for bins of size `binSize`. The
#' function ensures that the full data range is included in the bins, i.e. the
#' last value (upper boundary of the last bin) is always equal `toX`. This
#' however means that the size of the last bin will not always be equal to the
#' desired bin size.
#'
#' See examples for more details and a comparisom to R's `seq()` function.
#'
#' @title Generate breaks for binning using a defined bin size.
#'
#' @param fromX `numeric(1)` specifying the lowest value for the bins.
#'
#' @param toX `numeric(1)` specifying the largest value for the bins.
#'
#' @param binSize `numeric(1)` defining the size of a bin.
#'
#' @return A numeric vector defining the lower and upper bounds of the bins.
#'
#' @author Johannes Rainer
#'
#' @family functions to define bins
#'
#' @seealso [binYonX()] for a binning function.
#'
#' @md
#'
#' @examples
#' ## Define breaks with a size of 0.13 for a data range from 1 to 10:
#' breaks_on_binSize(1, 10, 0.13)
#' ## The size of the last bin is however larger than 0.13:
#' diff(breaks_on_binSize(1, 10, 0.13))
#' ## If we would use seq, the max value would not be included:
#' seq(1, 10, by = 0.13)
#'
#' ## In the next example we use binSize that leads to an additional last bin with
#' ## a smaller binSize:
#' breaks_on_binSize(1, 10, 0.51)
#' ## Again, the max value is included, but the size of the last bin is < 0.51.
#' diff(breaks_on_binSize(1, 10, 0.51))
#' ## Using just seq would result in the following bin definition:
#' seq(1, 10, by = 0.51)
#' ## Thus it defines one bin (break) less.
breaks_on_binSize <- function(fromX, toX, binSize) {
    if(missing(fromX) | missing(toX) | missing(binSize))
        stop("'fromX', 'toX' and 'binSize' are required!")
    if (!is.double(fromX)) fromX <- as.double(fromX)
    if (!is.double(toX)) toX <- as.double(toX)
    if (!is.integer(binSize)) binSize <- as.double(binSize)
    return(.Call("breaks_on_binSize", fromX, toX, binSize, PACKAGE = "xcms"))
}

############################################################
## profBinR
##
## "reference" implementation of profBin in R.
##
## @description This function is a reference implementation in R for the profBin
## C-function.
##
## @details Same as the `profBin` method, this function calculates
## centered breaks (i.e. shifted by half the bin-size, thus the center for the
## first bin is the value of `fromX`) and bins the input vector `x`
## into the bins defined by these breaks. The maximum `y` value corresponding
## to the `x` values per bin is then returned.
## Either `nBins` or `step` has to be defined.
## @title Reference implementation of the profBin C-function in R
## @param x Numeric vector of values used for binning.
## @param y Numeric vector of values to be binned.
## @param nBins Number of bins.
## @param step Bin size.
## @param fromX X value from which binning should start. Default to `min(x)`.
## @param toX Maximum `x` value that should be included in the binning.
## @param FUN Function to aggregate `y` values for `x` values falling
## within a bin.
## @return A numeric vector of the maximum `y` values for each bin.
## @author Johannes Rainer
profBinR <- function(x, y, nBins, step, fromX = min(x), toX = max(x),
                     FUN = max, shiftByHalfBinSize = TRUE) {
    if (missing(x) | missing(y))
        stop("'x' and 'y' are required arguments!")
    if (missing(nBins) & missing(step))
        stop("Either 'nBins' or 'step' has to be specified.")
    ## When we've got nBin
    if (!missing(nBins)) {
        if (shiftByHalfBinSize )
            binSize <- (toX - fromX) / (nBins - 1)
        else
            binSize <- (toX - fromX) / nBins
    }
    ## When we've got step
    if (!missing(step)) {
        binSize <- step
    }
    ## Calculate the bin size
    if (shiftByHalfBinSize)
        brks_xcms <- seq((fromX - binSize/2), (toX + binSize/2),
                         by = binSize)
    else
        brks_xcms <- seq(fromX, toX, by = binSize)
    ## Define output vector:
    res <- numeric(length(brks_xcms) - 1)
    ## Find elements in x falling within the breaks. The breaks span
    ## a little larger range than the values in x, thus we can be sure
    ## that all values are within.
    idx <- findInterval(x, brks_xcms, all.inside = !shiftByHalfBinSize)
    ## Split the y values by index and get the maximum value / bin
    yList <- split(y, idx)
    res[as.numeric(names(yList))] <- unlist(lapply(yList, FUN = FUN))
    return(res)
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
.imputeMethods <- c(0, 1, 2)
names(.imputeMethods) <- c("no", "lin", "linbase")
