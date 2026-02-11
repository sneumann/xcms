# Aggregate values in y for bins defined on x

This functions takes two same-sized numeric vectors `x` and `y`,
bins/cuts `x` into bins (either a pre-defined number of equal-sized bins
or bins of a pre-defined size) and aggregates values in `y`
corresponding to `x` values falling within each bin. By default (i.e.
`method = "max"`) the maximal `y` value for the corresponding `x` values
is identified. `x` is expected to be incrementally sorted and, if not,
it will be internally sorted (in which case also `y` will be ordered
according to the order of `x`).

## Usage

``` r
binYonX(
  x,
  y,
  breaks,
  nBins,
  binSize,
  binFromX,
  binToX,
  fromIdx = 1L,
  toIdx = length(x),
  method = "max",
  baseValue,
  sortedX = !is.unsorted(x),
  shiftByHalfBinSize = FALSE,
  returnIndex = FALSE,
  returnX = TRUE
)
```

## Arguments

- x:

  Numeric vector to be used for binning.

- y:

  Numeric vector (same length than `x`) from which the maximum values
  for each bin should be defined. If not provided, `x` will be used.

- breaks:

  Numeric vector defining the breaks for the bins, i.e. the lower and
  upper values for each bin. See examples below.

- nBins:

  `integer(1)` defining the number of desired bins.

- binSize:

  `numeric(1)` defining the desired bin size.

- binFromX:

  Optional `numeric(1)` allowing to manually specify the range of
  x-values to be used for binning. This will affect only the calculation
  of the breaks for the bins (i.e. if `nBins` or `binSize` is provided).
  If not provided the minimal value in the sub-set `fromIdx`-`toIdx` in
  input vector `x` will be used.

- binToX:

  Same as `binFromX`, but defining the maximum x-value to be used for
  binning.

- fromIdx:

  Integer vector defining the start position of one or multiple sub-sets
  of input vector `x` that should be used for binning.

- toIdx:

  Same as `toIdx`, but defining the maximum index (or indices) in x to
  be used for binning.

- method:

  A character string specifying the method that should be used to
  aggregate values in `y`. Allowed are `"max"`, `"min"`, `"sum"` and
  `"mean"` to identify the maximal or minimal value or to sum all values
  within a bin or calculate their mean value.

- baseValue:

  The base value for empty bins (i.e. bins into which either no values
  in `x` did fall, or to which only `NA` values in `y` were assigned).
  By default (i.e. if not specified), `NA` is assigned to such bins.

- sortedX:

  Whether `x` is sorted.

- shiftByHalfBinSize:

  Logical specifying whether the bins should be shifted by half the bin
  size to the left. Thus, the first bin will have its center at `fromX`
  and its lower and upper boundary are `fromX - binSize/2` and
  `fromX + binSize/2`. This argument is ignored if `breaks` are
  provided.

- returnIndex:

  Logical indicating whether the index of the max (if `method = "max"`)
  or min (if `method = "min"`) value within each bin in input vector `x`
  should also be reported. For methods other than `"max"` or `"min"`
  this argument is ignored.

- returnX:

  `logical` allowing to avoid returning `$x`, i.e. the mid-points of the
  bins. `returnX = FALSE` might be useful in cases where `breaks` are
  pre-defined as it considerably reduces the memory demand.

## Value

Returns a list of length 2, the first element (named `"x"`) contains the
bin mid-points, the second element (named `"y"`) the aggregated values
from input vector `y` within each bin. For `returnIndex = TRUE` the list
contains an additional element `"index"` with the index of the max or
min (depending on whether `method = "max"` or `method = "min"`) value
within each bin in input vector `x`.

## Details

The breaks defining the boundary of each bin can be either passed
directly to the function with the argument `breaks`, or are calculated
on the data based on arguments `nBins` or `binSize` along with
`fromIdx`, `toIdx` and optionally `binFromX` and `binToX`. Arguments
`fromIdx` and `toIdx` allow to specify subset(s) of the input vector `x`
on which bins should be calculated. The default the full `x` vector is
considered. Also, if not specified otherwise with arguments `binFromX`
and `binToX`, the range of the bins within each of the sub-sets will be
from `x[fromIdx]` to `x[toIdx]`. Arguments `binFromX` and `binToX` allow
to overwrite this by manually defining the a range on which the breaks
should be calculated. See examples below for more details.

    Calculation of breaks: for `nBins` the breaks correspond to
    `seq(min(x[fromIdx])), max(x[fromIdx], length.out = (nBins + 1))`.
    For `binSize` the breaks correspond to
    `seq(min(x[fromIdx]), max(x[toIdx]), by = binSize)` with the
    exception that the last break value is forced to be equal to
    `max(x[toIdx])`. This ensures that all values from the specified
    range are covered by the breaks defining the bins. The last bin could
    however in some instances be slightly larger than `binSize`. See
    [breaks_on_binSize()] and [breaks_on_nBins()] for
    more details.

## Note

The function ensures that all values within the range used to define the
breaks are considered in the binning (and assigned to a bin). This means
that for all bins except the last one values in `x` have to be
`>= xlower` and `< xupper` (with `xlower` and `xupper` being the lower
and upper boundary, respectively). For the last bin the condition is
`x >= xlower & x <= xupper`. Note also that if `shiftByHalfBinSize` is
`TRUE` the range of values that is used for binning is expanded by
`binSize` (i.e. the lower boundary will be `fromX - binSize/2`, the
upper `toX + binSize/2`). Setting this argument to `TRUE` resembles the
binning that is/was used in `profBin` function from *xcms* \< 1.51.

    `NA` handling: by default the function ignores `NA` values in
    `y` (thus inherently assumes `na.rm = TRUE`). No `NA`
    values are allowed in `x`.

## See also

[`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)

## Author

Johannes Rainer

## Examples

``` r
########
## Simple example illustrating the breaks and the binning.
##
## Define breaks for 5 bins:
brks <- seq(2, 12, length.out = 6)
## The first bin is then [2,4), the second [4,6) and so on.
brks
#> [1]  2  4  6  8 10 12
## Get the max value falling within each bin.
binYonX(x = 1:16, y = 1:16, breaks = brks)
#> $x
#> [1]  3  5  7  9 11
#> 
#> $y
#> [1]  3  5  7  9 12
#> 
## Thus, the largest value in x = 1:16 falling into the bin [2,4) (i.e. being
## >= 2 and < 4) is 3, the largest one falling into [4,6) is 5 and so on.
## Note however the function ensures that the minimal and maximal x-value
## (in this example 1 and 12) fall within a bin, i.e. 12 is considered for
## the last bin.

#######
## Performing the binning ons sub-set of x
##
X <- 1:16
## Bin X from element 4 to 10 into 5 bins.
X[4:10]
#> [1]  4  5  6  7  8  9 10
binYonX(X, X, nBins = 5L, fromIdx = 4, toIdx = 10)
#> $x
#> [1] 4.6 5.8 7.0 8.2 9.4
#> 
#> $y
#> [1]  5  6  7  8 10
#> 
## This defines breaks for 5 bins on the values from 4 to 10 and bins
## the values into these 5 bins. Alternatively, we could manually specify
## the range for the binning, i.e. the minimal and maximal value for the
## breaks:
binYonX(X, X, nBins = 5L, fromIdx = 4, toIdx = 10, binFromX = 1, binToX = 16)
#> $x
#> [1]  2.5  5.5  8.5 11.5 14.5
#> 
#> $y
#> [1] NA  6  9 10 NA
#> 
## In this case the breaks for 5 bins were defined from a value 1 to 16 and
## the values 4 to 10 were binned based on these breaks.

#######
## Bin values within a sub-set of x, second example
##
## This example illustrates how the fromIdx and toIdx parameters can be used.
## x defines 3 times the sequence form 1 to 10, while y is the sequence from
## 1 to 30. In this very simple example x is supposed to represent M/Z values
## from 3 consecutive scans and y the intensities measured for each M/Z in
## each scan. We want to get the maximum intensities for M/Z value bins only
## for the second scan, and thus we use fromIdx = 11 and toIdx = 20. The breaks
## for the bins are defined with the nBins, binFromX and binToX.
X <- rep(1:10, 3)
Y <- 1:30
## Bin the M/Z values in the second scan into 5 bins and get the maximum
## intensity for each bin. Note that we have to specify sortedX = TRUE as
## the x and y vectors would be sorted otherwise.
binYonX(X, Y, nBins = 5L, sortedX = TRUE, fromIdx = 11, toIdx = 20)
#> $x
#> [1] 1.9 3.7 5.5 7.3 9.1
#> 
#> $y
#> [1] 12 14 16 18 20
#> 

#######
## Bin in overlapping sub-sets of X
##
## In this example we define overlapping sub-sets of X and perform the binning
## within these.
X <- 1:30
## Define the start and end indices of the sub-sets.
fIdx <- c(2, 8, 21)
tIdx <- c(10, 25, 30)
binYonX(X, nBins = 5L, fromIdx = fIdx, toIdx = tIdx)
#> [[1]]
#> [[1]]$x
#> [1] 2.8 4.4 6.0 7.6 9.2
#> 
#> [[1]]$y
#> [1]  3  5  6  8 10
#> 
#> 
#> [[2]]
#> [[2]]$x
#> [1]  9.7 13.1 16.5 19.9 23.3
#> 
#> [[2]]$y
#> [1] 11 14 18 21 25
#> 
#> 
#> [[3]]
#> [[3]]$x
#> [1] 21.9 23.7 25.5 27.3 29.1
#> 
#> [[3]]$y
#> [1] 22 24 26 28 30
#> 
#> 
## The same, but pre-defining also the desired range of the bins.
binYonX(X, nBins = 5L, fromIdx = fIdx, toIdx = tIdx, binFromX = 4, binToX = 28)
#> [[1]]
#> [[1]]$x
#> [1]  6.4 11.2 16.0 20.8 25.6
#> 
#> [[1]]$y
#> [1]  8 10 NA NA NA
#> 
#> 
#> [[2]]
#> [[2]]$x
#> [1]  6.4 11.2 16.0 20.8 25.6
#> 
#> [[2]]$y
#> [1]  8 13 18 23 25
#> 
#> 
#> [[3]]
#> [[3]]$x
#> [1]  6.4 11.2 16.0 20.8 25.6
#> 
#> [[3]]$y
#> [1] NA NA NA 23 28
#> 
#> 
## The same bins are thus used for each sub-set.
```
