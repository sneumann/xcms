# Generate breaks for binning

Calculate breaks for same-sized bins for data values from `fromX` to
`toX`.

## Usage

``` r
breaks_on_nBins(fromX, toX, nBins, shiftByHalfBinSize = FALSE)
```

## Arguments

- fromX:

  `numeric(1)` specifying the lowest value for the bins.

- toX:

  `numeric(1)` specifying the largest value for the bins.

- nBins:

  `numeric(1)` defining the number of bins.

- shiftByHalfBinSize:

  Logical indicating whether the bins should be shifted left by half bin
  size. This results centered bins, i.e. the first bin being centered at
  `fromX` and the last around `toX`.

## Value

A numeric vector of length `nBins + 1` defining the lower and upper
bounds of the bins.

## Details

This generates bins such as a call to
`seq(fromX, toX, length.out = nBins)` would. The first and second
element in the result vector thus defines the lower and upper boundary
for the first bin, the second and third value for the second bin and so
on.

## See also

[`binYonX()`](https://sneumann.github.io/xcms/reference/binYonX.md) for
a binning function.

Other functions to define bins:
[`breaks_on_binSize()`](https://sneumann.github.io/xcms/reference/breaks_on_binSize.md)

## Author

Johannes Rainer

## Examples

``` r
## Create breaks to bin values from 3 to 20 into 20 bins
breaks_on_nBins(3, 20, nBins = 20)
#>  [1]  3.00  3.85  4.70  5.55  6.40  7.25  8.10  8.95  9.80 10.65 11.50 12.35
#> [13] 13.20 14.05 14.90 15.75 16.60 17.45 18.30 19.15 20.00
## The same call but using shiftByHalfBinSize
breaks_on_nBins(3, 20, nBins = 20, shiftByHalfBinSize = TRUE)
#>  [1]  2.552632  3.447368  4.342105  5.236842  6.131579  7.026316  7.921053
#>  [8]  8.815789  9.710526 10.605263 11.500000 12.394737 13.289474 14.184211
#> [15] 15.078947 15.973684 16.868421 17.763158 18.657895 19.552632 20.447368
```
