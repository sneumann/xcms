# Generate breaks for binning using a defined bin size.

Defines breaks for `binSize` sized bins for values ranging from `fromX`
to `toX`.

## Usage

``` r
breaks_on_binSize(fromX, toX, binSize)
```

## Arguments

- fromX:

  `numeric(1)` specifying the lowest value for the bins.

- toX:

  `numeric(1)` specifying the largest value for the bins.

- binSize:

  `numeric(1)` defining the size of a bin.

## Value

A numeric vector defining the lower and upper bounds of the bins.

## Details

This function creates breaks for bins of size `binSize`. The function
ensures that the full data range is included in the bins, i.e. the last
value (upper boundary of the last bin) is always equal `toX`. This
however means that the size of the last bin will not always be equal to
the desired bin size.

See examples for more details and a comparisom to R's
[`seq()`](https://rdrr.io/r/base/seq.html) function.

## See also

[`binYonX()`](https://sneumann.github.io/xcms/reference/binYonX.md) for
a binning function.

Other functions to define bins:
[`breaks_on_nBins()`](https://sneumann.github.io/xcms/reference/breaks_on_nBins.md)

## Author

Johannes Rainer

## Examples

``` r
## Define breaks with a size of 0.13 for a data range from 1 to 10:
breaks_on_binSize(1, 10, 0.13)
#>  [1]  1.00  1.13  1.26  1.39  1.52  1.65  1.78  1.91  2.04  2.17  2.30  2.43
#> [13]  2.56  2.69  2.82  2.95  3.08  3.21  3.34  3.47  3.60  3.73  3.86  3.99
#> [25]  4.12  4.25  4.38  4.51  4.64  4.77  4.90  5.03  5.16  5.29  5.42  5.55
#> [37]  5.68  5.81  5.94  6.07  6.20  6.33  6.46  6.59  6.72  6.85  6.98  7.11
#> [49]  7.24  7.37  7.50  7.63  7.76  7.89  8.02  8.15  8.28  8.41  8.54  8.67
#> [61]  8.80  8.93  9.06  9.19  9.32  9.45  9.58  9.71  9.84 10.00
## The size of the last bin is however larger than 0.13:
diff(breaks_on_binSize(1, 10, 0.13))
#>  [1] 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13
#> [16] 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13
#> [31] 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13
#> [46] 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13
#> [61] 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.13 0.16
## If we would use seq, the max value would not be included:
seq(1, 10, by = 0.13)
#>  [1] 1.00 1.13 1.26 1.39 1.52 1.65 1.78 1.91 2.04 2.17 2.30 2.43 2.56 2.69 2.82
#> [16] 2.95 3.08 3.21 3.34 3.47 3.60 3.73 3.86 3.99 4.12 4.25 4.38 4.51 4.64 4.77
#> [31] 4.90 5.03 5.16 5.29 5.42 5.55 5.68 5.81 5.94 6.07 6.20 6.33 6.46 6.59 6.72
#> [46] 6.85 6.98 7.11 7.24 7.37 7.50 7.63 7.76 7.89 8.02 8.15 8.28 8.41 8.54 8.67
#> [61] 8.80 8.93 9.06 9.19 9.32 9.45 9.58 9.71 9.84 9.97

## In the next example we use binSize that leads to an additional last bin with
## a smaller binSize:
breaks_on_binSize(1, 10, 0.51)
#>  [1]  1.00  1.51  2.02  2.53  3.04  3.55  4.06  4.57  5.08  5.59  6.10  6.61
#> [13]  7.12  7.63  8.14  8.65  9.16  9.67 10.00
## Again, the max value is included, but the size of the last bin is < 0.51.
diff(breaks_on_binSize(1, 10, 0.51))
#>  [1] 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51 0.51
#> [16] 0.51 0.51 0.33
## Using just seq would result in the following bin definition:
seq(1, 10, by = 0.51)
#>  [1] 1.00 1.51 2.02 2.53 3.04 3.55 4.06 4.57 5.08 5.59 6.10 6.61 7.12 7.63 8.14
#> [16] 8.65 9.16 9.67
## Thus it defines one bin (break) less.
```
