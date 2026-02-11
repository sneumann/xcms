# Subset an xcmsRaw object by scans

Subset an `xcmsRaw` object by scans. The returned `xcmsRaw` object
contains values for all scans specified with argument `i`. Note that the
`scanrange` slot of the returned `xcmsRaw` will be
`c(1, length(object@scantime))` and hence not `range(i)`.

## Usage

``` r
# S4 method for class 'xcmsRaw,logicalOrNumeric,missing,missing'
x[i, j, drop]
```

## Arguments

- x:

  The `xcmsRaw` object that should be sub-setted.

- i:

  Integer or logical vector specifying the scans/spectra to which `x`
  should be sub-setted.

- j:

  Not supported.

- drop:

  Not supported.

## Value

The sub-setted `xcmsRaw` object.

## Details

Only subsetting by scan index in increasing order or by a logical vector
are supported. If not ordered, argument `i` is sorted automatically.
Indices which are larger than the total number of scans are discarded.

## Author

Johannes Rainer

## Examples

``` r
## Load a test file
file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
xraw <- xcmsRaw(file, profstep = 0)
## The number of scans/spectra:
length(xraw@scantime)
#> [1] 1278

## Subset the object to scans with a scan time from 3500 to 4000.
xsub <- xraw[xraw@scantime >= 3500 & xraw@scantime <= 4000]
range(xsub@scantime)
#> [1] 3501.383 3999.039
## The number of scans:
length(xsub@scantime)
#> [1] 319
## The number of values of the subset:
length(xsub@env$mz)
#> [1] 130569
```
