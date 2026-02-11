# DEPRECATED: Extract a `data.frame` containing MS data

**UPDATE**: the `extractMsData` and `plotMsData` functions are
deprecated and `as(x, "data.frame")` and `plot(x, type = "XIC")` (`x`
being an `OnDiskMSnExp` or `XCMSnExp` object) should be used instead.
See examples below. Be aware that filtering the raw object might however
drop the adjusted retention times. In such cases it is advisable to use
the
[`applyAdjustedRtime()`](https://sneumann.github.io/xcms/reference/applyAdjustedRtime.md)
function prior to filtering.

Extract a `data.frame` of retention time, mz and intensity values from
each file/sample in the provided rt-mz range (or for the full data range
if `rt` and `mz` are not defined).

## Usage

``` r
# S4 method for class 'OnDiskMSnExp'
extractMsData(object, rt, mz, msLevel = 1L)

# S4 method for class 'XCMSnExp'
extractMsData(
  object,
  rt,
  mz,
  msLevel = 1L,
  adjustedRtime = hasAdjustedRtime(object)
)
```

## Arguments

- object:

  A `XCMSnExp` or `OnDiskMSnExp` object.

- rt:

  `numeric(2)` with the retention time range from which the data should
  be extracted.

- mz:

  `numeric(2)` with the mz range.

- msLevel:

  `integer` defining the MS level(s) to which the data should be
  sub-setted prior to extraction; defaults to `msLevel = 1L`.

- adjustedRtime:

  (for `extractMsData,XCMSnExp`): `logical(1)` specifying if adjusted or
  raw retention times should be reported. Defaults to adjusted retention
  times, if these are present in `object`.

## Value

A `list` of length equal to the number of samples/files in `object`.
Each element being a `data.frame` with columns `"rt"`, `"mz"` and `"i"`
with the retention time, mz and intensity tuples of a file. If no data
is available for the mz-rt range in a file a `data.frame` with 0 rows is
returned for that file.

## See also

`XCMSnExp` for the data object.

## Author

Johannes Rainer

## Examples

``` r

## Load a test data set with detected peaks
library(MSnbase)
data(faahko_sub)
## Update the path to the files for the local system
dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")

## Disable parallel processing for this example
register(SerialParam())

## Extract the full MS data for a certain retention time range
## as a data.frame
tmp <- filterRt(faahko_sub, rt = c(2800, 2900))
ms_all <- as(tmp, "data.frame")
head(ms_all)
#>   file       rt    mz     i
#> 1    1 2800.284 200.1   968
#> 2    1 2800.284 201.0  1273
#> 3    1 2800.284 201.9  1729
#> 4    1 2800.284 203.0   859
#> 5    1 2800.284 205.0 24136
#> 6    1 2800.284 206.0  3913
nrow(ms_all)
#> [1] 79960
```
