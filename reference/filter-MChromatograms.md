# Filtering sets of chromatographic data

These functions allow to filter (subset)
[`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
or
[`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
objects, i.e. sets of chromatographic data, without changing the data
(intensity and retention times) within the individual chromatograms
([`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
objects).

- `filterColumnsIntensityAbove`: subsets a `MChromatograms` objects
  keeping only columns (samples) for which `value` is larger than the
  provided `threshold` in `which` rows (i.e. if `which = "any"` a column
  is kept if **any** of the chromatograms in that column have a `value`
  larger than `threshold` or with `which = "all"` **all** chromatograms
  in that column fulfill this criteria). Parameter `value` allows to
  define on which value the comparison should be performed, with
  `value = "bpi"` the maximum intensity of each chromatogram is compared
  to `threshold`, with
  `value = "tic" the total sum of intensities of each chromatogram is compared to `threshold`. For `XChromatograms`object,`value
  = "maxo"`and`value =
  "into"`are supported which compares the largest intensity of all identified chromatographic peaks in the chromatogram with`threshold\`,
  or the integrated peak area, respectively.

- `filterColumnsKeepTop`: subsets a `MChromatograms` object keeping the
  top `n` columns sorted by the value specified with `sortBy`. In
  detail, for each column the value defined by `sortBy` is extracted
  from each chromatogram and aggregated using the `aggregationFun`.
  Thus, by default, for each chromatogram the maximum intensity is
  determined (`sortBy = "bpi"`) and these values are summed up for
  chromatograms in the same column (`aggregationFun = sum`). The columns
  are then sorted by these values and the top `n` columns are retained
  in the returned `MChromatograms`. Similar to the
  `filterColumnsIntensityAbove` function, this function allows to use
  for `XChromatograms` objects to sort the columns by column
  `sortBy = "maxo"` or `sortBy = "into"` of the `chromPeaks` matrix.

## Usage

``` r
# S4 method for class 'MChromatograms'
filterColumnsIntensityAbove(
  object,
  threshold = 0,
  value = c("bpi", "tic"),
  which = c("any", "all")
)

# S4 method for class 'MChromatograms'
filterColumnsKeepTop(
  object,
  n = 1L,
  sortBy = c("bpi", "tic"),
  aggregationFun = sum
)

# S4 method for class 'XChromatograms'
filterColumnsIntensityAbove(
  object,
  threshold = 0,
  value = c("bpi", "tic", "maxo", "into"),
  which = c("any", "all")
)

# S4 method for class 'XChromatograms'
filterColumnsKeepTop(
  object,
  n = 1L,
  sortBy = c("bpi", "tic", "maxo", "into"),
  aggregationFun = sum
)
```

## Arguments

- object:

  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  or
  [`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  object.

- threshold:

  for `filterColumnsIntensityAbove`: `numeric(1)` with the threshold
  value to compare against.

- value:

  `character(1)` defining which value should be used in the comparison
  or sorting. Can be `value = "bpi"` (default) to use the maximum
  intensity per chromatogram or `value = "tic"` to use the sum of
  intensities per chromatogram. For
  [`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  objects also `value = "maxo"` and `value = "into"` is supported to use
  the maximum intensity or the integrated area of identified
  chromatographic peaks in each chromatogram.

- which:

  for `filterColumnsIntensityAbove`: `character(1)` defining whether
  **any** (`which = "any"`, default) or **all** (`which = "all"`)
  chromatograms in a column have to fulfill the criteria for the column
  to be kept.

- n:

  for `filterColumnsKeepTop`: `integer(1)` specifying the number of
  columns that should be returned. `n` will be rounded to the closest
  (larger) integer value.

- sortBy:

  for `filterColumnsKeepTop`: the value by which columns should be
  ordered to determine the top n columns. Can be either `sortBy = "bpi"`
  (the default), in which case the maximum intensity of each column's
  chromatograms is used, or `sortBy = "tic"` to use the total intensity
  sum of all chromatograms. For
  [`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  objects also `value = "maxo"` and `value = "into"` is supported to use
  the maximum intensity or the integrated area of identified
  chromatographic peaks in each chromatogram.

- aggregationFun:

  for `filterColumnsKeepTop`: function to be used to aggregate (combine)
  the values from all chromatograms in each column. Defaults to
  `aggregationFun = sum` in which case the sum of the values is used to
  rank the columns. Alternatively the `mean`, `median` or similar
  function can be used.

## Value

a filtered `MChromatograms` (or `XChromatograms`) object with the same
number of rows (EICs) but eventually a lower number of columns
(samples).

## Author

Johannes Rainer

## Examples

``` r

library(MSnbase)
chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
    intensity = c(53, 80, 130, 15, 5, 3, 2))

chrs <- MChromatograms(list(chr1, chr2, chr1, chr3, chr2, chr3),
    ncol = 3, byrow = FALSE)
chrs
#> MChromatograms with 2 rows and 3 columns
#>                   1              2              3
#>      <Chromatogram> <Chromatogram> <Chromatogram>
#> [1,]     length: 10     length: 10     length: 10
#> [2,]     length: 10      length: 7      length: 7
#> phenoData with 0 variables
#> featureData with 0 variables

#### filterColumnsIntensityAbove
##
## Keep all columns with for which the maximum intensity of any of its
## chromatograms is larger 90
filterColumnsIntensityAbove(chrs, threshold = 90)
#> MChromatograms with 2 rows and 3 columns
#>                   1              2              3
#>      <Chromatogram> <Chromatogram> <Chromatogram>
#> [1,]     length: 10     length: 10     length: 10
#> [2,]     length: 10      length: 7      length: 7
#> phenoData with 0 variables
#> featureData with 0 variables

## Require that ALL chromatograms in a column have a value larger 90
filterColumnsIntensityAbove(chrs, threshold = 90, which = "all")
#> MChromatograms with 2 rows and 1 column
#>                   2
#>      <Chromatogram>
#> [1,]     length: 10
#> [2,]      length: 7
#> phenoData with 0 variables
#> featureData with 0 variables

## If none of the columns fulfills the criteria no columns are returned
filterColumnsIntensityAbove(chrs, threshold = 900)
#> MChromatograms with 2 rows and 0 columns
#> phenoData with 0 variables
#> featureData with 0 variables

## Filtering XChromatograms allow in addition to filter on the columns
## "maxo" or "into" of the identified chromatographic peaks within each
## chromatogram.

#### filterColumnsKeepTop
##
## Keep the 2 columns with the highest sum of maximal intensities in their
## chromatograms
filterColumnsKeepTop(chrs, n = 1)
#> MChromatograms with 2 rows and 1 column
#>                   2
#>      <Chromatogram>
#> [1,]     length: 10
#> [2,]      length: 7
#> phenoData with 0 variables
#> featureData with 0 variables

## Keep the 50 percent of columns with the highest total sum of signal. Note
## that n will be rounded to the next larger integer value
filterColumnsKeepTop(chrs, n = 0.5 * ncol(chrs), sortBy = "tic")
#> MChromatograms with 2 rows and 2 columns
#>                   2              3
#>      <Chromatogram> <Chromatogram>
#> [1,]     length: 10     length: 10
#> [2,]      length: 7      length: 7
#> phenoData with 0 variables
#> featureData with 0 variables
```
