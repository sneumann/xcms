# Correlate chromatograms

**For `xcms` \>= 3.15.3 please use
[`MSnbase::compareChromatograms()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
instead of `correlate`**

Correlate intensities of two chromatograms with each other. If the two
`Chromatogram` objects have different retention times they are first
*aligned* to match data points in the first to data points in the second
chromatogram. See help on `alignRt` in
[`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
for more details.

If `correlate` is called on a single
[`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
object a pairwise correlation of each chromatogram with each other is
performed and a `matrix` with the correlation coefficients is returned.

Note that the correlation of two chromatograms depends also on their
order, e.g. `correlate(chr1, chr2)` might not be identical to
`correlate(chr2, chr1)`. The lower and upper triangular part of the
correlation matrix might thus be different.

## Usage

``` r
# S4 method for class 'Chromatogram,Chromatogram'
correlate(
  x,
  y,
  use = "pairwise.complete.obs",
  method = c("pearson", "kendall", "spearman"),
  align = c("closest", "approx"),
  ...
)

# S4 method for class 'MChromatograms,missing'
correlate(
  x,
  y = NULL,
  use = "pairwise.complete.obs",
  method = c("pearson", "kendall", "spearman"),
  align = c("closest", "approx"),
  ...
)

# S4 method for class 'MChromatograms,MChromatograms'
correlate(
  x,
  y = NULL,
  use = "pairwise.complete.obs",
  method = c("pearson", "kendall", "spearman"),
  align = c("closest", "approx"),
  ...
)
```

## Arguments

- x:

  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
  or
  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  object.

- y:

  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
  or
  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  object.

- use:

  `character(1)` passed to the `cor` function. See
  [`cor()`](https://rdrr.io/r/stats/cor.html) for details.

- method:

  `character(1)` passed to the `cor` function. See
  [`stats::cor()`](https://rdrr.io/r/stats/cor.html) for details.

- align:

  `character(1)` defining the alignment method to be used. See help on
  `alignRt` in
  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
  for details. The value of this parameter is passed to the `method`
  parameter of `alignRt`.

- ...:

  optional parameters passed along to the `alignRt` method such as
  `tolerance` that, if set to `0` requires the retention times to be
  identical.

## Value

`numeric(1)` or `matrix` (if called on `MChromatograms` objects) with
the correlation coefficient. If a `matrix` is returned, the rows
represent the chromatograms in `x` and the columns the chromatograms in
`y`.

## Author

Michael Witting, Johannes Rainer

## Examples

``` r

library(MSnbase)
chr1 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))
chr2 <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(80, 50, 20, 10, 9, 4, 3, 4, 1, 3))
chr3 <- Chromatogram(rtime = 3:9 + rnorm(7, sd = 0.3),
    intensity = c(53, 80, 130, 15, 5, 3, 2))

chrs <- MChromatograms(list(chr1, chr2, chr3))

## Using `compareChromatograms` instead of `correlate`.
compareChromatograms(chr1, chr2)
#> [1] -0.02314823
compareChromatograms(chr2, chr1)
#> [1] 0.2630029

compareChromatograms(chrs, chrs)
#>           [,1]        [,2]      [,3]
#> [1,] 1.0000000 -0.02314823 0.9954772
#> [2,] 0.2630029  1.00000000 0.5341594
#> [3,] 0.9954772  0.53415936 1.0000000
```
