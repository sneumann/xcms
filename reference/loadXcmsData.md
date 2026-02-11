# LC-MS preprocessing result test data sets

Data sets with `xcms` preprocessing results are provided within the
`xcms` package and can be loaded with the `loadXcmsData` function. The
available Test data sets are:

- `xdata`: an
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with the results from a `xcms`-based pre-processing of an LC-MS
  untargeted metabolomics data set. The raw data files are provided in
  the `faahKO` R package.

- `xmse`: an
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  object with the results from an `xcms`-based pre-processing of an
  LC-MS untargeted metabolomics data set (same original data set and
  pre-processing settings as for the `xdata` data set). The
  pre-processing of this data set is described in detail in the *xcms*
  vignette of the `xcms` package.

- `faahko_sub`: an
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with identified chromatographic peaks in 3 samples from the
  data files in the `faahKO` R package.

- `faahko_sub2`: an
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  object with identified chromatographic peaks in 3 samples from the
  data files in the `faahKO` R package.

Data sets can also be loaded using `data`, which would however require
to update objects to point to the location of the raw data files. The
`loadXcmsData` loads the data and ensures that all paths are updated
accordingly.

## Usage

``` r
loadXcmsData(x = c("xmse", "xdata", "faahko_sub", "faahko_sub2"))
```

## Arguments

- x:

  For `loadXcmsData`: `character(1)` with the name of the data file
  (object) to load.

## Examples

``` r

library(xcms)
xdata <- loadXcmsData()
```
