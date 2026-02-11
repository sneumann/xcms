# Chromatographic peak summaries

The `chromPeakSummary()` method calculates summary statistics or other
metrics for each of the identified chromatographic peaks in an *xcms*
result object, such as the
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md).
Different metrics can be calculated, depending upon (and configured by)
using dedicated *parameter* classes. As a result, the method returns a
`matrix` or `data.frame` with one row per chromatographic peak. Each
column contains calculated values, depending on the used
method/parameter class.

Currently implemented methods/parameter classes are:

- `BetaDistributionParam`: calculates the *beta_cor* and *beta_snr*
  quality metrics as described in Kumler 2023 representing the result
  from a (correlation) test of similarity (using Pearson's correlation
  coefficient) to a bell curve and the signal-to-noise ratio calculated
  on the residuals of this test.

## Usage

``` r
chromPeakSummary(object, param, ...)

# S4 method for class 'XcmsExperiment,BetaDistributionParam'
chromPeakSummary(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

BetaDistributionParam()
```

## Arguments

- object:

  an *xcms* result object containing information on identified
  chromatographic peaks.

- param:

  a parameter object defining the method/summaries that should be
  calculated (see description above for supported parameter classes).

- ...:

  additional arguments passed to the method implementation.

- msLevel:

  `integer(1)` with the MS level of the chromatographic peaks on which
  the metric should be calculated.

- chunkSize:

  `integer(1)` defining the number of samples from which data should be
  loaded and processed at a time.

- BPPARAM:

  Parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

## Value

A `matrix` or `data.frame` with the same number of rows as there are
chromatographic peaks. Columns contain the calculated values. The number
of columns, their names and content depend on the used parameter object.
See the respective documentation above for more details.

## References

Kumler W, Hazelton B J and Ingalls A E (2023) "Picky with peakpicking:
assessing chromatographic peak quality with simple metrics in
metabolomics" *BMC Bioinformatics* 24(1):404. doi:
[10.1186/s12859-023-05533-4](https://doi.org/10.1186/s12859-023-05533-4)

## Author

Pablo Vangeenderhuysen, Johannes Rainer, William Kumler
