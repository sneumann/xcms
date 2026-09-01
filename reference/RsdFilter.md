# Filter features based on their coefficient of variation

The `RsdFilter` class and methods enable users to filter features from
an `XcmsExperiment` or `SummarizedExperiment` object based on their
relative standard deviation (coefficient of variation) for a specified
threshold.

This `filter` is part of the possible dispatch of the generic function
`filterFeatures`. Features *above* (`>`) the user-input threshold will
be removed from the entire dataset.

## Usage

``` r
RsdFilter(threshold = 0.3, qcIndex = integer(), na.rm = TRUE, mad = FALSE)

# S4 method for class 'XcmsResult,RsdFilter'
filterFeatures(object, filter, ...)

# S4 method for class 'SummarizedExperiment,RsdFilter'
filterFeatures(object, filter, assay = 1)
```

## Arguments

- threshold:

  `numeric` value representing the threshold. Features with a
  coefficient of variation *strictly higher* (`>`) than this will be
  removed from the entire dataset.

- qcIndex:

  `integer` (or `logical`) vector corresponding to the indices of QC
  samples.

- na.rm:

  `logical` indicates whether missing values (`NA`) should be removed
  prior to the calculations.

- mad:

  `logical` indicates whether the *Median Absolute Deviation* (MAD)
  should be used instead of the standard deviation. This is suggested
  for non-gaussian distributed data.

- object:

  `XcmsExperiment` or `SummarizedExperiment`. For an `XcmsExperiment`
  object, the `featureValues(object)` will be evaluated, and for
  `Summarizedesxperiment` the `assay(object, assay)`. The object will be
  filtered.

- filter:

  The parameter object selecting and configuring the type of filtering.
  It can be one of the following classes: `RsdFilter`,
  [`DratioFilter`](https://sneumann.github.io/xcms/reference/DratioFilter.md),
  [`PercentMissingFilter`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md)
  or
  [`BlankFlag`](https://sneumann.github.io/xcms/reference/BlankFlag.md).

- ...:

  Optional parameters. For `object` being an `XcmsExperiment`:
  parameters for the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  call.

- assay:

  For filtering of `SummarizedExperiment` objects only. Indicates which
  assay the filtering will be based on. Note that the features for the
  entire object will be removed, but the computations are performed on a
  single assay. Default is 1, which means the first assay of the
  `object` will be evaluated.

## Value

For `RsdFilter`: a `RsdFilter` class. `filterFeatures` return the input
object minus the features that did not met the user input threshold.

## Note

It is assumed that the abundance values are in natural scale. Abundances
in log scale should be first transformed to natural scale before
calculating the RSD.

## See also

Other Filter features in xcms:
[`BlankFlag()`](https://sneumann.github.io/xcms/reference/BlankFlag.md),
[`DratioFilter()`](https://sneumann.github.io/xcms/reference/DratioFilter.md),
[`PercentMissingFilter()`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md)

## Author

Philippine Louail
