# Filter features based on the percentage of missing data

The `PercentMissingFilter` class and method enable users to filter
features from an `XcmsExperiment` or `SummarizedExperiment` object based
on the percentage (values from 1 to 100) of missing values for each
features in different sample groups and filters them according to a
provided threshold.

This `filter` is part of the possible dispatch of the generic function
`filterFeatures`. Features with a percentage of missing values *higher*
(`>`) than the user input threshold in all sample groups will be removed
(i.e. features for which the proportion of missing values is below
(`<=`) the threshold in at least one sample group will be retained).

## Usage

``` r
PercentMissingFilter(threshold = 30, f = factor())

# S4 method for class 'XcmsResult,PercentMissingFilter'
filterFeatures(object, filter, ...)

# S4 method for class 'SummarizedExperiment,PercentMissingFilter'
filterFeatures(object, filter, assay = 1)
```

## Arguments

- threshold:

  `numeric` percentage (between 0 and 100) of accepted missing values
  for a feature in one sample group.

- f:

  `vector` of the same length as the `object`, specifying the sample
  type for each sample in the dataset. The percentage of missing values
  per feature will be computed within each of these sample groups.
  Parameter `f`, if not already a `factor`, will be converted to one
  using the factor function. Samples with an `NA` as their value in `f`
  will be excluded from calculation.

- object:

  `XcmsExperiment` or `SummarizedExperiment`. For an `XcmsExperiment`
  object, the `featureValues(object)` will be evaluated, and for
  `Summarizedesxperiment` the `assay(object, assay)`. The object will be
  filtered.

- filter:

  The parameter object selecting and configuring the type of filtering.
  It can be one of the following classes:
  [`RsdFilter`](https://sneumann.github.io/xcms/reference/RsdFilter.md),
  [`DratioFilter`](https://sneumann.github.io/xcms/reference/DratioFilter.md),
  `PercentMissingFilter` or
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

For `PercentMissingFilter`: a `PercentMissingFilter` class.
`filterFeatures` return the input object minus the features that did not
met the user input threshold

## See also

Other Filter features in xcms:
[`BlankFlag`](https://sneumann.github.io/xcms/reference/BlankFlag.md),
[`DratioFilter`](https://sneumann.github.io/xcms/reference/DratioFilter.md),
[`RsdFilter`](https://sneumann.github.io/xcms/reference/RsdFilter.md)

## Author

Philippine Louail
