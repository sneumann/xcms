# Flag features based on the intensity in blank samples

The `BlankFlag` class and method enable users to flag features of an
`XcmsExperiment` or `SummarizedExperiment` object based on the
relationship between the intensity of a feature in blanks compared to
the intensity in the samples.

This class and method are part of the possible dispatch of the generic
function `filterFeatures`. Features *below* (`<`) the user-input
threshold will be flagged by calling the `filterFeatures` function. This
means that an extra column will be created in `featureDefinitions` or
`rowData` called `possible_contaminants` with a logical value for each
feature.

## Usage

``` r
BlankFlag(
  threshold = 2,
  blankIndex = integer(),
  qcIndex = integer(),
  na.rm = TRUE
)

# S4 method for class 'XcmsResult,BlankFlag'
filterFeatures(object, filter, ...)

# S4 method for class 'SummarizedExperiment,BlankFlag'
filterFeatures(object, filter, assay = 1)
```

## Arguments

- threshold:

  `numeric` indicates the minimum difference required between the mean
  abundance of a feature in samples compared to the mean abundance of
  the same feature in blanks for it to not be considered a possible
  contaminant. For example, the default threshold of 2 signifies that
  the mean abundance of the features in samples has to be at least twice
  the mean abundance in blanks for it to not be flagged as a possible
  contaminant.

- blankIndex:

  `integer` (or `logical`) vector corresponding to the indices of blank
  samples.

- qcIndex:

  `integer` (or `logical`) vector corresponding to the indices of
  quality control (QC) samples.

- na.rm:

  `logical` indicates whether missing values (`NA`) should be removed
  prior to the calculations.

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
  [`PercentMissingFilter`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md)
  or `BlankFlag`.

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

For `BlankFlag`: a `BlankFlag` class. `filterFeatures` returns the input
object with an added column in the features metadata called
`possible_contaminants` with a logical value for each feature. This is
added to `featureDefinitions` for `XcmsExperiment` objects and `rowData`
for `SummarizedExperiment` objects.

## See also

Other Filter features in xcms:
[`DratioFilter()`](https://sneumann.github.io/xcms/reference/DratioFilter.md),
[`PercentMissingFilter()`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md),
[`RsdFilter()`](https://sneumann.github.io/xcms/reference/RsdFilter.md)

## Author

Philippine Louail
