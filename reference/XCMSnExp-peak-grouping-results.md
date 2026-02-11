# Accessing mz-rt feature data values

`featureValues,XCMSnExp()` : extract a `matrix` for feature values with
rows representing features and columns samples. Parameter `value` allows
to define which column from the
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
matrix should be returned. Multiple chromatographic peaks from the same
sample can be assigned to a feature. Parameter `method` allows to
specify the method to be used in such cases to chose from which of the
peaks the value should be returned. Parameter `msLevel` allows to choose
a specific MS level for which feature values should be returned (given
that features have been defined for that MS level).

`quantify,XCMSnExp()`: return the preprocessing results as an
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
object containing the feature abundances as assay matrix, the feature
definitions (returned by
[`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md))
as `rowData` and the phenotype information as `colData`. This is an
ideal container for further processing of the data. Internally, the
`featureValues()` method is used to extract the feature abundances,
parameters for that method can be passed to `quantify` with `...`.

## Usage

``` r
# S4 method for class 'XCMSnExp'
quantify(object, ...)

# S4 method for class 'XCMSnExp'
featureValues(
  object,
  method = c("medret", "maxint", "sum"),
  value = "into",
  intensity = "into",
  filled = TRUE,
  missing = NA,
  msLevel = integer()
)
```

## Arguments

- object:

  A
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object providing the feature definitions.

- ...:

  For
  [`quantify()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  additional parameters to be passed on to the \[featureValues()\`
  method.

- method:

  `character` specifying the method to resolve multi-peak mappings
  within the same sample, i.e. to define the *representative* peak for a
  feature in samples where more than one peak was assigned to the
  feature. If `"medret"`: select the peak closest to the median
  retention time of the feature. If `"maxint"`: select the peak yielding
  the largest signal. If `"sum"`: sum the values (only if `value` is
  `"into"` or `"maxo"`.

- value:

  `character` specifying the name of the column in `chromPeaks(object)`
  that should be returned. Defaults to `"into"` in which case the
  integrated peak area is returned. To get the index of the peak in the
  `chromPeaks(object)` matrix use `"index"`.

- intensity:

  `character` specifying the name of the column in the
  `chromPeaks(objects)` matrix containing the intensity value of the
  peak that should be used for the conflict resolution if
  `method = "maxint"`.

- filled:

  `logical(1)` specifying whether values for filled-in peaks should be
  returned or not. If `filled = FALSE`, an `NA` is returned in the
  matrix for the respective peak. See
  [`fillChromPeaks()`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)
  for details on peak filling.

- missing:

  how missing values should be reported. Allowed values are `NA` (the
  default), a `numeric` or `missing = "rowmin_half"`. The latter
  replaces any `NA` with half of the row's minimal (non-missing) value.

- msLevel:

  for `featureValues()`: `integer` defining the MS level(s) for which
  feature values should be returned. By default, values for features
  defined for all MS levels are returned.

## Value

For `featureValues()`: a `matrix` with feature values, columns
representing samples, rows features. The order of the features matches
the order found in the `featureDefinitions(object)` `DataFrame`. The
rownames of the `matrix` are the same than those of the
`featureDefinitions` `DataFrame`. `NA` is reported for features without
corresponding chromatographic peak in the respective sample(s).

For
[`quantify()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
a
[`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html)
representing the preprocessing results.

## See also

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for information on the data object.

[`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
to extract the `DataFrame` with the feature definitions.

[`featureChromatograms()`](https://sneumann.github.io/xcms/reference/featureChromatograms.md)
to extract ion chromatograms for each feature.

[`hasFeatures()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
to evaluate whether the `XCMSnExp` provides feature definitions.

## Author

Johannes Rainer
