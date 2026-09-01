# Export data for use in MetaboAnalyst

Export the feature table for further analysis in the MetaboAnalyst
software (or the `MetaboAnalystR` R package).

## Usage

``` r
exportMetaboAnalyst(
  x,
  file = NULL,
  label,
  value = "into",
  digits = NULL,
  groupnames = FALSE,
  ...
)
```

## Arguments

- x:

  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with identified chromatographic peaks grouped across samples.

- file:

  `character(1)` defining the file name. If not specified, the `matrix`
  with the content is returned.

- label:

  either `character(1)` specifying the phenodata column in `x` defining
  the sample grouping or a vector with the same length than samples in
  `x` defining the group assignment of the samples.

- value:

  `character(1)` specifying the value to be returned for each feature.
  See
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  for more details.

- digits:

  `integer(1)` defining the number of significant digits to be used for
  numeric. The default `NULL` uses `getOption("digits")`. See
  [`format()`](https://rdrr.io/r/base/format.html) for more information.

- groupnames:

  `logical(1)` whether row names of the resulting matrix should be the
  feature IDs (`groupnames = FALSE`; default) or IDs that are composed
  of the m/z and retention time of the features (in the format
  `M<m/z>T<rt>` (`groupnames = TRUE`). See help of the
  [groupnames](https://sneumann.github.io/xcms/reference/groupnames-methods.md)
  function for details.

- ...:

  additional parameters to be passed to the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  function.

## Value

If `file` is not specified, the function returns the `matrix` in the
format supported by MetaboAnalyst.

## Author

Johannes Rainer
