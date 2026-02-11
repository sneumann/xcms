# Extract spectra associated with features

This function returns spectra associated with the identified features in
the input object. By default, spectra are returned for all features
(from all MS levels), but parameter `features` allows to specify/select
features for which spectra should be returned. Parameter `msLevel`
allows to define whether MS level 1 or 2 spectra should be returned. For
`msLevel = 1L` MS1 spectra within the retention time range of each
chromatographic peak (in that respective data file) associated with a
feature are returned. For `msLevel = 2L` MS2 spectra with a retention
time within the retention time range and their precursor m/z within the
m/z range of any chromatographic peak of a feature are returned. Thus,
only MS2 spectra for chromatographic peaks associated with the feature
and also measured in the sample in which the chromatographic was
identified are reported. By default, all spectra fulfilling the above
described condition are reported. This can be adapted with parameter
`method`. See the description of `method` in the
[`chromPeakSpectra()`](https://sneumann.github.io/xcms/reference/chromPeakSpectra.md)
documentation for more information. Internally, `featureSpectra()` uses
[`chromPeakSpectra()`](https://sneumann.github.io/xcms/reference/chromPeakSpectra.md)
to extract the feature's chromatographic peaks' spectra, thus any other
parameter for this function can be passed through `...`.

Note that with the default for parameter `skipFilled`
(`skipFilled = FALSE`) also gap-filled chromatographic peaks are
considered. Use `skipFilled = TRUE` to report only spectra for
**detected** peaks.

The information from `featureDefinitions` for each feature can be
included in the returned
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object using the `featureColumns` parameter. This is useful for keeping
details such as the median retention time (`rtmed`) or median m/z
(`mzmed`). The columns will retain their names as specified in the
`featureDefinitions` data, prefixed by `"feature_"` (e.g.,
`"feature_mzmed"`). Additionally, the *feature ID* (i.e., the row name
of the feature in the `featureDefinitions` data frame) is always added
as a metadata column named `"feature_id"`.

See also
[`chromPeakSpectra()`](https://sneumann.github.io/xcms/reference/chromPeakSpectra.md),
as it supports a similar parameter for including columns from the
chromatographic peaks in the returned spectra object. These parameters
can be used in combination to include information from both the
chromatographic peaks and the features in the returned
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html).
The *peak ID* (i.e., the row name of the peak in the `chromPeaks`
matrix) is added as a metadata column named `"chrom_peak_id"`.

## Usage

``` r
featureSpectra(object, ...)

# S4 method for class 'XcmsExperiment'
featureSpectra(
  object,
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  skipFilled = FALSE,
  return.type = c("Spectra", "List"),
  features = character(),
  featureColumns = c("rtmed", "mzmed"),
  ...
)

# S4 method for class 'XCMSnExp'
featureSpectra(
  object,
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  skipFilled = FALSE,
  return.type = c("MSpectra", "Spectra", "list", "List"),
  features = character(),
  ...
)
```

## Arguments

- object:

  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with feature defitions.

- ...:

  additional arguments to be passed along to
  [`chromPeakSpectra()`](https://sneumann.github.io/xcms/reference/chromPeakSpectra.md),
  such as `method` or `chromPeakColumns`.

- msLevel:

  `integer(1)` defining the MS level of the spectra that should be
  returned.

- expandRt:

  `numeric(1)` to expand the retention time range of each peak by a
  constant value on each side.

- expandMz:

  `numeric(1)` to expand the m/z range of each peak by a constant value
  on each side.

- ppm:

  `numeric(1)` to expand the m/z range of each peak (on each side) by a
  value dependent on the peak's m/z.

- skipFilled:

  `logical(1)` whether spectra for filled-in peaks should be reported or
  not. Defaults to `skipFilled = FALSE` thus also spectra for gap-filled
  chromatographic peaks are returned. Set to `skipFilled = TRUE` to get
  only spectra for detected chromatographic peaks.

- return.type:

  `character(1)` defining the type of result object that should be
  returned.

- features:

  `character`, `logical` or `integer` allowing to specify a subset of
  features in `featureDefinitions` for which spectra should be returned
  (providing either their ID, a logical vector same length than
  `nrow(featureDefinitions(x))` or their index in
  `featureDefinitions(x)`). This parameter is only supported for
  `return.type` being either `"Spectra"` or `"List"`.

- featureColumns:

  `character` vector with the names of the columns from
  `featureDefinitions` that should be added to the returned spectra
  object. The columns will be named as they are written in the
  `featureDefinitions` object with the prefix `"feature_`. Defaults to
  `c("mzmed", "rtmed")`.

## Value

The function returns either a
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
(for `return.type = "Spectra"`) or a `List` of `Spectra` (for
`return.type = "List"`). For the latter, the order and the length
matches parameter `features` (or if no `features` is defined the order
of the features in `featureDefinitions(object)`).

Spectra variables `"chrom_peak_id"` and `"feature_id"` define to which
chromatographic peak or feature each individual spectrum is associated
with.

## Author

Johannes Rainer
