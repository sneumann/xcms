# Filter features based on the dispersion ratio

The `DratioFilter` class and method enable users to filter features from
an `XcmsExperiment` or `SummarizedExperiment` object based on the
D-ratio or *dispersion ratio*. This is defined as the standard deviation
for QC samples divided by the standard deviation for biological test
samples, for each feature of the object (Broadhurst et al.).

This `filter` is part of the possible dispatch of the generic function
`filterFeatures`. Features *above* (`>`) the user-input threshold will
be removed from the entire dataset.

## Usage

``` r
DratioFilter(
  threshold = 0.5,
  qcIndex = integer(),
  studyIndex = integer(),
  na.rm = TRUE,
  mad = FALSE
)

# S4 method for class 'XcmsResult,DratioFilter'
filterFeatures(object, filter, ...)

# S4 method for class 'SummarizedExperiment,DratioFilter'
filterFeatures(object, filter, assay = 1)
```

## Arguments

- threshold:

  `numeric` value representing the threshold. Features with a D-ratio
  *strictly higher* (`>`) than this will be removed from the entire
  dataset.

- qcIndex:

  `integer` (or `logical`) vector corresponding to the indices of QC
  samples.

- studyIndex:

  `integer` (or `logical`) vector corresponding of the indices of study
  samples.

- na.rm:

  `logical` Indicates whether missing values (`NA`) should be removed
  prior to the calculations.

- mad:

  `logical` Indicates whether the *Median Absolute Deviation* (MAD)
  should be used instead of the standard deviation. This is suggested
  for non-gaussian distributed data.

- object:

  `XcmsExperiment` or `SummarizedExperiment`. For an `XcmsExperiment`
  object, the `featureValues(object)` will be evaluated, and for
  `Summarizedesxperiment` the `assay(object, assay)`. The object will be
  filtered.

- filter:

  The parameter object selecting and configuring the type of filtering.
  It can be one of the following classes:
  [`RsdFilter`](https://sneumann.github.io/xcms/reference/RsdFilter.md),
  `DratioFilter`,
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

For `DratioFilter`: a `DratioFilter` class. `filterFeatures` return the
input object minus the features that did not met the user input
threshold

## References

Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR,
Dunn WB. Guidelines and considerations for the use of system suitability
and quality control samples in mass spectrometry assays applied in
untargeted clinical metabolomic studies. Metabolomics. 2018;14(6):72.
doi:
[10.1007/s11306-018-1367-3](https://doi.org/10.1007/s11306-018-1367-3).
Epub 2018 May 18. PMID: 29805336; PMCID: PMC5960010.

## See also

Other Filter features in xcms:
[`BlankFlag()`](https://sneumann.github.io/xcms/reference/BlankFlag.md),
[`PercentMissingFilter()`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md),
[`RsdFilter()`](https://sneumann.github.io/xcms/reference/RsdFilter.md)

## Author

Philippine Louail
