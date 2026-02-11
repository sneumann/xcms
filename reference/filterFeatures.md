# Filtering of features based on conventional quality assessment

When dealing with metabolomics results, it is often necessary to filter
features based on certain criteria. These criteria are typically derived
from statistical formulas applied to full rows of data, where each row
represents a feature and its abundance of signal in each samples. The
`filterFeatures` function filters features based on these conventional
quality assessment criteria. Multiple types of filtering are implemented
and can be defined by the `filter` argument.

Supported `filter` arguments are:

- [`RsdFilter`](https://sneumann.github.io/xcms/reference/RsdFilter.md):
  Calculates the relative standard deviation (i.e. coefficient of
  variation) in abundance for each feature in QC (Quality Control)
  samples and filters them in the input object according to a provided
  threshold.

- [`DratioFilter`](https://sneumann.github.io/xcms/reference/DratioFilter.md):
  Computes the D-ratio or *dispersion ratio*, defined as the standard
  deviation in abundance for QC samples divided by the standard
  deviation for biological test samples, for each feature and filters
  them according to a provided threshold.

- [`PercentMissingFilter`](https://sneumann.github.io/xcms/reference/PercentMissingFilter.md):
  Determines the percentage of missing values for each feature in the
  various sample groups and filters them according to a provided
  threshold.

- [`BlankFlag`](https://sneumann.github.io/xcms/reference/BlankFlag.md):
  Identifies features where the mean abundance in test samples is lower
  than a specified multiple of the mean abundance of blank samples. This
  can be used to flag features that result from contamination in the
  solvent of the samples. A new column `possible_contaminants` is added
  to the `featureDefinitions` (`XcmsExperiment` object) or `rowData`
  (`SummarizedExperiment` object) reflecting this.

For specific examples, see the help pages of the individual parameter
classes listed above.

## Arguments

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
  or
  [`BlankFlag`](https://sneumann.github.io/xcms/reference/BlankFlag.md).

- assay:

  For filtering of `SummarizedExperiment` objects only. Indicates which
  assay the filtering will be based on. Note that the features for the
  entire object will be removed, but the computations are performed on a
  single assay. Default is 1, which means the first assay of the
  `object` will be evaluated.

- ...:

  Optional parameters. For `object` being an `XcmsExperiment`:
  parameters for the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  call.

## References

Broadhurst D, Goodacre R, Reinke SN, Kuligowski J, Wilson ID, Lewis MR,
Dunn WB. Guidelines and considerations for the use of system suitability
and quality control samples in mass spectrometry assays applied in
untargeted clinical metabolomic studies. Metabolomics. 2018;14(6):72.
doi:
[10.1007/s11306-018-1367-3](https://doi.org/10.1007/s11306-018-1367-3).
Epub 2018 May 18. PMID: 29805336; PMCID: PMC5960010.

## Author

Philippine Louail

## Examples

``` r
## See the vignettes for more detailed examples
library(MsExperiment)

## Load a test data set with features defined.
test_xcms <- loadXcmsData()
## Set up parameter to filter based on coefficient of variation. By setting
## the filter such as below, features that have a coefficient of variation
## superior to 0.3 in QC samples will be removed from the object `test_xcms`
## when calling the `filterFeatures` function.

rsd_filter <- RsdFilter(threshold = 0.3,
                        qcIndex = sampleData(test_xcms)$sample_type == "QC")

filtered_data_rsd <- filterFeatures(object = test_xcms, filter = rsd_filter)
#> 255 features were removed

## Set up parameter to filter based on D-ratio. By setting the filter such
## as below, features that have a D-ratio computed based on their abundance
## between QC and study samples superior to 0.5 will be removed from the
## object `test_xcms`.

dratio_filter <- DratioFilter(threshold = 0.5,
                 qcIndex = sampleData(test_xcms)$sample_type == "QC",
                 studyIndex = sampleData(test_xcms)$sample_type == "study")

filtered_data_dratio <- filterFeatures(object = test_xcms,
                                       filter = dratio_filter)
#> 237 features were removed

## Set up parameter to filter based on the percent of missing data.
## Parameter f should represent the sample group of samples, for which the
## percentage of missing values will be evaluated. As the setting is defined
## bellow, if a feature as less (or equal) to 30% missing values in one
## sample group, it will be kept in the `test_xcms` object.

missing_data_filter <- PercentMissingFilter(threshold = 30,
                                       f = sampleData(test_xcms)$sample_type)

filtered_data_missing <- filterFeatures(object = test_xcms,
                                        filter = missing_data_filter)
#> 2 features were removed

## Set up parameter to flag possible contaminants based on blank samples'
## abundance. By setting the filter such as below, features that have mean
## abundance ratio between blank(here use study as an example) and QC
## samples less than 2 will be marked as `TRUE` in an extra column named
## `possible_contaminants` in the `featureDefinitions` table of the object
## `test_xcms`.

filter <- BlankFlag(threshold = 2,
                    qcIndex = sampleData(test_xcms)$sample_type == "QC",
                    blankIndex = sampleData(test_xcms)$sample_type == "study")
filtered_xmse <- filterFeatures(test_xcms, filter)
#> -10 features were flagged
```
