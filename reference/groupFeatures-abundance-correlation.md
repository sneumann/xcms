# Compounding/feature grouping based on similarity of abundances across samples

Features from the same originating compound are expected to have similar
intensities across samples. This method thus groups features based on
similarity of abundances (i.e. *feature values*) across samples in a
data set. See also
[`MsFeatures::AbundanceSimilarityParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-abundance.html)
for additional information and details.

This help page lists parameters specific for `xcms` result objects (i.e.
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
and
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
objects). Documentation of the parameters for the similarity calculation
is available in the
[`MsFeatures::AbundanceSimilarityParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-abundance.html)
help page in the *MsFeatures* package.

## Usage

``` r
# S4 method for class 'XcmsResult,AbundanceSimilarityParam'
groupFeatures(
  object,
  param,
  msLevel = 1L,
  method = c("medret", "maxint", "sum"),
  value = "into",
  intensity = "into",
  filled = TRUE,
  ...
)
```

## Arguments

- object:

  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object containing LC-MS pre-processing results.

- param:

  `AbudanceSimilarityParam` object with the settings for the method. See
  [`MsFeatures::AbundanceSimilarityParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-abundance.html)
  for details on the grouping method and its parameters.

- msLevel:

  `integer(1)` defining the MS level on which the features should be
  grouped.

- method:

  `character(1)` passed to the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  call. See
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  for details. Defaults to `method = "medret"`.

- value:

  `character(1)` passed to the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  call. See
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  for details. Defaults to `value = "into"`.

- intensity:

  `character(1)` passed to the
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  call. See
  [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  for details. Defaults to `intensity = "into"`.

- filled:

  `logical(1)` whether filled-in values should be included in the
  correlation analysis. Defaults to `filled = TRUE`.

- ...:

  additional parameters passed to the
  [`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
  method for `matrix`.

## Value

input object with feature group definitions added to (or updated in) a
column `"feature_group"` in its `featureDefinitions` data frame.

## See also

feature-grouping for a general overview.

Other feature grouping methods:
[`groupFeatures-eic-similarity`](https://sneumann.github.io/xcms/reference/groupFeatures-eic-similarity.md),
[`groupFeatures-similar-rtime`](https://sneumann.github.io/xcms/reference/groupFeatures-similar-rtime.md)

## Author

Johannes Rainer

## Examples

``` r

library(MsFeatures)
library(MsExperiment)
## Load a test data set with detected peaks
faahko_sub <- loadXcmsData("faahko_sub2")

## Disable parallel processing for this example
register(SerialParam())

## Group chromatographic peaks across samples
xodg <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = rep(1, 3)))

## Group features based on correlation of feature values (integrated
## peak area) across samples. Note that there are many missing values
## in the feature value which influence grouping of features in the present
## data set.
xodg_grp <- groupFeatures(xodg,
    param = AbundanceSimilarityParam(threshold = 0.8))
table(featureDefinitions(xodg_grp)$feature_group)
#> 
#> FG.001 FG.002 FG.003 FG.004 FG.005 FG.006 FG.007 FG.008 
#>      8      8      8      6     12      3      1      1 

## Group based on the maximal peak intensity per feature
xodg_grp <- groupFeatures(xodg,
    param = AbundanceSimilarityParam(threshold = 0.8, value = "maxo"))
table(featureDefinitions(xodg_grp)$feature_group)
#> 
#> FG.001 FG.002 FG.003 FG.004 FG.005 FG.006 FG.007 FG.008 
#>      8      8      8      6     12      3      1      1 
```
