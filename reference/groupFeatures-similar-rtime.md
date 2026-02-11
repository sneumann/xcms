# Compounding/feature grouping based on similar retention times

Group features based on similar retention time. This method is supposed
to be used as an initial *crude* grouping of features based on the
median retention time of all their chromatographic peaks. All features
with a difference in their retention time which is `<=` parameter
`diffRt` of the parameter object are grouped together. If a column
`"feature_group"` is found in
[`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
this is further sub-grouped by this method.

See
[`MsFeatures::SimilarRtimeParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-rtime.html)
in `MsFeatures` for more details.

## Usage

``` r
# S4 method for class 'XcmsResult,SimilarRtimeParam'
groupFeatures(object, param, msLevel = 1L, ...)
```

## Arguments

- object:

  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object containing also correspondence results.

- param:

  `SimilarRtimeParam` object with the settings for the method. See
  [`MsFeatures::SimilarRtimeParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-rtime.html)
  for details and options.

- msLevel:

  `integer(1)` defining the MS level on which the features should be
  grouped.

- ...:

  passed to the
  [`groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
  function on numeric values.

## Value

the input object with feature groups added (i.e. in column
`"feature_group"` of its `featureDefinitions` data frame.

## See also

Other feature grouping methods:
[`groupFeatures-abundance-correlation`](https://sneumann.github.io/xcms/reference/groupFeatures-abundance-correlation.md),
[`groupFeatures-eic-similarity`](https://sneumann.github.io/xcms/reference/groupFeatures-eic-similarity.md)

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

## Group features based on similar retention time (i.e. difference <= 2 seconds)
xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 2))

## Feature grouping get added to the featureDefinitions in column "feature_group"
head(featureDefinitions(xodg_grp)$feature_group)
#> [1] "FG.004" "FG.009" "FG.003" "FG.004" "FG.010" "FG.011"

table(featureDefinitions(xodg_grp)$feature_group)
#> 
#> FG.001 FG.002 FG.003 FG.004 FG.005 FG.006 FG.007 FG.008 FG.009 FG.010 FG.011 
#>      2      2      3      3      3      2      2      2      1      1      1 
#> FG.012 FG.013 FG.014 FG.015 FG.016 FG.017 FG.018 FG.019 FG.020 FG.021 FG.022 
#>      1      1      1      1      1      1      1      1      1      1      1 
#> FG.023 FG.024 FG.025 FG.026 FG.027 FG.028 FG.029 FG.030 FG.031 FG.032 FG.033 
#>      1      1      1      1      1      1      1      1      1      1      1 
#> FG.034 FG.035 FG.036 
#>      1      1      1 
length(unique(featureDefinitions(xodg_grp)$feature_group))
#> [1] 36

## Using an alternative groupiing method that creates larger groups
xodg_grp <- groupFeatures(xodg,
    param = SimilarRtimeParam(diffRt = 2, groupFun = MsCoreUtils::group))

length(unique(featureDefinitions(xodg_grp)$feature_group))
#> [1] 35
```
