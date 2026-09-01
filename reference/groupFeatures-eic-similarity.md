# Compounding/feature grouping based on similarity of extracted ion chromatograms

Features from the same originating compound are expected to share their
elution pattern (i.e. chromatographic peak shape) with it. Thus, this
methods allows to group features based on similarity of their extracted
ion chromatograms (EICs). The similarity calculation is performed
separately for each sample with the similarity score being aggregated
across samples for the final generation of the similarity matrix on
which the grouping (considering parameter `threshold`) will be
performed.

The
[`MSnbase::compareChromatograms()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function is used for similarity calculation which by default calculates
the Pearson's correlation coefficient. The settings for
[`compareChromatograms()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
can be specified with parameters `ALIGNFUN`, `ALIGNFUNARGS`, `FUN` and
`FUNARGS`. `ALIGNFUN` defaults to `alignRt` and is the function used to
*align* the chromatograms before comparison. For information and
parameters of
[`alignRt()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
see the documentation for
[`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html).
`ALIGNFUNARGS` allows to specify additional arguments for the `ALIGNFUN`
function. It defaults to
`ALIGNFUNARGS = list(tolerance = 0, method = "closest")` which ensures
that data points from the same spectrum (scan, i.e. with the same
retention time) are compared between the EICs from the same sample.
Parameter `FUN` defines the function to calculate the similarity score
and defaults to `FUN = cor` and `FUNARGS` allows to pass additional
arguments to this function (defaults to
`FUNARGS = list(use = "pairwise.complete.obs")`. See also
[`MSnbase::compareChromatograms()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
for more information.

The grouping of features based on the EIC similarity matrix is performed
with the function specified with parameter `groupFun` which defaults to
`groupFun = groupSimilarityMatrix` which groups all rows (features) in
the similarity matrix with a similarity score larger than `threshold`
into the same cluster. This creates clusters of features in which
**all** features have a similarity score `>= threshold` with **any**
other feature in that cluster. See
[`MsFeatures::groupSimilarityMatrix()`](https://rdrr.io/pkg/MsFeatures/man/groupSimilarityMatrix.html)
for details. Additional parameters to that function can be passed with
the `...` argument.

This feature grouping should be called **after** an initial feature
grouping by retention time (see
[`MsFeatures::SimilarRtimeParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-rtime.html)).
The feature groups defined in columns `"feature_group"` of
`featureDefinitions(object)` (for features matching `msLevel`) will be
used and refined by this method. Features with a value of `NA` in
`featureDefinitions(object)$feature_group` will be skipped/not
considered for feature grouping.

## Usage

``` r
EicSimilarityParam(
  threshold = 0.9,
  n = 1,
  onlyPeak = TRUE,
  value = c("maxo", "into"),
  groupFun = groupSimilarityMatrix,
  ALIGNFUN = alignRt,
  ALIGNFUNARGS = list(tolerance = 0, method = "closest"),
  FUN = cor,
  FUNARGS = list(use = "pairwise.complete.obs"),
  ...
)

# S4 method for class 'XcmsResult,EicSimilarityParam'
groupFeatures(object, param, msLevel = 1L)
```

## Arguments

- threshold:

  `numeric(1)` with the minimal required similarity score to group
  featues. This is passed to the `groupFun` function.

- n:

  `numeric(1)` defining the total number of samples per feature group on
  which this similarity calculation should be performed. This value is
  rounded up to the next larger integer value.

- onlyPeak:

  `logical(1)` whether the correlation should be performed only on the
  signals within the identified chromatographic peaks
  (`onlyPeak = TRUE`, default) or all the signal from the extracted ion
  chromatogram.

- value:

  `character(1)` defining whether samples should be grouped based on the
  sum of the maximal peak intensity (`value = "maxo"`, the default) or
  the integrated peak area (`value = "into"`) for a feature.

- groupFun:

  `function` defining the function to be used to group rows based on a
  pairwise similarity matrix. Defaults to
  [`MsFeatures::groupSimilarityMatrix()`](https://rdrr.io/pkg/MsFeatures/man/groupSimilarityMatrix.html).

- ALIGNFUN:

  `function` defining the function to be used to *align* chromatograms
  prior similarity calculation. Defaults to `ALIGNFUN = alignRt`. See
  documentation of
  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
  and
  [`MSnbase::compareChromatograms()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  for more information.

- ALIGNFUNARGS:

  **named** `list` with arguments for `ALIGNFUN`. Defaults to
  `ALIGNFUNARGS = list(tolerance = 0, method = "closest")`.

- FUN:

  `function` defining the function to be used to calculate a similarity
  between (aligned) chromatograms. Defaults to `FUN = cor`. See
  [`cor()`](https://rdrr.io/r/stats/cor.html) and
  [`MSnbase::compareChromatograms()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  for more information.

- FUNARGS:

  **named** `list` with arguments for `FUN`. Defaults to
  `FUN = list(use = "pairwise.complete.obs")`.

- ...:

  for `EicSimilarityParam`: additional arguments to be passed to
  `groupFun` and `featureChromatograms` (such as `expandRt` to expand
  the retention time range of each feature).

- object:

  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with LC-MS pre-processing results.

- param:

  `EicSimilarityParam` object with the settings for the method.

- msLevel:

  `integer(1)` defining the MS level on which the features should be
  grouped.

## Value

input object with feature groups added (i.e. in column `"feature_group"`
of its `featureDefinitions` data frame.

## Note

At present the
[`featureChromatograms()`](https://sneumann.github.io/xcms/reference/featureChromatograms.md)
function is used to extract the EICs for each feature, which does
however use one m/z and rt range for each feature and the EICs do thus
not exactly represent the identified chromatographic peaks of each
sample (i.e. their specific m/z and retention time ranges).

While being possible to be performed on the full data set without prior
feature grouping, this is not suggested for the following reasons: I)
the selection of the top `n` samples with the highest signal for the
*feature group* will be biased by very abundant compounds as this is
performed on the full data set (i.e. the samples with the highest
overall intensities are used for correlation of all features) and II) it
is computationally much more expensive because a pairwise correlation
between all features has to be performed.

It is also suggested to perform the correlation on a subset of samples
per feature with the highest intensities of the peaks (for that feature)
although it would also be possible to run the correlation on all samples
by setting `n` equal to the total number of samples in the data set. EIC
correlation should however be performed ideally on samples in which the
original compound is highly abundant to avoid correlation of missing
values or noisy peak shapes as much as possible.

By default also the signal which is outside identified chromatographic
peaks is excluded from the correlation.

## See also

feature-grouping for a general overview.

Other feature grouping methods:
[`groupFeatures-abundance-correlation`](https://sneumann.github.io/xcms/reference/groupFeatures-abundance-correlation.md),
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

## Performing a feature grouping based on EIC similarities on a single
## sample
xodg_grp <- groupFeatures(xodg, param = EicSimilarityParam(n = 1))

table(featureDefinitions(xodg_grp)$feature_group)
#> 
#> FG.001 FG.002 FG.003 FG.004 FG.005 FG.006 FG.007 FG.008 FG.009 FG.010 FG.011 
#>      2      2      2      3      4      2      3      2      1      1      1 
#> FG.012 FG.013 FG.014 FG.015 FG.016 FG.017 FG.018 FG.019 FG.020 FG.021 FG.022 
#>      1      1      1      1      1      1      1      1      1      1      1 
#> FG.023 FG.024 FG.025 FG.026 FG.027 FG.028 FG.029 FG.030 FG.031 FG.032 FG.033 
#>      1      1      1      1      1      1      1      1      1      1      1 
#> FG.034 FG.035 
#>      1      1 

## Usually it is better to perform this correlation on pre-grouped features
## e.g. based on similar retention time.
xodg_grp <- groupFeatures(xodg, param = SimilarRtimeParam(diffRt = 4))
xodg_grp <- groupFeatures(xodg_grp, param = EicSimilarityParam(n = 1))

table(featureDefinitions(xodg_grp)$feature_group)
#> 
#> FG.001.001 FG.002.001 FG.003.001 FG.003.002 FG.004.001 FG.004.002 FG.005.001 
#>          2          2          2          1          3          1          3 
#> FG.006.001 FG.006.002 FG.007.001 FG.008.001 FG.009.001 FG.009.002 FG.010.001 
#>          1          1          2          2          1          1          2 
#> FG.011.001 FG.012.001 FG.013.001 FG.014.001 FG.015.001 FG.016.001 FG.017.001 
#>          1          1          1          1          1          1          1 
#> FG.018.001 FG.019.001 FG.020.001 FG.021.001 FG.022.001 FG.023.001 FG.024.001 
#>          1          1          1          1          1          1          1 
#> FG.025.001 FG.026.001 FG.027.001 FG.028.001 FG.029.001 FG.030.001 FG.031.001 
#>          1          1          1          1          1          1          1 
#> FG.032.001 FG.033.001 
#>          1          1 
```
