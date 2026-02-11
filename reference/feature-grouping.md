# Compounding of LC-MS features

Feature *compounding* aims at identifying and grouping LC-MS features
representing different ions or adducts (including isotopes) of the same
originating compound. The
[MsFeatures](https://bioconductor.org/packages/MsFeatures) package
provides a general framework and functionality to group features based
on different properties. The `groupFeatures` methods for
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
or
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
objects implemented in `xcms` extend these to enable the *compounding*
of LC-MS data considering also e.g. feature peak shaped. Note that these
functions simply define feature groups but don't actually *aggregate* or
combine the features.

See
[`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html)
for an overview on the general feature grouping concept as well as
details on the individual settings and parameters.

The available options for `groupFeatures` on `xcms` preprocessing
results (i.e. on `XcmsExperiment` or `XCMSnExp` objects after
correspondence analysis with
[`groupChromPeaks()`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md))
are:

- Grouping by similar retention times:
  [`groupFeatures-similar-rtime()`](https://sneumann.github.io/xcms/reference/groupFeatures-similar-rtime.md).

- Grouping by similar feature values across samples:
  [`MsFeatures::AbundanceSimilarityParam()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures-similar-abundance.html).

- Grouping by similar peak shape of extracted ion chromatograms:
  [`EicSimilarityParam()`](https://sneumann.github.io/xcms/reference/groupFeatures-eic-similarity.md).

An ideal workflow grouping features should sequentially perform the
above methods (in the listed order).

Compounded feature groups can be accessed with the `featureGroups`
function.

## Usage

``` r
# S4 method for class 'XcmsResult'
featureGroups(object)

# S4 method for class 'XcmsResult'
featureGroups(object) <- value
```

## Arguments

- object:

  an
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with LC-MS pre-processing results.

- value:

  for `featureGroups<-`: replacement for the feature groups in `object`.
  Has to be of length 1 or length equal to the number of features in
  `object`.

## See also

[`plotFeatureGroups()`](https://sneumann.github.io/xcms/reference/plotFeatureGroups.md)
for visualization of grouped features.

## Author

Johannes Rainer, Mar Garcia-Aloy, Vinicius Veri Hernandes
