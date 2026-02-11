# Identify overlapping features

`overlappingFeatures` identifies features that are overlapping or close
in the m/z - rt space.

## Usage

``` r
overlappingFeatures(x, expandMz = 0, expandRt = 0, ppm = 0)
```

## Arguments

- x:

  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with the features.

- expandMz:

  `numeric(1)` with the value to expand each feature (on each side) in
  m/z dimension before identifying overlapping features. The resulting
  `"mzmin"` for the feature is thus `mzmin - expandMz` and the `"mzmax"`
  `mzmax + expandMz`.

- expandRt:

  `numeric(1)` with the value to expand each feature (on each side) in
  retention time dimension before identifying overlapping features. The
  resulting `"rtmin"` for the feature is thus `rtmin - expandRt` and the
  `"rtmax"` `rtmax + expandRt`.

- ppm:

  `numeric(1)` to grow the m/z width of the feature by a relative value:
  `mzmin - mzmin * ppm / 2e6`, `mzmax + mzmax * ppm / 2e6`. Each feature
  is thus expanded in m/z dimension by ppm/2 on each side before
  identifying overlapping features.

## Value

`list` with indices of features (in
[`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md))
that are overlapping.

## Author

Johannes Rainer

## Examples

``` r

## Load a test data set with detected peaks
library(MSnbase)
data(faahko_sub)
## Update the path to the files for the local system
dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")

## Disable parallel processing for this example
register(SerialParam())

## Correspondence analysis
xdata <- groupChromPeaks(faahko_sub, param = PeakDensityParam(sampleGroups = c(1, 1, 1)))

## Identify overlapping features
overlappingFeatures(xdata)
#> list()

## Identify features that are separated on retention time by less than
## 2 minutes
overlappingFeatures(xdata, expandRt = 60)
#> [[1]]
#> [1] 5 6
#> 
#> [[2]]
#> [1] 31 32
#> 
```
