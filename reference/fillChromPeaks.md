# Gap Filling

Gap filling integrate signal in the m/z-rt area of a feature (i.e., a
chromatographic peak group) for samples in which no chromatographic peak
for this feature was identified and add it to the
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
matrix. Such *filled-in* peaks are indicated with a `TRUE` in column
`"is_filled"` in the result object's
[`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
data frame.

The method for gap filling along with its settings can be defined with
the `param` argument. Two different approaches are available:

- `param = FillChromPeaksParam()`: the default of the original `xcms`
  code. Signal is integrated from the m/z and retention time range as
  defined in the
  [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  data frame, i.e. from the `"rtmin"`, `"rtmax"`, `"mzmin"` and
  `"mzmax"`. This method is not suggested as it underestimates the
  actual peak area and it is also not available for `object` being an
  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  object. See details below for more information and settings for this
  method.

- `param = ChromPeakAreaParam()`: the area from which the signal for a
  feature is integrated is defined based on the feature's
  chromatographic peak areas. The m/z range is by default defined as the
  the lower quartile of chromatographic peaks' `"mzmin"` value to the
  upper quartile of the chromatographic peaks' `"mzmax"` values. The
  retention time range for the area is defined analogously.
  Alternatively, by setting `mzmin = median`, `mzmax = median`,
  `rtmin = median` and `rtmax = median` in `ChromPeakAreaParam`, the
  median `"mzmin"`, `"mzmax"`, `"rtmin"` and `"rtmax"` values from all
  detected chromatographic peaks of a feature would be used instead.
  Parameter `minMzWidthPpm` allows in addition to define a minimal
  guaranteed m/z width expressed in ppm of the features' m/z and
  centered around it. The default is `minMzWidthPpm = 0.0`. With a
  `minMzWidthPpm` \> 0, the lower m/z boundary for a feature is defined
  as the smaller value from the m/z derived from its chromatographic
  peaks' `"mzmin"`, and the feature's m/z *minus* `minMzWidthPpm / 2`
  ppm of its m/z. The upper m/z boundary is determined in the same way.
  In contrast to the `FillChromPeaksParam` approach this method uses
  (all) identified chromatographic peaks of a feature to define the area
  from which the signal should be integrated.

## Usage

``` r
fillChromPeaks(object, param, ...)

# S4 method for class 'XcmsExperiment,ChromPeakAreaParam'
fillChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

FillChromPeaksParam(
  expandMz = 0,
  expandRt = 0,
  ppm = 0,
  fixedMz = 0,
  fixedRt = 0
)

ChromPeakAreaParam(
  mzmin = function(z, na.rm = TRUE) quantile(z, probs = 0.25, names = FALSE, na.rm =
    na.rm),
  mzmax = function(z, na.rm = TRUE) quantile(z, probs = 0.75, names = FALSE, na.rm =
    na.rm),
  rtmin = function(z, na.rm = TRUE) quantile(z, probs = 0.25, names = FALSE, na.rm =
    na.rm),
  rtmax = function(z, na.rm = TRUE) quantile(z, probs = 0.75, names = FALSE, na.rm =
    na.rm),
  minMzWidthPpm = 0
)

# S4 method for class 'XCMSnExp,FillChromPeaksParam'
fillChromPeaks(object, param, msLevel = 1L, BPPARAM = bpparam())

# S4 method for class 'XCMSnExp,ChromPeakAreaParam'
fillChromPeaks(object, param, msLevel = 1L, BPPARAM = bpparam())

# S4 method for class 'XCMSnExp,missing'
fillChromPeaks(object, param, BPPARAM = bpparam(), msLevel = 1L)
```

## Arguments

- object:

  `XcmsExperiment` or `XCMSnExp` object with identified and grouped
  chromatographic peaks.

- param:

  `ChromPeakAreaParam` or `FillChromPeaksParam` object defining which
  approach should be used (see details section).

- ...:

  currently ignored.

- msLevel:

  `integer(1)` defining the MS level on which peak filling should be
  performed (defaults to `msLevel = 1L`). Only peak filling on one MS
  level at a time is supported, to fill in peaks for MS level 1 and 2
  run first using `msLevel = 1` and then (on the returned result object)
  again with `msLevel = 2`.

- chunkSize:

  For `fillChromPeaks` if `object` is an `XcmsExperiment`: `integer(1)`
  defining the number of files (samples) that should be loaded into
  memory and processed at the same time. This setting thus allows to
  balance between memory demand and speed (due to parallel processing).
  Because parallel processing can only performed on the subset of data
  currently loaded into memory in each iteration, the value for
  `chunkSize` should match the defined parallel setting setup. Using a
  parallel processing setup using 4 CPUs (separate processes) but using
  `chunkSize = `1`will not perform any parallel processing, as only the data from one sample is loaded in memory at a time. On the other hand, setting`chunkSize\`
  to the total number of samples in an experiment will load the full MS
  data into memory and will thus in most settings cause an out-of-memory
  error.

- BPPARAM:

  Parallel processing settings.

- expandMz:

  for `FillChromPeaksParam`: `numeric(1)` defining the value by which
  the mz width of peaks should be expanded. Each peak is expanded in mz
  direction by `expandMz *` their original m/z width. A value of `0`
  means no expansion, a value of `1` grows each peak by `1 *` the m/z
  width of the peak resulting in peaks with twice their original size in
  m/z direction (expansion by half m/z width to both sides).

- expandRt:

  for `FillChromPeaksParam`: `numeric(1)`, same as `expandMz` but for
  the retention time width.

- ppm:

  for `FillChromPeaksParam`: `numeric(1)` optionally specifying a *ppm*
  by which the m/z width of the peak region should be expanded. For
  peaks with an m/z width smaller than
  `mean(c(mzmin, mzmax)) * ppm / 1e6`, the `mzmin` will be replaced by
  `mean(c(mzmin, mzmax)) - (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`
  `mzmax` by
  `mean(c(mzmin, mzmax)) + (mean(c(mzmin, mzmax)) * ppm / 2 / 1e6)`.
  This is applied before eventually expanding the m/z width using the
  `expandMz` parameter.

- fixedMz:

  for `FillChromPeaksParam`: `numeric(1)` defining a constant factor by
  which the m/z width of each feature is to be expanded. The m/z width
  is expanded on both sides by `fixedMz` (i.e. `fixedMz` is subtracted
  from the lower m/z and added to the upper m/z). This expansion is
  applied *after* `expandMz` and `ppm`.

- fixedRt:

  for `FillChromPeaksParam`: `numeric(1)` defining a constant factor by
  which the retention time width of each factor is to be expanded. The
  rt width is expanded on both sides by `fixedRt` (i.e. `fixedRt` is
  subtracted from the lower rt and added to the upper rt). This
  expansion is applied *after* `expandRt`.

- mzmin:

  `function` to be applied to values in the `"mzmin"` column of all
  chromatographic peaks of a feature to define the lower m/z value of
  the area from which signal for the feature should be integrated.
  Defaults to `mzmin = function(z) quantile(z, probs = 0.25)` hence
  using the 25% quantile of all values.

- mzmax:

  `function` to be applied to values in the `"mzmax"` column of all
  chromatographic peaks of a feature to define the upper m/z value of
  the area from which signal for the feature should be integrated.
  Defaults to `mzmax = function(z) quantile(z, probs = 0.75)` hence
  using the 75% quantile of all values.

- rtmin:

  `function` to be applied to values in the `"rtmin"` column of all
  chromatographic peaks of a feature to define the lower rt value of the
  area from which signal for the feature should be integrated. Defaults
  to `rtmin = function(z) quantile(z, probs = 0.25)` hence using the 25%
  quantile of all values.

- rtmax:

  `function` to be applied to values in the `"rtmax"` column of all
  chromatographic peaks of a feature to define the upper rt value of the
  area from which signal for the feature should be integrated. Defaults
  to `rtmax = function(z) quantile(z, probs = 0.75)` hence using the 75%
  quantile of all values.

- minMzWidthPpm:

  For `ChromPeakAreaParam()`: `numeric(1)` defining the minimal
  guaranteed m/z width (expressed in ppm of the feature's m/z) that will
  be used to integrate signal from (default `minMzWidthPpm = 0.0`). See
  documentation of `ChromPeakAreaParam()` for more information.

## Value

An
[XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
or `XCMSnExp` object with previously missing chromatographic peaks for
features filled into its
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
matrix.

The `FillChromPeaksParam()` function returns a `FillChromPeaksParam`
object.

## Details

After correspondence (i.e. grouping of chromatographic peaks across
samples) there will always be features (peak groups) that do not include
peaks from every sample. The `fillChromPeaks` method defines intensity
values for such features in the missing samples by integrating the
signal in the m/z-rt region of the feature. Two different approaches to
define this region are available: with `ChromPeakAreaParam` the region
is defined based on the detected **chromatographic peaks** of a feature,
while with `FillChromPeaksParam` the region is defined based on the m/z
and retention times of the **feature** (which represent the m/z and
retentention times of the apex position of the associated
chromatographic peaks). For the latter approach various parameters are
available to increase the area from which signal is to be integrated,
either by a constant value (`fixedMz` and `fixedRt`) or by a
feature-relative amount (`expandMz` and `expandRt`).

Adjusted retention times will be used if available.

Based on the peak finding algorithm that was used to identify the
(chromatographic) peaks, different internal functions are used to
guarantee that the integrated peak signal matches as much as possible
the peak signal integration used during the peak detection. For peaks
identified with the
[`matchedFilter()`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
method, signal integration is performed on the *profile matrix*
generated with the same settings used also during peak finding (using
the same `bin` size for example). For direct injection data and peaks
identified with the `MSW` algorithm signal is integrated only along the
mz dimension. For all other methods the complete (raw) signal within the
area is used.

## Note

The reported `"mzmin"`, `"mzmax"`, `"rtmin"` and `"rtmax"` for the
filled peaks represents the actual MS area from which the signal was
integrated.

No peak is filled in if no signal was present in a file/sample in the
respective mz-rt area. These samples will still show a `NA` in the
matrix returned by the
[`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
method.

## See also

[`groupChromPeaks()`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md)
for methods to perform the correspondence.

[featureArea](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
for the function to define the m/z-retention time region for each
feature.

## Author

Johannes Rainer

## Examples

``` r

## Load a test data set with identified chromatographic peaks
library(xcms)
library(MsExperiment)
res <- loadXcmsData("faahko_sub2")

## Disable parallel processing for this example
register(SerialParam())

## Perform the correspondence. We assign all samples to the same group.
res <- groupChromPeaks(res,
    param = PeakDensityParam(sampleGroups = rep(1, length(res))))

## For how many features do we lack an integrated peak signal?
sum(is.na(featureValues(res)))
#> [1] 26

## Filling missing peak data using the peak area from identified
## chromatographic peaks.
res <- fillChromPeaks(res, param = ChromPeakAreaParam())

## Alternatively, force a minimal guaranteed m/z width for the regions
## to integrate signal from.
res <- fillChromPeaks(res, param = ChromPeakAreaParam(minMzWidthPpm = 10))

## How many missing values do we have after peak filling?
sum(is.na(featureValues(res)))
#> [1] 2

## Get the peaks that have been filled in:
fp <- chromPeaks(res)[chromPeakData(res)$is_filled, ]
head(fp)
#>          mz mzmin mzmax       rt    rtmin    rtmax     into intb   maxo sn
#> CP249 286.2 286.2 286.2 3252.556 3236.907 3274.857  1288521   NA  78048 NA
#> CP250 380.1 380.1 380.1 3193.087 3132.054 3216.561  2094879   NA  35928 NA
#> CP251 447.2 447.2 447.2 3858.193 3828.068 3914.140  2156784   NA  61256 NA
#> CP252 497.2 497.2 497.2 3382.447 3347.236 3438.785 10344783   NA 310912 NA
#> CP253 510.2 510.2 510.2 3750.211 3741.213 3799.898  1927180   NA  61368 NA
#> CP254 531.2 531.2 531.2 3340.193 3328.848 3383.620  7589092   NA 254464 NA
#>       sample
#> CP249      1
#> CP250      1
#> CP251      1
#> CP252      1
#> CP253      1
#> CP254      1

## Get the process history step along with the parameters used to perform
## The peak filling:
ph <- processHistory(res, type = "Missing peak filling")[[1]]
ph
#> Object of class "XProcessHistory"
#>  type: Missing peak filling 
#>  date: Thu Mar 19 15:02:32 2026 
#>  info:  
#>  fileIndex: 1,2,3 
#>  Parameter class: ChromPeakAreaParam 
#>  MS level(s) 1 

## The parameter class:
ph@param
#> Object of class:  ChromPeakAreaParam 
#>  Parameters:
#>  - rtmin: function (z, na.rm = TRUE) 
#> quantile(z, probs = 0.25, names = FALSE, na.rm = na.rm)
#> <environment: 0x56193c3e4a48>
#>  - rtmax: function (z, na.rm = TRUE) 
#> quantile(z, probs = 0.75, names = FALSE, na.rm = na.rm)
#> <environment: 0x56193c3e4a48>
#>  - mzmin: function (z, na.rm = TRUE) 
#> quantile(z, probs = 0.25, names = FALSE, na.rm = na.rm)
#> <environment: 0x56193c3e4a48>
#>  - mzmax: function (z, na.rm = TRUE) 
#> quantile(z, probs = 0.75, names = FALSE, na.rm = na.rm)
#> <environment: 0x56193c3e4a48>
#>  - minMzWidthPpm: [1] 0

## It is also possible to remove filled-in peaks:
res <- dropFilledChromPeaks(res)

sum(is.na(featureValues(res)))
#> [1] 26
```
