# Two-step centWave peak detection considering also isotopes

This method performs a two-step centWave-based chromatographic peak
detection: in a first centWave run peaks are identified for which then
the location of their potential isotopes in the mz-retention time is
predicted. A second centWave run is then performed on these *regions of
interest* (ROIs). The final list of chromatographic peaks comprises all
non-overlapping peaks from both centWave runs.

The `findChromPeaks,OnDiskMSnExp,CentWavePredIsoParam()` method performs
a two-step centWave-based chromatographic peak detection on all samples
from an `OnDiskMSnExp` object. `OnDiskMSnExp` objects encapsule all
experiment specific data and load the spectra data (mz and intensity
values) on the fly from the original files applying also all eventual
data manipulations.

## Usage

``` r
CentWavePredIsoParam(
  ppm = 25,
  peakwidth = c(20, 50),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1L,
  mzdiff = -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = numeric(),
  extendLengthMSW = FALSE,
  verboseBetaColumns = FALSE,
  snthreshIsoROIs = 6.25,
  maxCharge = 3,
  maxIso = 5,
  mzIntervalExtension = TRUE,
  polarity = "unknown"
)

# S4 method for class 'OnDiskMSnExp,CentWavePredIsoParam'
findChromPeaks(
  object,
  param,
  BPPARAM = bpparam(),
  return.type = "XCMSnExp",
  msLevel = 1L,
  ...
)
```

## Arguments

- ppm:

  `numeric(1)` defining the maximal tolerated m/z deviation in
  consecutive scans in parts per million (ppm) for the initial ROI
  definition.

- peakwidth:

  `numeric(2)` with the expected approximate peak width in
  chromatographic space. Given as a range (min, max) in seconds.

- snthresh:

  `numeric(1)` defining the signal to noise ratio cutoff.

- prefilter:

  `numeric(2)`: `c(k, I)` specifying the prefilter step for the first
  analysis step (ROI detection). Mass traces are only retained if they
  contain at least `k` peaks with intensity `>= I`.

- mzCenterFun:

  Name of the function to calculate the m/z center of the
  chromatographic peak. Allowed are: `"wMean"`: intensity weighted mean
  of the peak's m/z values, `"mean"`: mean of the peak's m/z values,
  `"apex"`: use the m/z value at the peak apex, `"wMeanApex3"`:
  intensity weighted mean of the m/z value at the peak apex and the m/z
  values left and right of it and `"meanApex3"`: mean of the m/z value
  of the peak apex and the m/z values left and right of it.

- integrate:

  Integration method. For `integrate = 1` peak limits are found through
  descent on the mexican hat filtered data, for `integrate = 2` the
  descent is done on the real data. The latter method is more accurate
  but prone to noise, while the former is more robust, but less exact.

- mzdiff:

  `numeric(1)` representing the minimum difference in m/z dimension
  required for peaks with overlapping retention times; can be negative
  to allow overlap. During peak post-processing, peaks defined to be
  overlapping are reduced to the one peak with the largest signal.

- fitgauss:

  `logical(1)` whether or not a Gaussian should be fitted to each peak.
  This affects mostly the retention time position of the peak.

- noise:

  `numeric(1)` allowing to set a minimum intensity required for
  centroids to be considered in the first analysis step (centroids with
  intensity `< noise` are omitted from ROI detection).

- verboseColumns:

  `logical(1)` whether additional peak meta data columns should be
  returned.

- roiList:

  An optional list of regions-of-interest (ROI) representing detected
  mass traces. If ROIs are submitted the first analysis step is omitted
  and chromatographic peak detection is performed on the submitted ROIs.
  Each ROI is expected to have the following elements specified: `scmin`
  (start scan index), `scmax` (end scan index), `mzmin` (minimum m/z),
  `mzmax` (maximum m/z), `length` (number of scans), `intensity` (summed
  intensity). Each ROI should be represented by a `list` of elements or
  a single row `data.frame`.

- firstBaselineCheck:

  `logical(1)`. If `TRUE` continuous data within regions of interest is
  checked to be above the first baseline. In detail, a first rough
  estimate of the noise is calculated and peak detection is performed
  only in regions in which multiple sequential signals are higher than
  this first estimated baseline/noise level.

- roiScales:

  Optional numeric vector with length equal to `roiList` defining the
  scale for each region of interest in `roiList` that should be used for
  the centWave-wavelets.

- extendLengthMSW:

  Option to force centWave to use all scales when running centWave
  rather than truncating with the EIC length. Uses the "open" method to
  extend the EIC to a integer base-2 length prior to being passed to
  `convolve` rather than the default "reflect" method. See
  https://github.com/sneumann/xcms/issues/445 for more information.

- verboseBetaColumns:

  Option to calculate two additional metrics of peak quality via
  comparison to an idealized bell curve. Adds `beta_cor` and `beta_snr`
  to the `chromPeaks` output, corresponding to a Pearson correlation
  coefficient to a bell curve with several degrees of skew as well as an
  estimate of signal-to-noise using the residuals from the best-fitting
  bell curve. See https://github.com/sneumann/xcms/pull/685 and
  https://doi.org/10.1186/s12859-023-05533-4 for more information.

- snthreshIsoROIs:

  `numeric(1)` defining the signal to noise ratio cutoff to be used in
  the second centWave run to identify peaks for predicted isotope ROIs.

- maxCharge:

  `integer(1)` defining the maximal isotope charge. Isotopes will be
  defined for charges `1:maxCharge`.

- maxIso:

  `integer(1)` defining the number of isotope peaks that should be
  predicted for each peak identified in the first centWave run.

- mzIntervalExtension:

  `logical(1)` whether the mz range for the predicted isotope ROIs
  should be extended to increase detection of low intensity peaks.

- polarity:

  `character(1)` specifying the polarity of the data. Currently not
  used, but has to be `"positive"`, `"negative"` or `"unknown"` if
  provided.

- object:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  an
  [`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
  object containing the MS- and all other experiment-relevant data.

      For all other methods: a parameter object.

- param:

  An `CentWavePredIsoParam` object with the settings for the
  chromatographic peak detection algorithm.

- BPPARAM:

  A parameter class specifying if and how parallel processing should be
  performed. It defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  See documentation of the *BiocParallel* package for more details. If
  parallel processing is enabled, peak detection is performed in
  parallel on several of the input samples.

- return.type:

  Character specifying what type of object the method should return. Can
  be either `"XCMSnExp"` (default), `"list"` or `"xcmsSet"`.

- msLevel:

  `integer(1)` defining the MS level on which the peak detection should
  be performed. Defaults to `msLevel = 1`.

- ...:

  ignored.

## Value

The `CentWavePredIsoParam()` function returns a `CentWavePredIsoParam`
class instance with all of the settings specified for the two-step
centWave-based peak detection considering also isotopes.

For
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
if `return.type = "XCMSnExp"` an `XCMSnExp` object with the results of
the peak detection. If `return.type = "list"` a list of length equal to
the number of samples with matrices specifying the identified peaks. If
`return.type = "xcmsSet"` an `xcmsSet` object with the results of the
peak detection.

## Details

See
[`centWave()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md)
for details on the centWave method.

Parallel processing (one process per sample) is supported and can be
configured either by the `BPPARAM` parameter or by globally defining the
parallel processing mode using the
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
method from the *BiocParallel* package.

## See also

The
[`do_findChromPeaks_centWaveWithPredIsoROIs()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWaveWithPredIsoROIs.md)
core API function.

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for the object containing the results of the peak detection.

Other peak detection methods:
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md),
[`findChromPeaks-centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
[`findChromPeaks-massifquant`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
[`findChromPeaks-matchedFilter`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md),
[`findPeaks-MSW`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)

## Author

Hendrik Treutler, Johannes Rainer

## Examples

``` r

## Create a param object
p <- CentWavePredIsoParam(maxCharge = 4, snthresh = 25)
p
#> Object of class:  CentWavePredIsoParam 
#>  Parameters:
#>  - snthreshIsoROIs: [1] 6.25
#>  - maxCharge: [1] 4
#>  - maxIso: [1] 5
#>  - mzIntervalExtension: [1] TRUE
#>  - polarity: [1] "unknown"
#>  - ppm: [1] 25
#>  - peakwidth: [1] 20 50
#>  - snthresh: [1] 25
#>  - prefilter: [1]   3 100
#>  - mzCenterFun: [1] "wMean"
#>  - integrate: [1] 1
#>  - mzdiff: [1] -0.001
#>  - fitgauss: [1] FALSE
#>  - noise: [1] 0
#>  - verboseColumns: [1] FALSE
#>  - roiList: list()
#>  - firstBaselineCheck: [1] TRUE
#>  - roiScales: numeric(0)
#>  - extendLengthMSW: [1] FALSE
#>  - verboseBetaColumns: [1] FALSE
```
