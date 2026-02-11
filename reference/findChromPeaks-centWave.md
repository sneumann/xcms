# Chromatographic peak detection using the centWave method

The centWave algorithm perform peak density and wavelet based
chromatographic peak detection for high resolution LC/MS data in
centroid mode *Tautenhahn 2008*.

The `findChromPeaks,OnDiskMSnExp,CentWaveParam()` method performs
chromatographic peak detection using the *centWave* algorithm on all
samples from an `OnDiskMSnExp` object. `OnDiskMSnExp` objects encapsule
all experiment specific data and load the spectra data (mz and intensity
values) on the fly from the original files applying also all eventual
data manipulations.

## Usage

``` r
CentWaveParam(
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
  verboseBetaColumns = FALSE
)

# S4 method for class 'OnDiskMSnExp,CentWaveParam'
findChromPeaks(
  object,
  param,
  BPPARAM = bpparam(),
  return.type = "XCMSnExp",
  msLevel = 1L,
  ...
)

# S4 method for class 'CentWaveParam'
as.list(x, ...)
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

- object:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  an
  [`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
  object containing the MS- and all other experiment-relevant data.

      For all other methods: a parameter object.

- param:

  An `CentWaveParam()` object containing all settings for the centWave
  algorithm.

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

- x:

  The parameter object.

## Value

The `CentWaveParam()` function returns a `CentWaveParam` class instance
with all of the settings specified for chromatographic peak detection by
the centWave method.

For
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
if `return.type = "XCMSnExp"` an
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object with the results of the peak detection. If `return.type = "list"`
a list of length equal to the number of samples with matrices specifying
the identified peaks. If `return.type = "xcmsSet"` an `xcmsSet` object
with the results of the peak detection.

## Details

The centWave algorithm is most suitable for high resolution
LC/{TOF,OrbiTrap,FTICR}-MS data in centroid mode. In the first phase the
method identifies *regions of interest* (ROIs) representing mass traces
that are characterized as regions with less than `ppm` m/z deviation in
consecutive scans in the LC/MS map. In detail, starting with a single
m/z, a ROI is extended if a m/z can be found in the next scan (spectrum)
for which the difference to the mean m/z of the ROI is smaller than the
user defined `ppm` of the m/z. The mean m/z of the ROI is then updated
considering also the newly included m/z value.

These ROIs are then, after some cleanup, analyzed using continuous
wavelet transform (CWT) to locate chromatographic peaks on different
scales. The first analysis step is skipped, if regions of interest are
passed *via* the `param` parameter.

Parallel processing (one process per sample) is supported and can be
configured either by the `BPPARAM` parameter or by globally defining the
parallel processing mode using the
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
method from the *BiocParallel* package.

## Note

These methods and classes are part of the updated and modernized *xcms*
user interface which will eventually replace the
[`findPeaks()`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
methods.

## References

Ralf Tautenhahn, Christoph Böttcher, and Steffen Neumann "Highly
sensitive feature detection for high resolution LC/MS" *BMC
Bioinformatics* 2008, 9:504 doi:
[10.1186/1471-2105-9-504](https://doi.org/10.1186/1471-2105-9-504)

## See also

The
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md)
core API function and
[`findPeaks.centWave()`](https://sneumann.github.io/xcms/reference/findPeaks.centWave-methods.md)
for the old user interface.

[`peaksWithCentWave()`](https://sneumann.github.io/xcms/reference/peaksWithCentWave.md)
for functions to perform centWave peak detection in purely
chromatographic data.

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for the object containing the results of the peak detection.

Other peak detection methods:
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md),
[`findChromPeaks-centWaveWithPredIsoROIs`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md),
[`findChromPeaks-massifquant`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
[`findChromPeaks-matchedFilter`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md),
[`findPeaks-MSW`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)

## Author

Ralf Tautenhahn, Johannes Rainer

## Examples

``` r

## Create a CentWaveParam object. Note that the noise is set to 10000 to
## speed up the execution of the example - in a real use case the default
## value should be used, or it should be set to a reasonable value.
cwp <- CentWaveParam(ppm = 25, noise = 10000, prefilter = c(3, 10000))
cwp
#> Object of class:  CentWaveParam 
#>  Parameters:
#>  - ppm: [1] 25
#>  - peakwidth: [1] 20 50
#>  - snthresh: [1] 10
#>  - prefilter: [1]     3 10000
#>  - mzCenterFun: [1] "wMean"
#>  - integrate: [1] 1
#>  - mzdiff: [1] -0.001
#>  - fitgauss: [1] FALSE
#>  - noise: [1] 10000
#>  - verboseColumns: [1] FALSE
#>  - roiList: list()
#>  - firstBaselineCheck: [1] TRUE
#>  - roiScales: numeric(0)
#>  - extendLengthMSW: [1] FALSE
#>  - verboseBetaColumns: [1] FALSE

## Perform the peak detection using centWave on some of the files from the
## faahKO package. Files are read using the `readMsExperiment` function
## from the MsExperiment package
library(faahKO)
library(xcms)
library(MsExperiment)
fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
           full.names = TRUE)
raw_data <- readMsExperiment(fls[1])

## Perform the peak detection using the settings defined above.
res <- findChromPeaks(raw_data, param = cwp)
head(chromPeaks(res))
#>          mz mzmin mzmax       rt    rtmin    rtmax      into      intb   maxo
#> CP001 453.2 453.2 453.2 2506.073 2501.378 2527.982 1007409.0 1007380.8  38152
#> CP002 307.0 307.0 307.0 2618.750 2592.145 2645.354  284782.4  268039.8  16872
#> CP003 302.0 302.0 302.0 2617.185 2595.275 2640.659  687146.6  671297.8  30552
#> CP004 360.0 360.0 360.0 2682.913 2668.828 2698.562 5641322.3 5420634.7 317568
#> CP005 361.1 361.1 361.1 2684.478 2665.698 2698.562 1158340.2 1116522.0  72272
#> CP006 416.1 416.1 416.1 2682.913 2635.964 2709.517  487698.6  446552.1  12036
#>          sn sample
#> CP001 38151      1
#> CP002    20      1
#> CP003    46      1
#> CP004    11      1
#> CP005    11      1
#> CP006    11      1
```
