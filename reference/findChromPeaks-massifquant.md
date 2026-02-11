# Chromatographic peak detection using the massifquant method

Massifquant is a Kalman filter (KF)-based chromatographic peak detection
for XC-MS data in centroid mode. The identified peaks can be further
refined with the *centWave* method (see
[`findChromPeaks-centWave()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md)
for details on centWave) by specifying `withWave = TRUE`.

The `findChromPeaks,OnDiskMSnExp,MassifquantParam()` method performs
chromatographic peak detection using the *massifquant* algorithm on all
samples from an `OnDiskMSnExp` object. `OnDiskMSnExp` objects encapsule
all experiment specific data and load the spectra data (mz and intensity
values) on the fly from the original files applying also all eventual
data manipulations.

## Usage

``` r
MassifquantParam(
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
  criticalValue = 1.125,
  consecMissedLimit = 2,
  unions = 1,
  checkBack = 0,
  withWave = FALSE
)

# S4 method for class 'OnDiskMSnExp,MassifquantParam'
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

  `numeric(2)`. Only the first element is used by massifquant, which
  specifices the minimum peak length in time scans. For
  `withWave = TRUE` the second argument represents the maximum peak
  length subject to being greater than the mininum peak length (see also
  documentation of
  [`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md)).

- snthresh:

  `numeric(1)` defining the signal to noise ratio cutoff.

- prefilter:

  `numeric(2)`. The first argument is only used if (`withWave = TRUE`);
  see
  [`findChromPeaks-centWave()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md)
  for details. The second argument specifies the minimum threshold for
  the maximum intensity of a chromatographic peak that must be met.

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

- criticalValue:

  `numeric(1)`. Suggested values: (`0.1-3.0`). This setting helps
  determine the the Kalman Filter prediciton margin of error. A real
  centroid belonging to a bonafide peak must fall within the KF
  prediction margin of error. Much like in the construction of a
  confidence interval, `criticalVal` loosely translates to be a
  multiplier of the standard error of the prediction reported by the
  Kalman Filter. If the peak in the XC-MS sample have a small mass
  deviance in ppm error, a smaller critical value might be better and
  vice versa.

- consecMissedLimit:

  `integer(1)` Suggested values: (`1,2,3`). While a peak is in the
  proces of being detected by a Kalman Filter, the Kalman Filter may not
  find a predicted centroid in every scan. After 1 or more consecutive
  failed predictions, this setting informs Massifquant when to stop a
  Kalman Filter from following a candidate peak.

- unions:

  `integer(1)` set to `1` if apply t-test union on segmentation; set to
  `0` if no t-test to be applied on chromatographically continous peaks
  sharing same m/z range. Explanation: With very few data points,
  sometimes a Kalman Filter stops tracking a peak prematurely. Another
  Kalman Filter is instantiated and begins following the rest of the
  signal. Because tracking is done backwards to forwards, this
  algorithmic defect leaves a real peak divided into two segments or
  more. With this option turned on, the program identifies segmented
  peaks and combines them (merges them) into one with a two sample
  t-test. The potential danger of this option is that some truly
  distinct peaks may be merged.

- checkBack:

  `integer(1)` set to `1` if turned on; set to `0` if turned off. The
  convergence of a Kalman Filter to a peak's precise m/z mapping is very
  fast, but sometimes it incorporates erroneous centroids as part of a
  peak (especially early on). The `scanBack` option is an attempt to
  remove the occasional outlier that lies beyond the converged bounds of
  the Kalman Filter. The option does not directly affect identification
  of a peak because it is a postprocessing measure; it has not shown to
  be a extremely useful thus far and the default is set to being turned
  off.

- withWave:

  `logical(1)` if `TRUE`, the peaks identified first with Massifquant
  are subsequently filtered with the second step of the centWave
  algorithm, which includes wavelet estimation.

- object:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  an `OnDiskMSnExp` object containing the MS- and all other
  experiment-relevant data.

      For all other methods: a parameter object.

- param:

  An `MassifquantParam` object containing all settings for the
  massifquant algorithm.

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

The `MassifquantParam()` function returns a `MassifquantParam` class
instance with all of the settings specified for chromatographic peak
detection by the *massifquant* method.

For
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
if `return.type = "XCMSnExp"` an `XCMSnExp` object with the results of
the peak detection. If `return.type = "list"` a list of length equal to
the number of samples with matrices specifying the identified peaks. If
`return.type = "xcmsSet"` an `xcmsSet` object with the results of the
peak detection.

## Details

This algorithm's performance has been tested rigorously on high
resolution LC/(OrbiTrap, TOF)-MS data in centroid mode. Simultaneous
kalman filters identify chromatographic peaks and calculate their area
under the curve. The default parameters are set to operate on a complex
LC-MS Orbitrap sample. Users will find it useful to do some simple
exploratory data analysis to find out where to set a minimum intensity,
and identify how many scans an average peak spans. The
`consecMissedLimit` parameter has yielded good performance on Orbitrap
data when set to (`2`) and on TOF data it was found best to be at (`1`).
This may change as the algorithm has yet to be tested on many samples.
The `criticalValue` parameter is perhaps most dificult to dial in
appropriately and visual inspection of peak identification is the best
suggested tool for quick optimization. The `ppm` and `checkBack`
parameters have shown less influence than the other parameters and exist
to give users flexibility and better accuracy.

Parallel processing (one process per sample) is supported and can be
configured either by the `BPPARAM` parameter or by globally defining the
parallel processing mode using the
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
method from the *BiocParallel* package.

## References

Conley CJ, Smith R, Torgrip RJ, Taylor RM, Tautenhahn R and Prince JT
"Massifquant: open-source Kalman filter-based XC-MS isotope trace
feature detection" *Bioinformatics* 2014, 30(18):2636-43. doi:
[10.1093/bioinformatics/btu359](https://doi.org/10.1093/bioinformatics/btu359)

## See also

The
[`do_findChromPeaks_massifquant()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_massifquant.md)
core API function and
[`findPeaks.massifquant()`](https://sneumann.github.io/xcms/reference/findPeaks.massifquant-methods.md)
for the old user interface.

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for the object containing the results of the peak detection.

Other peak detection methods:
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md),
[`findChromPeaks-centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
[`findChromPeaks-centWaveWithPredIsoROIs`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md),
[`findChromPeaks-matchedFilter`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md),
[`findPeaks-MSW`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)

## Author

Christopher Conley, Johannes Rainer

## Examples

``` r

## Create a MassifquantParam object.
mqp <- MassifquantParam(snthresh = 30, prefilter = c(6, 10000))
mqp
#> Object of class:  MassifquantParam 
#>  Parameters:
#>  - ppm: [1] 25
#>  - peakwidth: [1] 20 50
#>  - snthresh: [1] 30
#>  - prefilter: [1]     6 10000
#>  - mzCenterFun: [1] "wMean"
#>  - integrate: [1] 1
#>  - mzdiff: [1] -0.001
#>  - fitgauss: [1] FALSE
#>  - noise: [1] 0
#>  - verboseColumns: [1] FALSE
#>  - criticalValue: [1] 1.125
#>  - consecMissedLimit: [1] 2
#>  - unions: [1] 1
#>  - checkBack: [1] 0
#>  - withWave: [1] FALSE

## Perform the peak detection using massifquant on the files from the
## faahKO package. Files are read using the readMSData from the MSnbase
## package
library(faahKO)
library(MSnbase)
fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
           full.names = TRUE)
raw_data <- readMSData(fls[1], mode = "onDisk")
#> Polarity can not be extracted from netCDF files, please set manually the polarity with the 'polarity' method.
## Perform the peak detection using the settings defined above.
res <- findChromPeaks(raw_data, param = mqp)
#> 
#>  Massifquant, Copyright (C) 2013 Brigham Young University.
#>  Massifquant comes with ABSOLUTELY NO WARRANTY. See LICENSE for details.
#> 
#>  Detecting  mass traces at 25ppm ... 
#> OK
#>  334 Peaks.
head(chromPeaks(res))
#>          mz mzmin mzmax    rtmin    rtmax       rt    into   maxo sample
#> CP001 578.4 578.4 578.4 4103.890 4141.449 4119.540 1748295 187328      1
#> CP002 577.3 577.3 577.3 4100.760 4130.495 4114.845 3854210 457024      1
#> CP003 367.2 367.2 367.2 4099.196 4164.924 4122.670  140753  11330      1
#> CP004 459.1 459.1 459.1 4099.196 4128.930 4113.280  144739  15874      1
#> CP005 401.2 401.2 401.2 4083.546 4221.262 4147.709  472407  11783      1
#> CP006 412.3 412.3 412.3 4072.591 4111.715 4089.806 2192350 183936      1
```
