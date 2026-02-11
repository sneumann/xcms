# Core API function for massifquant peak detection

Massifquant is a Kalman filter (KF)-based chromatographic peak detection
for XC-MS data in centroid mode. The identified peaks can be further
refined with the *centWave* method (see
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md)
for details on centWave) by specifying `withWave = TRUE`.

## Usage

``` r
do_findChromPeaks_massifquant(
  mz,
  int,
  scantime,
  valsPerSpect,
  ppm = 10,
  peakwidth = c(20, 50),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1,
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
```

## Arguments

- mz:

  Numeric vector with the individual m/z values from all scans/ spectra
  of one file/sample.

- int:

  Numeric vector with the individual intensity values from all
  scans/spectra of one file/sample.

- scantime:

  Numeric vector of length equal to the number of spectra/scans of the
  data representing the retention time of each scan.

- valsPerSpect:

  Numeric vector with the number of values for each spectrum.

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

## Value

A matrix, each row representing an identified chromatographic peak, with
columns:

- `"mz"`: Intensity weighted mean of m/z values of the peaks across
  scans.

- `"mzmin"`: Minumum m/z of the peak.

- `"mzmax"`: Maximum m/z of the peak.

- `"rtmin"`: Minimum retention time of the peak.

- `"rtmax"`: Maximum retention time of the peak.

- `"rt"`: Retention time of the peak's midpoint.

- `"into"`: Integrated (original) intensity of the peak.

- `"maxo"`: Maximum intensity of the peak.

If `withWave` is set to `TRUE`, the result is the same as returned by
the
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md)
method.

## Details

This algorithm's performance has been tested rigorously on high
resolution LC/(OrbiTrap, TOF)-MS data in centroid mode. Simultaneous
kalman filters identify peaks and calculate their area under the curve.
The default parameters are set to operate on a complex LC-MS Orbitrap
sample. Users will find it useful to do some simple exploratory data
analysis to find out where to set a minimum intensity, and identify how
many scans an average peak spans. The `consecMissedLimit` parameter has
yielded good performance on Orbitrap data when set to (`2`) and on TOF
data it was found best to be at (`1`). This may change as the algorithm
has yet to be tested on many samples. The `criticalValue` parameter is
perhaps most dificult to dial in appropriately and visual inspection of
peak identification is the best suggested tool for quick optimization.
The `ppm` and `checkBack` parameters have shown less influence than the
other parameters and exist to give users flexibility and better
accuracy.

## References

Conley CJ, Smith R, Torgrip RJ, Taylor RM, Tautenhahn R and Prince JT
"Massifquant: open-source Kalman filter-based XC-MS isotope trace
feature detection" *Bioinformatics* 2014, 30(18):2636-43. doi:
[10.1093/bioinformatics/btu359](https://doi.org/10.1093/bioinformatics/btu359)

## See also

Other core peak detection functions:
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md),
[`do_findChromPeaks_centWaveWithPredIsoROIs()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWaveWithPredIsoROIs.md),
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md),
[`do_findPeaks_MSW()`](https://sneumann.github.io/xcms/reference/do_findPeaks_MSW.md)

## Author

Christopher Conley

## Examples

``` r

## Load the test file
faahko_sub <- loadXcmsData("faahko_sub")

## Subset to one file and restrict to a certain retention time range
data <- filterRt(filterFile(faahko_sub, 1), c(2500, 3000))

## Get m/z and intensity values
mzs <- mz(data)
ints <- intensity(data)

## Define the values per spectrum:
valsPerSpect <- lengths(mzs)

## Perform the peak detection using massifquant - setting prefilter to
## a high value to speed up the call for the example
res <- do_findChromPeaks_massifquant(mz = unlist(mzs), int = unlist(ints),
    scantime = rtime(data), valsPerSpect = valsPerSpect,
    prefilter = c(3, 10000))
#> 
#>  Massifquant, Copyright (C) 2013 Brigham Young University.
#>  Massifquant comes with ABSOLUTELY NO WARRANTY. See LICENSE for details.
#> 
#>  Detecting  mass traces at 10ppm ... 
#> OK
#>  69 Peaks.
head(res)
#>         mz mzmin mzmax    rtmin    rtmax       rt    into  maxo
#> [1,] 426.1 426.1 426.1 2963.039 2999.033 2977.124  129958 29408
#> [2,] 590.3 590.3 590.3 2963.039 2999.033 2978.689  255317 19072
#> [3,] 454.1 454.1 454.1 2934.870 2975.559 2953.650  372685 20432
#> [4,] 309.1 309.1 309.1 2905.136 2967.734 2936.435  652140 38888
#> [5,] 475.2 475.2 475.2 2880.097 2999.033 2934.870 2415224 37160
#> [6,] 532.2 532.2 532.2 2853.493 2906.701 2878.532  232493 16480
```
