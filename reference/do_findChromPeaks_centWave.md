# Core API function for centWave peak detection

This function performs peak density and wavelet based chromatographic
peak detection for high resolution LC/MS data in centroid mode
*Tautenhahn 2008*.

## Usage

``` r
do_findChromPeaks_centWave(
  mz,
  int,
  scantime,
  valsPerSpect,
  ppm = 25,
  peakwidth = c(20, 50),
  snthresh = 10,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1,
  mzdiff = -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  roiList = list(),
  firstBaselineCheck = TRUE,
  roiScales = NULL,
  sleep = 0,
  extendLengthMSW = FALSE,
  verboseBetaColumns = FALSE
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

- sleep:

  `numeric(1)` defining the number of seconds to wait between
  iterations. Defaults to `sleep = 0`. If `> 0` a plot is generated
  visualizing the identified chromatographic peak. Note: this argument
  is for backward compatibility only and will be removed in future.

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

## Value

A matrix, each row representing an identified chromatographic peak, with
columns:

- `"mz"`: Intensity weighted mean of m/z values of the peak across
  scans.

- `"mzmin"`: Minimum m/z of the peak.

- `"mzmax"`: Maximum m/z of the peak.

- `"rt"`: Retention time of the peak's midpoint.

- `"rtmin"`: Minimum retention time of the peak.

- \`"rtmax: Maximum retention time of the peak.

- `"into"`: Integrated (original) intensity of the peak.

- `"intb"`: Per-peak baseline corrected integrated peak intensity.

- `"maxo"`: Maximum intensity of the peak.

- `"sn"`: Signal to noise ratio, defined as `(maxo - baseline)/sd`, `sd`
  being the standard deviation of local chromatographic noise.

- `"egauss"`: RMSE of Gaussian fit.

Additional columns for `verboseColumns = TRUE`:

- `"mu"`: Gaussian parameter mu.

- `"sigma"`: Gaussian parameter sigma.

- `"h"`: Gaussian parameter h.

- `"f"`: Region number of the m/z ROI where the peak was localized.

- `"dppm"`: m/z deviation of mass trace across scans in ppm.

- `"scale"`: Scale on which the peak was localized.

- `"scpos"`: Peak position found by wavelet analysis (scan number).

- `"scmin"`: Left peak limit found by wavelet analysis (scan number).

- `"scmax"`: Right peak limit found by wavelet analysis (scan numer).

Additional columns for `verboseBetaColumns = TRUE`:

- `"beta_cor"`: Correlation between an "ideal" bell curve and the raw
  data.

- `"beta_snr"`: Signal-to-noise residuals calculated from the beta_cor
  fit.

## Details

This algorithm is most suitable for high resolution
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
passed with the `roiList` parameter.

## Note

The *centWave* was designed to work on centroided mode, thus it is
expected that such data is presented to the function.

    This function exposes core chromatographic peak detection functionality
    of the *centWave* method. While this function can be called
    directly, users will generally call the corresponding method for the
    data object instead.

## References

Ralf Tautenhahn, Christoph Böttcher, and Steffen Neumann "Highly
sensitive feature detection for high resolution LC/MS" *BMC
Bioinformatics* 2008, 9:504 doi:
[10.1186/1471-2105-9-504](https://doi.org/10.1186/1471-2105-9-504)

## See also

Other core peak detection functions:
[`do_findChromPeaks_centWaveWithPredIsoROIs()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWaveWithPredIsoROIs.md),
[`do_findChromPeaks_massifquant()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_massifquant.md),
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md),
[`do_findPeaks_MSW()`](https://sneumann.github.io/xcms/reference/do_findPeaks_MSW.md)

## Author

Ralf Tautenhahn, Johannes Rainer

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

## Calling the function. We're using a large value for noise and prefilter
## to speed up the call in the example - in a real use case we would either
## set the value to a reasonable value or use the default value.
res <- do_findChromPeaks_centWave(mz = unlist(mzs), int = unlist(ints),
    scantime = rtime(data), valsPerSpect = valsPerSpect, noise = 10000,
    prefilter = c(3, 10000))
#> Detecting mass traces at 25 ppm ... 
#> OK
#> Detecting chromatographic peaks in 186 regions of interest ...
#>  OK: 47 found.
head(res)
#>         mz mzmin mzmax       rt    rtmin    rtmax      into      intb   maxo
#> [1,] 453.2 453.2 453.2 2506.073 2501.378 2527.982 1007409.0 1007380.8  38152
#> [2,] 307.0 307.0 307.0 2618.750 2592.145 2645.354  284782.4  268039.8  16872
#> [3,] 302.0 302.0 302.0 2617.185 2595.275 2640.659  687146.6  671297.8  30552
#> [4,] 360.0 360.0 360.0 2682.913 2668.828 2698.562 5641322.3 5420634.7 317568
#> [5,] 361.1 361.1 361.1 2684.478 2665.698 2698.562 1158340.2 1116522.0  72272
#> [6,] 416.1 416.1 416.1 2682.913 2635.964 2709.517  487698.6  446552.1  12036
#>         sn
#> [1,] 38151
#> [2,]    20
#> [3,]    46
#> [4,]    11
#> [5,]    11
#> [6,]    11
```
