# Core API function for two-step centWave peak detection with isotopes

The `do_findChromPeaks_centWaveWithPredIsoROIs` performs a two-step
centWave based peak detection: chromatographic peaks are identified
using centWave followed by a prediction of the location of the
identified peaks' isotopes in the mz-retention time space. These
locations are fed as *regions of interest* (ROIs) to a subsequent
centWave run. All non overlapping peaks from these two peak detection
runs are reported as the final list of identified peaks.

The `do_findChromPeaks_centWaveAddPredIsoROIs` performs centWave based
peak detection based in regions of interest (ROIs) representing
predicted isotopes for the peaks submitted with argument `peaks`. The
function returns a matrix with the identified peaks consisting of all
input peaks and peaks representing predicted isotopes of these (if found
by the centWave algorithm).

## Usage

``` r
do_findChromPeaks_centWaveWithPredIsoROIs(
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
  snthreshIsoROIs = 6.25,
  maxCharge = 3,
  maxIso = 5,
  mzIntervalExtension = TRUE,
  polarity = "unknown",
  extendLengthMSW = FALSE,
  verboseBetaColumns = FALSE
)

do_findChromPeaks_addPredIsoROIs(
  mz,
  int,
  scantime,
  valsPerSpect,
  ppm = 25,
  peakwidth = c(20, 50),
  snthresh = 6.25,
  prefilter = c(3, 100),
  mzCenterFun = "wMean",
  integrate = 1,
  mzdiff = -0.001,
  fitgauss = FALSE,
  noise = 0,
  verboseColumns = FALSE,
  peaks. = NULL,
  maxCharge = 3,
  maxIso = 5,
  mzIntervalExtension = TRUE,
  polarity = "unknown"
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

  For `do_findChromPeaks_addPredIsoROIs`: `numeric(1)` defining the
  signal to noise threshold for the centWave algorithm. For
  `do_findChromPeaks_centWaveWithPredIsoROIs`: `numeric(1)` defining the
  signal to noise threshold for the initial (first) centWave run.

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

- peaks.:

  A matrix such as one returned by a call to
  [`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md)
  (with `verboseColumns = TRUE`) with the peaks for which isotopes
  should be predicted and used for an additional peak detectoin using
  the centWave method. Required columns are: `"mz"`, `"mzmin"`,
  `"mzmax"`, `"scmin"`, `"scmax"`, `"scale"` and `"into"`.

## Value

A matrix, each row representing an identified chromatographic peak. All
non-overlapping peaks identified in both centWave runs are reported. The
matrix columns are:

- `"mz"`: Intensity weighted mean of m/z values of the peaks across
  scans.

- `"mzmin"`: Minimum m/z of the peaks.

- `"mzmax"`: Maximum m/z of the peaks.

- `"rt"`: Retention time of the peak's midpoint.

- `"rtmin"`: Minimum retention time of the peak.

- `"rtmax"`: Maximum retention time of the peak.

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

For more details on the centWave algorithm see
[`centWave()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md).

## See also

Other core peak detection functions:
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md),
[`do_findChromPeaks_massifquant()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_massifquant.md),
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md),
[`do_findPeaks_MSW()`](https://sneumann.github.io/xcms/reference/do_findPeaks_MSW.md)

## Author

Hendrik Treutler, Johannes Rainer
