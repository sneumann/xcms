# Data independent acquisition (DIA): peak detection in isolation windows

The `findChromPeaksIsolationWindow` function allows to perform a
chromatographic peak detection in MS level \> 1 spectra of certain
isolation windows (e.g. SWATH pockets). The function performs a peak
detection, separately for all spectra belonging to the same isolation
window and adds them to the
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
matrix of the result object. Information about the isolation window in
which they were detected is added to
[`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
data frame.

Note that peak detection with this method does not remove previously
identified chromatographic peaks (e.g. on MS1 level using the
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
function but adds newly identified peaks to the existing
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
matrix.

Isolation windows can be defined with the `isolationWindow` parameter,
that by default uses the definition of
[`isolationWindowTargetMz()`](https://sneumann.github.io/xcms/reference/isolationWindowTargetMz-OnDiskMSnExp-method.md),
i.e. chromatographic peak detection is performed for all spectra with
the same isolation window target m/z (seprarately for each file). The
parameter `param` allows to define and configure the peak detection
algorithm (see
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
for more information).

## Usage

``` r
findChromPeaksIsolationWindow(object, ...)

# S4 method for class 'MsExperiment'
findChromPeaksIsolationWindow(
  object,
  param,
  msLevel = 2L,
  isolationWindow = isolationWindowTargetMz(spectra(object)),
  chunkSize = 2L,
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'OnDiskMSnExp'
findChromPeaksIsolationWindow(
  object,
  param,
  msLevel = 2L,
  isolationWindow = isolationWindowTargetMz(object),
  ...
)
```

## Arguments

- object:

  `MsExperiment`, `XcmsExperiment`, `OnDiskMSnExp` or `XCMSnExp` object
  with the DIA data.

- ...:

  currently not used.

- param:

  Peak detection parameter object, such as a
  [CentWaveParam](https://sneumann.github.io/xcms/reference/hidden_aliases.md)
  object defining and configuring the chromographic peak detection
  algorithm. See also
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
  for more details.

- msLevel:

  `integer(1)` specifying the MS level in which the peak detection
  should be performed. By default `msLevel = 2L`.

- isolationWindow:

  `factor` or similar defining the isolation windows in which the peak
  detection should be performed with length equal to the number of
  spectra in `object`.

- chunkSize:

  if `object` is an `MsExperiment` or `XcmsExperiment`: `integer(1)`
  defining the number of files (samples) that should be loaded into
  memory and processed at a time. See
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
  for more information.

- BPPARAM:

  if `object` is an `MsExperiment` or `XcmsExperiment`: parallel
  processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

## Value

An `XcmsExperiment` or `XCMSnExp` object with the chromatographic peaks
identified in spectra of each isolation window from each file added to
the `chromPeaks` matrix. Isolation window definition for each identified
peak are stored as additional columns in
[`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md).

## See also

[`reconstructChromPeakSpectra()`](https://sneumann.github.io/xcms/reference/reconstructChromPeakSpectra.md)
for the function to reconstruct MS2 spectra for each MS1 chromatographic
peak.

## Author

Johannes Rainer, Michael Witting
