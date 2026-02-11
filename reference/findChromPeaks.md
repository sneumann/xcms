# Chromatographic Peak Detection

The `findChromPeaks` method performs chromatographic peak detection on
LC/GC-MS data. The peak detection algorithm can be selected, and
configured, using the `param` argument.

Supported `param` objects are:

- [`CentWaveParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md):
  chromatographic peak detection using the *centWave* algorithm.

- [`CentWavePredIsoParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md):
  *centWave* with predicted isotopes. Peak detection uses a two-step
  centWave-based approach considering also feature isotopes.

- [`MatchedFilterParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md):
  peak detection using the *matched filter* algorithm.

- [`MassifquantParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md):
  peak detection using the Kalman filter-based *massifquant* method.

- [`MSWParam()`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md):
  single-spectrum non-chromatography MS data peak detection.

For specific examples see the help pages of the individual parameter
classes listed above.

## Usage

``` r
findChromPeaks(object, param, ...)

# S4 method for class 'MsExperiment,Param'
findChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  hdf5File = character(),
  force.overwrite = FALSE,
  ...,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperiment,Param'
findChromPeaks(
  object,
  param,
  msLevel = 1L,
  chunkSize = 2L,
  add = FALSE,
  ...,
  BPPARAM = bpparam()
)
```

## Arguments

- object:

  The data object on which to perform the peak detection. Can be an
  [`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html),
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md),
  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  or
  [`MsExperiment::MsExperiment()`](https://rdrr.io/pkg/MsExperiment/man/MsExperiment.html)
  object.

- param:

  The parameter object selecting and configuring the algorithm.

- ...:

  Optional parameters.

- msLevel:

  `integer(1)` defining the MS level on which the chromatographic peak
  detection should be performed.

- chunkSize:

  `integer(1)` for `object` being an `MsExperiment` or
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  defines the number of files (samples) for which the full peaks data
  (m/z and intensity values) should be loaded into memory at the same
  time. Peak detection is then performed in parallel (per sample) on
  this subset of loaded data. This setting thus allows to balance
  between memory demand and speed (due to parallel processing) of the
  peak detection. Because parallel processing can only performed on the
  subset of data loaded currently into memory (in each iteration), the
  value for `chunkSize` should be match the defined parallel setting
  setup. Using a parallel processing setup using 4 CPUs (separate
  processes) but using
  `chunkSize = `1`will not perform any parallel processing, as only the data from one sample is loaded in memory at a time. On the other hand, setting`chunkSize`to the total number of samples in an experiment will load the full MS data into memory and will thus in most settings cause an out-of-memory error. By setting`chunkSize
  =
  -1`the peak detection will be performed separately, and in parallel, for each sample. This will however not work for all`Spectra\`
  *backends* (see eventually
  [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  for details).

- hdf5File:

  For `object` being an `MsExperiment`: `character(1)` specifying the
  name (inclusive path) of a file that should be used for on-disk
  storage of preprocessing results. This option is suggested for very
  large data sets since it significantly reduces the memory demand. See
  [XcmsExperimentHdf5](https://sneumann.github.io/xcms/reference/XcmsExperimentHdf5.md)
  for more information. Note that an error is thrown if the file already
  exists. Overwriting an existing result file can be forced using
  `force.overwrite = TRUE`.

- force.overwrite:

  For `object` being an `MsExperiment` and parameter `hdf5File` being
  defined (see below): `logical(1)` whether an eventually existing
  result file should be overwritten.

- BPPARAM:

  Parallel processing setup. Uses by default the system-wide default
  setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more details.

- add:

  `logical(1)` (if `object` contains already chromatographic peaks, i.e.
  is either an `XCMSnExp` or `XcmsExperiment`) whether chromatographic
  peak detection results should be **added** to existing results. By
  default (`add = FALSE`) any additional `findChromPeaks` call on a
  result object will remove previous results.

## See also

[`plotChromPeaks()`](https://sneumann.github.io/xcms/reference/plotChromPeaks.md)
to plot identified chromatographic peaks for one file.

[`refineChromPeaks()`](https://sneumann.github.io/xcms/reference/refineChromPeaks.md)
for methods to *refine* or clean identified chromatographic peaks.

[`manualChromPeaks()`](https://sneumann.github.io/xcms/reference/manualChromPeaks.md)
to manually add/define chromatographic peaks.

Other peak detection methods:
[`findChromPeaks-centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
[`findChromPeaks-centWaveWithPredIsoROIs`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md),
[`findChromPeaks-massifquant`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
[`findChromPeaks-matchedFilter`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md),
[`findPeaks-MSW`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)

## Author

Johannes Rainer
