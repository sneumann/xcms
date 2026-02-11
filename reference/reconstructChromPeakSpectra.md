# Data independent acquisition (DIA): reconstruct MS2 spectra

*Reconstructs* MS2 spectra for each MS1 chromatographic peak (if
possible) for data independent acquisition (DIA) data (such as SWATH).
See the *LC-MS/MS analysis* vignette for more details and examples.

## Usage

``` r
reconstructChromPeakSpectra(object, ...)

# S4 method for class 'XcmsExperiment'
reconstructChromPeakSpectra(
  object,
  expandRt = 0,
  diffRt = 2,
  minCor = 0.8,
  intensity = "maxo",
  peakId = rownames(chromPeaks(object, msLevel = 1L)),
  BPPARAM = bpparam()
)

# S4 method for class 'XCMSnExp'
reconstructChromPeakSpectra(
  object,
  expandRt = 0,
  diffRt = 2,
  minCor = 0.8,
  intensity = "maxo",
  peakId = rownames(chromPeaks(object, msLevel = 1L)),
  BPPARAM = bpparam(),
  return.type = c("Spectra", "MSpectra")
)
```

## Arguments

- object:

  `XCMSnExp` with identified chromatographic peaks.

- ...:

  ignored.

- expandRt:

  `numeric(1)` allowing to expand the retention time range for extracted
  ion chromatograms by a constant value (for the peak shape
  correlation). Defaults to `expandRt = 0` hence correlates only the
  signal included in the identified chromatographic peaks.

- diffRt:

  `numeric(1)` defining the maximal allowed difference between the
  retention time of the chromatographic peak (apex) and the retention
  times of MS2 chromatographic peaks (apex) to consider them as
  representing candidate fragments of the original ion.

- minCor:

  `numeric(1)` defining the minimal required correlation coefficient for
  MS2 chromatographic peaks to be considered for MS2 spectrum
  reconstruction.

- intensity:

  `character(1)` defining the column in the `chromPeaks` matrix that
  should be used for the intensities of the reconstructed spectra's
  peaks. The same value from the MS1 chromatographic peaks will be used
  as `precursorIntensity` of the resulting spectra.

- peakId:

  optional `character` vector with peak IDs (i.e. rownames of
  `chromPeaks`) of MS1 peaks for which MS2 spectra should be
  reconstructed. By default they are reconstructed for all MS1
  chromatographic peaks.

- BPPARAM:

  parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more information.

- return.type:

  `character(1)` defining the type of the returned object. Only
  `return.type = "Spectra"` is supported, `return.type = "MSpectra"` is
  deprecated.

## Value

- [`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
  object (defined in the `Spectra` package) with the reconstructed MS2
  spectra for all MS1 peaks in `object`. Contains empty spectra (i.e.
  without m/z and intensity values) for MS1 peaks for which
  reconstruction was not possible (either no MS2 signal was recorded or
  the correlation of the MS2 chromatographic peaks with the MS1
  chromatographic peak was below threshold `minCor`. Spectra variables
  `"ms2_peak_id"` and `"ms2_peak_cor"` (of type
  [`IRanges::CharacterList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  and
  [`IRanges::NumericList()`](https://rdrr.io/pkg/IRanges/man/AtomicList-class.html)
  with length equal to the number of peaks per reconstructed MS2
  spectrum) providing the IDs and the correlation of the MS2
  chromatographic peaks from which the MS2 spectrum was reconstructed.
  As retention time the median retention times of all MS2
  chromatographic peaks used for the spectrum reconstruction is
  reported. The MS1 chromatographic peak intensity is reported as the
  reconstructed spectrum's `precursorIntensity` value (see parameter
  `intensity` above).

## Details

In detail, the function performs for each MS1 chromatographic peak:

- Identify all MS2 chromatographic peaks from the isolation window
  containing the m/z of the ion (i.e. the MS1 chromatographic peak) with
  approximately the same retention time than the MS1 peak (accepted rt
  shift can be specified with the `diffRt` parameter).

- Correlate the peak shapes of the candidate MS2 chromatographic peaks
  with the peak shape of the MS1 peak retaining only MS2 chromatographic
  peaks for which the correlation is `> minCor`.

- Reconstruct the MS2 spectrum using the m/z of all above selected MS2
  chromatographic peaks and their intensity (either `"maxo"` or
  `"into"`). Each MS2 chromatographic peak selected for an MS1 peak will
  thus represent one **mass peak** in the reconstructed spectrum.

The resulting
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object provides also the peak IDs of the MS2 chromatographic peaks for
each spectrum as well as their correlation value with spectra variables
*ms2_peak_id* and *ms2_peak_cor*.

## See also

[`findChromPeaksIsolationWindow()`](https://sneumann.github.io/xcms/reference/findChromPeaksIsolationWindow.md)
for the function to perform MS2 peak detection in DIA isolation windows
and for examples.

## Author

Johannes Rainer, Michael Witting
