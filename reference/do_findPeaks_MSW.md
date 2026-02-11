# Core API function for single-spectrum non-chromatography MS data peak detection

This function performs peak detection in mass spectrometry direct
injection spectrum using a wavelet based algorithm.

## Usage

``` r
do_findPeaks_MSW(
  mz,
  int,
  snthresh = 3,
  verboseColumns = FALSE,
  scantime = numeric(),
  valsPerSpect = integer(),
  ...
)
```

## Arguments

- mz:

  Numeric vector with the individual m/z values from all scans/ spectra
  of one file/sample.

- int:

  Numeric vector with the individual intensity values from all
  scans/spectra of one file/sample.

- snthresh:

  `numeric(1)` defining the signal to noise ratio cutoff.

- verboseColumns:

  `logical(1)` whether additional peak meta data columns should be
  returned.

- scantime:

  ignored.

- valsPerSpect:

  ignored.

- ...:

  Additional parameters to be passed to the `peakDetectionCWT` function.

## Value

A matrix, each row representing an identified peak, with columns:

- `"mz"`: m/z value of the peak at the centroid position.

- `"mzmin"`: Minimum m/z of the peak.

- `"mzmax"`: Maximum m/z of the peak.

- `"rt"`: Always `-1`.

- `"rtmin"`: Always `-1`.

- `"rtmax"`: Always `-1`.

- `"into"`: Integrated (original) intensity of the peak.

- `"maxo"`: Maximum intensity of the peak.

- `"intf"`: Always `NA`.

- `"maxf"`: Maximum MSW-filter response of the peak.

- `"sn"`: Signal to noise ratio.

## Details

This is a wrapper around the peak picker in Bioconductor's
*MassSpecWavelet* package calling `peakDetectionCWT()` and
`tuneInPeakInfo()` functions. See the *xcmsDirect* vignette for more
information.

## See also

Other core peak detection functions:
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md),
[`do_findChromPeaks_centWaveWithPredIsoROIs()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWaveWithPredIsoROIs.md),
[`do_findChromPeaks_massifquant()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_massifquant.md),
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md)

## Author

Joachim Kutzera, Steffen Neumann, Johannes Rainer
