# Peak detection for single-spectrum non-chromatography MS data

This method performs peak detection in mass spectrometry direct
injection spectrum using a wavelet based algorithm.

## Usage

``` r
# S4 method for class 'xcmsRaw'
findPeaks.MSW(object, snthresh = 3, verbose.columns = FALSE, ...)
```

## Arguments

- object:

  The `xcmsRaw` object on which peak detection should be performed.

- snthresh:

  `numeric(1)` defining the signal to noise ratio cutoff.

- verbose.columns:

  Logical whether additional peak meta data columns should be returned.

- ...:

  Additional parameters to be passed to the `peakDetectionCWT()` and
  `identifyMajorPeaks()` functions from the *MassSpecWavelet* package.

## Value

    A matrix, each row representing an intentified peak.

## Details

This is a wrapper around the peak picker in Bioconductor's
*MassSpecWavelet* package calling `peakDetectionCWT` and
`tuneInPeakInfo` functions.

## Author

Joachim Kutzera, Steffen Neumann, Johannes Rainer
