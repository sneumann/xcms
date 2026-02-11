# Peak detection in the chromatographic time domain

Find peaks in the chromatographic time domain of the profile matrix. For
more details see
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md).

## Usage

``` r
# S4 method for class 'xcmsRaw'
findPeaks.matchedFilter(
  object,
  fwhm = 30,
  sigma = fwhm/2.3548,
  max = 5,
  snthresh = 10,
  step = 0.1,
  steps = 2,
  mzdiff = 0.8 - step * steps,
  index = FALSE,
  sleep = 0,
  scanrange = numeric()
)
```

## Arguments

- object:

  The `xcmsRaw` object on which peak detection should be performed.

- fwhm:

  `numeric(1)` specifying the full width at half maximum of matched
  filtration gaussian model peak. Only used to calculate the actual
  sigma, see below.

- sigma:

  `numeric(1)` specifying the standard deviation (width) of the matched
  filtration model peak.

- max:

  `numeric(1)` representing the maximum number of peaks that are
  expected/will be identified per slice.

- snthresh:

  `numeric(1)` defining the signal to noise cutoff to be used in the
  chromatographic peak detection step.

- step:

  numeric(1) specifying the width of the bins/slices in m/z dimension.

- steps:

  `numeric(1)` defining the number of bins to be merged before
  filtration (i.e. the number of neighboring bins that will be joined to
  the slice in which filtration and peak detection will be performed).

- mzdiff:

  `numeric(1)` defining the minimum difference in m/z for peaks with
  overlapping retention times

- index:

  `logical(1)` specifying whether indicies should be returned instead of
  values for m/z and retention times.

- sleep:

  (DEPRECATED). The use of this parameter is highly discouraged, as it
  could cause problems in parallel processing mode.

- scanrange:

  Numeric vector defining the range of scans to which the original
  `object` should be sub-setted before peak detection.

## Value

A matrix, each row representing an intentified chromatographic peak.

## References

Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
*Anal. Chem.* 2006, 78:779-787. doi:
[10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

## Author

Colin A. Smith
