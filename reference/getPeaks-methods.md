# Get peak intensities for specified regions

Integrate extracted ion chromatograms in pre-defined defined regions.
Return output similar to
[`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md).

## Methods

- object = "xcmsRaw":

  `getPeaks(object, peakrange, step = 0.1)`

## Arguments

- object:

  the `xcmsSet` object

- peakrange:

  matrix or data frame with 4 columns: `mzmin`, `mzmax`, `rtmin`,
  `rtmax` (they must be in that order or named)

- step:

  step size to use for profile generation

## Value

A matrix with columns:

- i:

  rank of peak identified in merged EIC (\<= `max`), always `NA`

- mz:

  weighted (by intensity) mean of peak m/z across scans

- mzmin:

  m/z of minimum step

- mzmax:

  m/z of maximum step

- ret:

  retention time of peak midpoint

- retmin:

  leading edge of peak retention time

- retmax:

  trailing edge of peak retention time

- into:

  integrated area of original (raw) peak

- intf:

  integrated area of filtered peak, always `NA`

- maxo:

  maximum intensity of original (raw) peak

- maxf:

  maximum intensity of filtered peak, always `NA`

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
