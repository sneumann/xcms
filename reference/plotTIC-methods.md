# Plot total ion count

Plot chromatogram of total ion count. Optionally allow identification of
target peaks and viewing/identification of individual spectra.

## Methods

- object = "xcmsRaw":

  `plotTIC(object, ident = FALSE, msident = FALSE)`

## Arguments

- object:

  the `xcmsRaw` object

- ident:

  logical, use mouse to identify and label chromatographic peaks

- msident:

  logical, use mouse to identify and label spectral peaks

## Value

If `ident == TRUE`, an integer vector with the indecies of the points
that were identified. Otherwise a two-column matrix with the plotted
points.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
