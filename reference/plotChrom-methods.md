# Plot extracted ion chromatograms from the profile matrix

Uses the pre-generated profile mode matrix to plot averaged or base peak
extracted ion chromatograms over a specified mass range.

## Methods

- object = "xcmsRaw":

  `plotChrom(object, base = FALSE, ident = FALSE, fitgauss = FALSE, vline = numeric(0), ...)`

## Arguments

- object:

  the `xcmsRaw` object

- base:

  logical, plot a base-peak chromatogram

- ident:

  logical, use mouse to identify and label peaks

- fitgauss:

  logical, fit a gaussian to the largest peak

- vline:

  numeric vector with locations of vertical lines

- ...:

  arguments passed to
  [`profRange`](https://sneumann.github.io/xcms/reference/profRange-methods.md)

## Value

If `ident == TRUE`, an integer vector with the indecies of the points
that were identified. If `fitgauss == TRUE`, a `nls` model with the
fitted gaussian. Otherwise a two-column matrix with the plotted points.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
