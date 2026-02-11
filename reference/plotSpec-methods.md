# Plot mass spectra from the profile matrix

Uses the pre-generated profile mode matrix to plot mass spectra over a
specified retention time range.

## Methods

- object = "xcmsRaw":

  `plotSpec(object, ident = FALSE, vline = numeric(0), ...)`

## Arguments

- object:

  the `xcmsRaw` object

- ident:

  logical, use mouse to identify and label peaks

- vline:

  numeric vector with locations of vertical lines

- ...:

  arguments passed to
  [`profRange`](https://sneumann.github.io/xcms/reference/profRange-methods.md)

## Value

If `ident == TRUE`, an integer vector with the indecies of the points
that were identified. Otherwise a two-column matrix with the plotted
points.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
