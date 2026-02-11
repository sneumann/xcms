# Plot a grid of a large number of peaks

Plot extracted ion chromatograms for many peaks simultaneously,
indicating peak integration start and end points with vertical grey
lines.

## Methods

- `signature(object = "xcmsSet")`:

  `plotPeaks(object, peaks, figs, width = 200)`

## Arguments

- object:

  the `xcmsRaw` object

- peaks:

  matrix with peak information as produced by
  [`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)

- figs:

  two-element vector describing the number of rows and the number of
  columns of peaks to plot, if missing then an approximately square grid
  that will fit the number of peaks supplied

- width:

  width of chromatogram retention time to plot for each peak

## Details

This function is intended to help graphically analyze the results of
peak picking. It can help estimate the number of false positives and
improper integration start and end points. Its output is very compact
and tries to waste as little space as possible. Each plot is labeled
with rounded m/z and retention time separated by a space.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md),
[`split.screen`](https://rdrr.io/r/graphics/screen.html)
