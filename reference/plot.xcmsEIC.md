# Plot extracted ion chromatograms from multiple files

Batch plot a list of extracted ion chromatograms to the current graphics
device.

## Methods

- x = "xcmsEIC":

  ` plot.xcmsEIC(x, y, groupidx = groupnames(x), sampleidx = sampnames(x), rtrange = x@rtrange, col = rep(1, length(sampleidx)), legtext = NULL, peakint = TRUE, sleep = 0, ...)`

## Value

A `xcmsSet` object.

## Arguments

- x:

  the `xcmsEIC` object

- y:

  optional `xcmsSet` object with peak integration data

- groupidx:

  either character vector with names or integer vector with indicies of
  peak groups for which to plot EICs

- sampleidx:

  either character vector with names or integer vector with indicies of
  samples for which to plot EICs

- rtrange:

  a two column matrix with minimum and maximum retention times between
  which to return EIC data points

  if it has the same number of rows as the number groups in the
  `xcmsEIC` object, then `sampleidx` is used to subset it. otherwise, it
  is repeated over the length of `sampleidx`

  it may also be a single number specifying the time window around the
  peak for which to plot EIC data

- col:

  color to use for plotting extracted ion chromatograms. if missing and
  `y` is specified, colors are taken from `unclass(sampclass(y))` and
  the default palette

  if it is the same length as the number groups in the `xcmsEIC` object,
  then `sampleidx` is used to subset it. otherwise, it is repeated over
  the length of `sampleidx`

- legtext:

  text to use for legend. if `NULL` and `y` is specified, legend text is
  taken from the sample class information found in the `xcmsSet`

- peakint:

  logical, plot integrated peak area with darkened lines (requires that
  `y` also be specified)

- sleep:

  seconds to pause between plotting EICs

- ...:

  other graphical parameters

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`xcmsEIC-class`](https://sneumann.github.io/xcms/reference/xcmsEIC-class.md),
[`png`](https://rdrr.io/r/grDevices/png.html),
[`pdf`](https://rdrr.io/r/grDevices/pdf.html),
[`postscript`](https://rdrr.io/r/grDevices/postscript.html),
