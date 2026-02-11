# Plot extracted ion chromatograms for specified m/z range

Plot extracted ion chromatogram for m/z values of interest. The raw data
is used in contrast to
[`plotChrom`](https://sneumann.github.io/xcms/reference/plotChrom-methods.md)
which uses data from the profile matrix.

## Methods

- object = "xcmsRaw":

  ` plotEIC(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric(), mzdec=2, type="l", add=FALSE, ...) `

## Arguments

- object:

  `xcmsRaw` object

- mzrange:

  m/z range for EIC. Uses the full m/z range by default.

- rtrange:

  retention time range for EIC. Uses the full retention time range by
  default.

- scanrange:

  scan range for EIC

- mzdec:

  Number of decimal places of title m/z values in the eic plot.

- type:

  Speficies how the data should be plotted (by default as a line).

- add:

  If the EIC should be added to an existing plot.

- ...:

  Additional parameters passed to the plotting function (e.g. `col`
  etc).

## Value

A two-column matrix with the plotted points.

## Author

Ralf Tautenhahn

## See also

[`rawEIC`](https://sneumann.github.io/xcms/reference/rawEIC-methods.md)`,`[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
