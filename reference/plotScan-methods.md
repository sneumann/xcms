# Plot a single mass scan

Plot a single mass scan using the impulse representation. Most useful
for centroided data.

## Methods

- object = "xcmsRaw":

  `plotScan(object, scan, mzrange = numeric(), ident = FALSE)`

## Arguments

- object:

  the `xcmsRaw` object

- scan:

  integer with number of scan to plot

- mzrange:

  numeric vector of length \>= 2 whose range will be used to select
  masses to plot

- ident:

  logical, use mouse to interactively identify and label individual
  masses

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
