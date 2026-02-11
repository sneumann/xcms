# Scatterplot of raw data points

Produce a scatterplot showing raw data point location in retention time
and m/z. This plot is more useful for centroided data than continuum
data.

## Methods

- object = "xcmsRaw":

  `plotRaw(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric(), log=FALSE, title='Raw Data')`

## Arguments

- object:

  the `xcmsRaw` object

- mzrange:

  numeric vector of length \>= 2 whose range will be used to select the
  masses to plot

- rtrange:

  numeric vector of length \>= 2 whose range will be used to select the
  retention times to plot

- scanrange:

  numeric vector of length \>= 2 whose range will be used to select
  scans to plot

- log:

  logical, log transform intensity

- title:

  main title of the plot

## Value

A matrix with the points plotted.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
