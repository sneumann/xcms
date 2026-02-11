# Get a raw data matrix

Returns a matrix with columns for time, m/z, and intensity that
represents the raw data from a chromatography mass spectrometry
experiment.

## Methods

- object = "xcmsRaw":

  `rawMat(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric(), log=FALSE)`

## Arguments

- object:

  The container of the raw data

- mzrange:

  Subset by m/z range

- rtrange:

  Subset by retention time range

- scanrange:

  Subset by scan index range

- log:

  Whether to log transform the intensities

## Value

A numeric matrix with three columns: time, mz and intensity.

## See also

[`plotRaw`](https://sneumann.github.io/xcms/reference/plotRaw-methods.md)
for plotting the raw intensities

## Author

Michael Lawrence
