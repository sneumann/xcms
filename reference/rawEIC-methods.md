# Get extracted ion chromatograms for specified m/z range

Generate extracted ion chromatogram for m/z values of interest. The raw
data is used in contrast to
[`getEIC`](https://sneumann.github.io/xcms/reference/getEIC-methods.md)
which uses data from the profile matrix (i.e. values binned along the
M/Z dimension).

## Methods

- object = "xcmsRaw":

  ` rawEIC(object, mzrange = numeric(), rtrange = numeric(), scanrange = numeric()) `

## Arguments

- object:

  `xcmsRaw` object

- mzrange:

  m/z range for EIC

- rtrange:

  retention time range for EIC

- scanrange:

  scan range for EIC

## Value

A list of :

- scan:

  scan number

- intensity:

  added intensity values

## Author

Ralf Tautenhahn

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
