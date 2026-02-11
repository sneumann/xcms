# Get m/z and intensity values for a single mass scan

Return the data from a single mass scan using the numeric index of the
scan as a reference.

## Methods

- object = "xcmsRaw":

  `getScan(object, scan, mzrange = numeric())`
  `getMsnScan(object, scan, mzrange = numeric())`

## Arguments

- object:

  the `xcmsRaw` object

- scan:

  integer index of scan. if negative, the index numbered from the end

- mzrange:

  limit data points returned to those between in the range,
  `range(mzrange)`

## Value

A matrix with two columns:

- mz:

  m/z values

- intensity:

  intensity values

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`getSpec`](https://sneumann.github.io/xcms/reference/getSpec-methods.md)
