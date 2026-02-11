# Median filtering of the profile matrix

Apply a median filter of given size to a profile matrix.

## Methods

- object = "xcmsRaw":

  `profMedFilt(object, massrad = 0, scanrad = 0)`

## Arguments

- object:

  the `xcmsRaw` object

- massrad:

  number of m/z grid points on either side to use for median calculation

- scanrad:

  number of scan grid points on either side to use for median
  calculation

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`medianFilter`](https://sneumann.github.io/xcms/reference/medianFilter.md)
