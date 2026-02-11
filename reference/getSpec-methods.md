# Get average m/z and intensity values for multiple mass scans

Return full-resolution averaged data from multiple mass scans.

## Methods

- object = "xcmsRaw":

  `getSpec(object, ...)`

## Arguments

- object:

  the `xcmsRaw` object

- ...:

  arguments passed to
  [`profRange`](https://sneumann.github.io/xcms/reference/profRange-methods.md)
  used to sepecify the spectral segments of interest for averaging

## Details

Based on the mass points from the spectra selected, a master unique list
of masses is generated. Every spectra is interpolated at those masses
and then averaged.

## Value

A matrix with two columns:

- mz:

  m/z values

- intensity:

  intensity values

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`profRange`](https://sneumann.github.io/xcms/reference/profRange-methods.md),
[`getScan`](https://sneumann.github.io/xcms/reference/getScan-methods.md)
