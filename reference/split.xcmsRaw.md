# Divide an xcmsRaw object

Divides the scans from a `xcmsRaw` object into a list of multiple
objects. MS\$^n\$ data is discarded.

## Methods

- xr = "xcmsRaw":

  ` split(x, f, drop = TRUE, ...) `

## Arguments

- x:

  `xcmsRaw` object

- f:

  factor such that `factor(f)` defines the scans which go into the new
  `xcmsRaw` objects

- drop:

  logical indicating if levels that do not occur should be dropped (if
  'f' is a 'factor' or a list).

- ...:

  further potential arguments passed to methods.

## Value

A list of `xcmsRaw` objects.

## Author

Steffen Neumann, <sneumann@ipb-halle.de>

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
