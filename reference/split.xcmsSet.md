# Divide an xcmsSet object

Divides the samples and peaks from a `xcmsSet` object into a list of
multiple objects. Group data is discarded.

## Methods

- xs = "xcmsSet":

  ` split(x, f, drop = TRUE, ...) `

## Arguments

- xs:

  `xcmsSet` object

- f:

  factor such that `factor(f)` defines the grouping

- drop:

  logical indicating if levels that do not occur should be dropped (if
  'f' is a 'factor' or a list).

- ...:

  further potential arguments passed to methods.

## Value

A list of `xcmsSet` objects.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
