# Combine xcmsSet objects

Combines the samples and peaks from multiple `xcmsSet` objects into a
single object. Group and retention time correction data are discarded.
The `profinfo` list is set to be equal to the first object.

## Methods

- xs1 = "xcmsRaw":

  ` c(xs1, ...) `

## Arguments

- xs1:

  `xcmsSet` object

- ...:

  `xcmsSet` objects

## Value

A `xcmsSet` object.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
