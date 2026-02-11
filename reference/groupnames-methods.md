# Generate unque names for peak groups

Allow linking of peak group data between classes using unique group
names that remain the same as long as no re-grouping occurs.

## Methods

- object = "xcmsSet":

  `(object, mzdec = 0, rtdec = 0, template = NULL)`

- object = "xcmsEIC":

  `(object)`

## Arguments

- object:

  the `xcmsSet` or `xcmsEIC` object

- mzdec:

  number of decimal places to use for m/z

- rtdec:

  number of decimal places to use for retention time

- template:

  a character vector with existing group names whose format should be
  emulated

## Value

A character vector with unique names for each peak group in the object.
The format is `M[m/z]T[time in seconds]`.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`xcmsEIC-class`](https://sneumann.github.io/xcms/reference/xcmsEIC-class.md)
