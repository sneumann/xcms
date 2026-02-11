# Group peaks from different samples together

A number of grouping (or alignment) methods exist in XCMS. `group` is
the generic method.

## Methods

- object = "xcmsSet":

  ` group(object, ...) `

## Arguments

- object:

  [`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  object

- method:

  Method to use for grouping. See details.

- ...:

  Optional arguments to be passed along

## Details

Different algorithms can be used by specifying them with the `method`
argument. For example to use the density-based approach described by
Smith et al (2006) one would use: `group(object, method="density")`.
This is also the default.

Further arguments given by `...` are passed through to the function
implementing the `method`.

A character vector of *nicknames* for the algorithms available is
returned by `getOption("BioC")$xcms$group.methods`. If the nickname of a
method is called "mzClust", the help page for that specific method can
be accessed with
[`?group.mzClust`](https://sneumann.github.io/xcms/reference/group.mzClust.md).

## Value

An `xcmsSet` object with peak group assignments and statistics.

## See also

[`group.density`](https://sneumann.github.io/xcms/reference/group.density.md)
[`group.mzClust`](https://sneumann.github.io/xcms/reference/group.mzClust.md)
[`group.nearest`](https://sneumann.github.io/xcms/reference/group.nearest.md)
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
