# Correct retention time from different samples

To correct differences between retention times between different
samples, a number of of methods exist in XCMS. `retcor` is the generic
method.

## Methods

- object = "xcmsSet":

  ` retcor(object, ...) `

## Arguments

- object:

  [`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  object

- method:

  Method to use for retention time correction. See details.

- ...:

  Optional arguments to be passed along

## Details

Different algorithms can be used by specifying them with the `method`
argument. For example to use the approach described by Smith et al
(2006) one would use: `retcor(object, method="loess")`. This is also the
default.

Further arguments given by `...` are passed through to the function
implementing the `method`.

A character vector of *nicknames* for the algorithms available is
returned by `getOption("BioC")$xcms$retcor.methods`. If the nickname of
a method is called "loess", the help page for that specific method can
be accessed with
[`?retcor.loess`](https://sneumann.github.io/xcms/reference/retcor.peakgroups-methods.md).

## Value

An `xcmsSet` object with corrected retntion times.

## See also

[`retcor.loess`](https://sneumann.github.io/xcms/reference/retcor.peakgroups-methods.md)
[`retcor.obiwarp`](https://sneumann.github.io/xcms/reference/retcor.obiwarp-methods.md)
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
