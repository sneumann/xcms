# Determine which peaks are absent / present in a sample class

Determine which peaks are absent / present in a sample class

## Methods

- object = "xcmsSet":

  ` absent(object, ...) present(object, ...) `

## Arguments

- object:

  [`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  object

- class:

  Name of a sample class from
  [`sampclass`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)

- minfrac:

  minimum fraction of samples necessary in the class to be
  absent/present

## Details

Determine which peaks are absent / present in a sample class The
functions treat peaks that are only present because of
[`fillPeaks`](https://sneumann.github.io/xcms/reference/fillPeaks-methods.md)
correctly, i.e. does not count them as present.

## Value

An logical vector with the same length as `nrow(groups(object))`.

## See also

[`group`](https://sneumann.github.io/xcms/reference/group-methods.md)
[`diffreport`](https://sneumann.github.io/xcms/reference/diffreport-methods.md)
