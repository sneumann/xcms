# Extract a matrix of peak values for each group

Generate a matrix of peak values with rows for every group and columns
for every sample. The value included in the matrix can be any of the
columns from the `xcmsSet` `peaks` slot matrix. Collisions where more
than one peak from a single sample are in the same group get resolved
with one of several user-selectable methods.

## Methods

- object = "xcmsSet":

  `groupval(object, method = c("medret", "maxint"), value = "index", intensity = "into")`

## Arguments

- object:

  the `xcmsSet` object

- method:

  conflict resolution method, `"medret"` to use the peak closest to the
  median retention time or `"maxint"` to use the peak with the highest
  intensity

- value:

  name of peak column to enter into returned matrix, or `"index"` for
  index to the corresponding row in the `peaks` slot matrix

- intensity:

  if `method == "maxint"`, name of peak column to use for intensity

## Value

A matrix with with rows for every group and columns for every sample.
Missing peaks have `NA` values.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
