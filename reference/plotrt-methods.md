# Plot retention time deviation profiles

Use corrected retention times for each sample to calculate retention
time deviation profiles and plot each on the same graph.

## Methods

- object = "xcmsSet":

  `plotrt(object, col = NULL, ty = NULL, leg = TRUE, densplit = FALSE)`

## Arguments

- object:

  the `xcmsSet` object

- col:

  vector of colors for plotting each sample

- ty:

  vector of line and point types for plotting each sample

- leg:

  logical plot legend with sample labels

- densplit:

  logical, also plot peak overall peak density

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`retcor`](https://sneumann.github.io/xcms/reference/retcor-methods.md)
