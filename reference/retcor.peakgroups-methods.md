# Align retention times across samples

These two methods use “well behaved” peak groups to calculate retention
time deviations for every time point of each sample. Use smoothed
deviations to align retention times.

## Methods

- object = "xcmsSet":

  `retcor(object, missing = 1, extra = 1, smooth = c("loess", "linear"), span = .2, family = c("gaussian", "symmetric"), plottype = c("none", "deviation", "mdevden"), col = NULL, ty = NULL)`

## Arguments

- object:

  the `xcmsSet` object

- missing:

  number of missing samples to allow in retention time correction groups

- extra:

  number of extra peaks to allow in retention time correction correction
  groups

- smooth:

  either `"loess"` for non-linear alignment or `"linear"` for linear
  alignment

- span:

  degree of smoothing for local polynomial regression fitting

- family:

  if `gaussian` fitting is by least-squares with no outlier removal, and
  if `symmetric` a re-descending M estimator is used with Tukey's
  biweight function, allowing outlier removal

- plottype:

  if `deviation` plot retention time deviation points and regression
  fit, and if `mdevden` also plot peak overall peak density and
  retention time correction peak density

- col:

  vector of colors for plotting each sample

- ty:

  vector of line and point types for plotting each sample

## Value

An `xcmsSet` object

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`loess`](https://rdrr.io/r/stats/loess.html)
[`retcor.obiwarp`](https://sneumann.github.io/xcms/reference/retcor.obiwarp-methods.md)
