# DEPRECATED: Create a plot that combines a XIC and a mz/rt 2D plot for one sample

**UPDATE**: please use [`plot()`](https://rdrr.io/r/base/plot.html) from
the `MsExperiment` or `plot(x, type = "XIC")` from the `MSnbase` package
instead. See examples in the vignette for more information.

The `plotMsData` creates a plot that combines an (base peak ) extracted
ion chromatogram on top (rt against intensity) and a plot of rt against
m/z values at the bottom.

## Usage

``` r
plotMsData(
  x,
  main = "",
  cex = 1,
  mfrow = c(2, 1),
  grid.color = "lightgrey",
  colramp = colorRampPalette(rev(brewer.pal(9, "YlGnBu")))
)
```

## Arguments

- x:

  `data.frame` such as returned by the
  [`extractMsData()`](https://sneumann.github.io/xcms/reference/extractMsData-method.md)
  function. Only a single `data.frame` is supported.

- main:

  `character(1)` specifying the title.

- cex:

  `numeric(1)` defining the size of points. Passed directly to the
  `plot` function.

- mfrow:

  `numeric(2)` defining the plot layout. This will be passed directly to
  `par(mfrow = mfrow)`. See `par` for more information. Setting
  `mfrow = NULL` avoids calling `par(mfrow = mfrow)` hence allowing to
  pre-define the plot layout.

- grid.color:

  a color definition for the grid line (or `NA` to skip creating them).

- colramp:

  a *color ramp palette* to be used to color the data points based on
  their intensity. See argument `col.regions` in
  [lattice::level.colors](https://rdrr.io/pkg/lattice/man/level.colors.html)
  documentation.

## Author

Johannes Rainer
