# Plot feature groups in the m/z-retention time space

`plotFeatureGroups()` visualizes defined feature groups in the m/z by
retention time space. Features are indicated by points with features
from the same feature group being connected by a line. See
[`MsFeatures::featureGroups()`](https://rdrr.io/pkg/MsFeatures/man/featureGroups.html)
for details on and options for feature grouping.

## Usage

``` r
plotFeatureGroups(
  x,
  xlim = numeric(),
  ylim = numeric(),
  xlab = "retention time",
  ylab = "m/z",
  pch = 4,
  col = "#00000060",
  type = "o",
  main = "Feature groups",
  featureGroups = character(),
  ...
)
```

## Arguments

- x:

  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with grouped features (i.e. after calling
  [`MsFeatures::groupFeatures()`](https://rdrr.io/pkg/MsFeatures/man/groupFeatures.html).

- xlim:

  `numeric(2)` with the lower and upper limit for the x-axis.

- ylim:

  `numeric(2)` with the lower and upper limit for the y-axis.

- xlab:

  `character(1)` with the label for the x-axis.

- ylab:

  `character(1)` with the label for the y-axis.

- pch:

  the plotting character. Defaults to `pch = 4` i.e. plotting features
  as crosses. See [`par()`](https://rdrr.io/r/graphics/par.html) for
  more information.

- col:

  color to be used to draw the features. At present only a single color
  is supported.

- type:

  plotting type (see [`par()`](https://rdrr.io/r/graphics/par.html)).
  Defaults to `type = "o"` which draws each feature as a point and
  connecting the features of the same feature group with a line.

- main:

  `character(1)` with the title of the plot.

- featureGroups:

  optional `character` of feature group IDs to draw only specified
  feature group(s). If not provided, all feature groups are drawn.

- ...:

  additional parameters to be passed to the `lines` function.

## Author

Johannes Rainer
