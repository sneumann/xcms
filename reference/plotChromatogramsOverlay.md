# Plot multiple chromatograms into the same plot

`plotOverlay` draws chromatographic peak data from multiple (different)
extracted ion chromatograms (EICs) into the same plot. This allows to
directly compare the peak shape of these EICs in the same sample. In
contrast to the `plot` function for
[`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
object, which draws the data from the same EIC across multiple samples
in the same plot, this function draws the different EICs from the same
sample into the same plot.

If `plotChromatogramsOverlay` is called on a `XChromatograms` object any
present chromatographic peaks will also be highlighted/drawn depending
on the parameters `peakType`, `peakCol`, `peakBg` and `peakPch` (see
also help on the `plot` function for
[`XChromatogram()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
object for details).

## Usage

``` r
# S4 method for class 'MChromatograms'
plotChromatogramsOverlay(
  object,
  col = "#00000060",
  type = "l",
  main = NULL,
  xlab = "rtime",
  ylab = "intensity",
  xlim = numeric(),
  ylim = numeric(),
  stacked = 0,
  transform = identity,
  ...
)

# S4 method for class 'XChromatograms'
plotChromatogramsOverlay(
  object,
  col = "#00000060",
  type = "l",
  main = NULL,
  xlab = "rtime",
  ylab = "intensity",
  xlim = numeric(),
  ylim = numeric(),
  peakType = c("polygon", "point", "rectangle", "none"),
  peakBg = NULL,
  peakCol = NULL,
  peakPch = 1,
  stacked = 0,
  transform = identity,
  ...
)
```

## Arguments

- object:

  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  or
  [`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  object.

- col:

  definition of the color in which the chromatograms should be drawn.
  Can be of length 1 or equal to `nrow(object)` to plot each overlayed
  chromatogram in a different color.

- type:

  `character(1)` defing the type of the plot. By default (`type = "l"`)
  each chromatogram is drawn as a line.

- main:

  optional title of the plot. If not defined, the range of m/z values is
  used.

- xlab:

  `character(1)` defining the x-axis label.

- ylab:

  `character(1)` defining the y-axis label.

- xlim:

  optional `numeric(2)` defining the x-axis limits.

- ylim:

  optional `numeric(2)` defining the y-axis limits.

- stacked:

  `numeric(1)` defining the part (proportion) of the y-axis to use to
  *stack* EICs depending on their m/z values. If `stacked = 0` (the
  default) no stacking is performed. With `stacked = 1` half of the
  y-axis is used for stacking and half for the intensity y-axis (i.e.
  the ratio between stacking and intensity y-axis is 1:1). Note that if
  `stacking` is different from 0 no y-axis and label are drawn.

- transform:

  `function` to transform the intensity values before plotting. Defaults
  to `transform = identity` which plots the data as it is. With
  `transform = log10` intensity values would be log10 transformed before
  plotting.

- ...:

  optional arguments to be passed to the plotting functions (see help on
  the base R `plot` function.

- peakType:

  if `object` is a `XChromatograms` object: how chromatographic peaks
  should be drawn: `peakType = "polygon"` (the default): label the full
  chromatographic peak area, `peakType = "rectangle"`: indicate the
  chromatographic peak by a rectangle and `peakType = "point"`: label
  the chromatographic peaks' apex position with a point.

- peakBg:

  if `object` is a `XChromatograms` object: definition of background
  color(s) for each chromatographic peak. Has to be either of length 1
  or equal to the number of peaks in `object`. If not specified, the
  peak will be drawn in the color defined by `col`.

- peakCol:

  if `object` is a `XChromatograms` object: definition of color(s) for
  each chromatographic peak. Has to be either of length 1 or equal to
  the number of peaks in `object`. If not specified, the peak will be
  drawn in the color defined by `col`.

- peakPch:

  if `object` is a `XChromatograms` object: *point character* to be used
  to label the apex position of the chromatographic peak if
  `peakType = "point"`.

## Value

silently returns a `list` (length equal to `ncol(object)` of `numeric`
(length equal to `nrow(object)`) with the y position of each EIC.

## Author

Johannes Rainer

## Examples

``` r

## Load preprocessed data and extract EICs for some features.
library(xcms)
library(MSnbase)
xdata <- loadXcmsData()
data(xdata)
## Update the path to the files for the local system
dirname(xdata) <- c(rep(system.file("cdf", "KO", package = "faahKO"), 4),
                    rep(system.file("cdf", "WT", package = "faahKO"), 4))
#> Error: unable to find an inherited method for function ‘path’ for signature ‘object = "XcmsExperiment"’
## Subset to the first 3 files.
xdata <- filterFile(xdata, 1:3, keepFeatures = TRUE)

## Define features for which to extract EICs
fts <- c("FT097", "FT163", "FT165")
chrs <- featureChromatograms(xdata, features = fts)
#> Processing chromatographic peaks for features

plotChromatogramsOverlay(chrs)


## plot the overlay of EICs in the first sample
plotChromatogramsOverlay(chrs[, 1])

## Define a different color for each feature (row in chrs). By default, also
## all chromatographic peaks of a feature is labeled in the same color.
plotChromatogramsOverlay(chrs[, 1],
    col = c("#ff000040", "#00ff0040", "#0000ff40"))

## Alternatively, we can define a color for each individual chromatographic
## peak and provide this with the `peakBg` and `peakCol` parameters.
chromPeaks(chrs[, 1])
#>              mz mzmin mzmax       rt    rtmin    rtmax     into     intb  maxo
#> CP2975 354.1357 354.1 354.2 2790.288 2755.843 2826.772 298657.9       NA  9700
#> CP0312 413.1000 413.1 413.1 3643.850 3620.732 3672.298 879851.5 852195.3 30920
#> CP0089 414.2000 414.2 414.2 3066.209 3038.494 3089.233 234308.5 218923.4  9669
#>        sn sample row column
#> CP2975 NA      1   1      1
#> CP0312 34      1   2      1
#> CP0089 15      1   3      1

## Use a color for each of the two identified peaks in that sample
plotChromatogramsOverlay(chrs[, 1],
    col = c("#ff000040", "#00ff0040", "#0000ff40"),
    peakBg = c("#ffff0020", "#00ffff20"))


## Plotting the data in all samples.
plotChromatogramsOverlay(chrs,
    col = c("#ff000040", "#00ff0040", "#0000ff40"))


## Creating a "stacked" EIC plot: the EICs are placed along the y-axis
## relative to their m/z value. With `stacked = 1` the y-axis is split in
## half, the lower half being used for the stacking of the EICs, the upper
## half being used for the *original* intensity axis.
res <- plotChromatogramsOverlay(chrs[, 1], stacked = 1,
    col = c("#ff000040", "#00ff0040", "#0000ff40"))
## add horizontal lines for the m/z values of each EIC
abline(h = res[[1]], col = "grey", lty = 2)

## Note that this type of visualization is different than the conventional
## plot function for chromatographic data, which will draw the EICs for
## multiple samples into the same plot
plot(chrs)


## Converting the object to a MChromatograms without detected peaks
chrs <- as(chrs, "MChromatograms")

plotChromatogramsOverlay(chrs,
    col = c("#ff000040", "#00ff0040", "#0000ff40"))

```
