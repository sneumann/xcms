# Plot log intensity image of a xcmsRaw object

Create an image of the raw (profile) data m/z against retention time,
with the intensity color coded.

## Methods

- x = "xcmsRaw":

  ` levelplot(x, log=TRUE, col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(256), ...) `

- x = "xcmsSet":

  ` levelplot(x, log=TRUE, col.regions=colorRampPalette(brewer.pal(9, "YlOrRd"))(256), rt="raw", ...) `

## Arguments

- x:

  xcmsRaw object.

- log:

  Whether the intensity should be log transformed.

- col.regions:

  The color ramp that should be used for encoding of the intensity.

- rt:

  wheter the original (`rt="raw"`) or the corrected (`rt="corrected"`)
  retention times should be used.

- ...:

  Arguments for `profRange`.

## Author

Johannes Rainer, <johannes.rainer@eurac.edu>

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
