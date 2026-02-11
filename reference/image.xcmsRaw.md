# Plot log intensity image of a xcmsRaw object

Create log intensity false-color image of a xcmsRaw object plotted with
m/z and retention time axes

## Methods

- x = "xcmsRaw":

  ` image(x, col = rainbow(256), ...) `

## Arguments

- x:

  xcmsRaw object

- col:

  vector of colors to use for for the image

- ...:

  arguments for `profRange`

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
