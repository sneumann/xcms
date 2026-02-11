# Plot profile matrix 3D surface using OpenGL

This method uses the rgl package to create interactive three dimensonal
representations of the profile matrix. It uses the terrain color scheme.

## Methods

- object = "xcmsRaw":

  `plotSurf(object, log = FALSE, aspect = c(1, 1, .5), ...)`

## Arguments

- object:

  the `xcmsRaw` object

- log:

  logical, log transform intensity

- aspect:

  numeric vector with aspect ratio of the m/z, retention time and
  intensity components of the plot

- ...:

  arguments passed to
  [`profRange`](https://sneumann.github.io/xcms/reference/profRange-methods.md)

## Details

The rgl package is still in development and imposes some limitations on
the output format. A bug in the axis label code means that the axis
labels only go from 0 to the aspect ratio constant of that axis.
Additionally the axes are not labeled with what they are.

It is important to only plot a small portion of the profile matrix.
Large portions can quickly overwhelm your CPU and memory.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
