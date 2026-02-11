# Distance methods for xcmsSet, xcmsRaw and xsAnnotate

There are several methods for calculating a distance between two sets of
peaks in xcms. `specDist` is the generic method.

## Methods

- object = "xcmsSet":

  ` specDist(object, peakIDs1, peakIDs2,...) `

&nbsp;

- object = "xsAnnotate":

  ` specDist(object, PSpec1, PSpec2,...) `

## Arguments

- object:

  a xcmsSet or xcmsRaw.

- method:

  Method to use for distance calculation. See details.

- ...:

  mzabs, mzppm and parameters for the distance function.

## Details

Different algorithms can be used by specifying them with the `method`
argument. For example to use the "meanMZmatch" approach with xcmsSet one
would use: `specDist(object, peakIDs1, peakIDs2, method="meanMZmatch")`.
This is also the default.

Further arguments given by `...` are passed through to the function
implementing the `method`.

A character vector of *nicknames* for the algorithms available is
returned by `getOption("BioC")$xcms$specDist.methods`. If the nickname
of a method is called "meanMZmatch", the help page for that specific
method can be accessed with
[`?specDist.meanMZmatch`](https://sneumann.github.io/xcms/reference/specDist.meanMZmatch-methods.md).

## Value

- mzabs:

  maximum absolute deviation for two matching peaks

- mzppm:

  relative deviations in ppm for two matching peaks

- symmetric:

  use symmetric pairwise m/z-matches only, or each match

## Author

Joachim Kutzera, <jkutzer@ipb-halle.de>
