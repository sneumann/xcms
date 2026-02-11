# Feature detection for GC/MS and LC/MS Data - methods

A number of peak pickers exist in XCMS. `findPeaks` is the generic
method.

## Methods

- object = "xcmsRaw":

  ` findPeaks(object, ...) `

## Arguments

- object:

  [`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
  object

- method:

  Method to use for peak detection. See details.

- ...:

  Optional arguments to be passed along

## Details

Different algorithms can be used by specifying them with the `method`
argument. For example to use the matched filter approach described by
Smith et al (2006) one would use:
`findPeaks(object, method="matchedFilter")`. This is also the default.

Further arguments given by `...` are passed through to the function
implementing the `method`.

A character vector of *nicknames* for the algorithms available is
returned by `getOption("BioC")$xcms$findPeaks.methods`. If the nickname
of a method is called "centWave", the help page for that specific method
can be accessed with
[`?findPeaks.centWave`](https://sneumann.github.io/xcms/reference/findPeaks.centWave-methods.md).

## Value

A matrix with columns:

- mz:

  weighted (by intensity) mean of peak m/z across scans

- mzmin:

  m/z of minimum step

- mzmax:

  m/z of maximum step

- rt:

  retention time of peak midpoint

- rtmin:

  leading edge of peak retention time

- rtmax:

  trailing edge of peak retention time

- into:

  integrated area of original (raw) peak

- maxo:

  maximum intensity of original (raw) peak

and additional columns depending on the choosen method.

## See also

[`findPeaks.matchedFilter`](https://sneumann.github.io/xcms/reference/findPeaks.matchedFilter-xcmsRaw-method.md)
[`findPeaks.centWave`](https://sneumann.github.io/xcms/reference/findPeaks.centWave-methods.md)
[`findPeaks.addPredictedIsotopeFeatures`](https://sneumann.github.io/xcms/reference/findPeaks.addPredictedIsotopeFeatures-methods.md)
[`findPeaks.centWaveWithPredictedIsotopeROIs`](https://sneumann.github.io/xcms/reference/findPeaks.centWaveWithPredictedIsotopeROIs-methods.md)
[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
