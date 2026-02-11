# Load the raw data for one or more files in the xcmsSet

Reads the raw data applies evential retention time corrections and
waters Lock mass correction and returns it as an `xcmsRaw` object (or
list of `xcmsRaw` objects) for one or more files of the `xcmsSet`
object.

## Methods

- object = "xcmsSet":

  `getXcmsRaw(object, sampleidx=1, profmethod=profinfo(object)$method, profstep=profinfo(object)$step, rt=c("corrected", "raw"), ... )`

## Arguments

- object:

  the `xcmsSet` object

- sampleidx:

  The index of the sample for which the raw data should be returned. Can
  be a single number or a numeric vector with the indices.
  Alternatively, the file name can be specified.

- profmethod:

  The profile method.

- profstep:

  The profile step.

- rt:

  Whether corrected or raw retention times should be returned.

- ...:

  Additional arguments submitted to the
  [`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md)
  function.

## Value

A single `xcmsRaw` object or a list of `xcmsRaw` objects.

## Author

Johannes Rainer, <johannes.rainer@eurac.edu>

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
