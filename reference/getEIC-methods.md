# Get extracted ion chromatograms for specified m/z ranges

Generate multiple extracted ion chromatograms for m/z values of
interest. For `xcmsSet` objects, reread original raw data and apply
precomputed retention time correction, if applicable.

Note that this method will *always* return profile, not raw data (with
profile data being the binned data along M/Z). See details for further
information.

## Methods

- object = "xcmsRaw":

  `getEIC(object, mzrange, rtrange = NULL, step = 0.1)`

- object = "xcmsSet":

  `getEIC(object, mzrange, rtrange = 200, groupidx, sampleidx = sampnames(object), rt = c("corrected", "raw"))`

## Arguments

- object:

  the `xcmsRaw` or `xcmsSet` object

- mzrange:

  Either a two column matrix with minimum or maximum m/z or a matrix of
  any dimensions containing columns `mzmin` and `mzmax`. If not
  specified, the method for `xcmsRaw` returns the base peak chromatogram
  (BPC, i.e. the most intense signal for each RT across all m/z).

  For `xcmsSet` objects the group data will be used if `mzrange` is not
  provided.

- rtrange:

  A two column matrix the same size as `mzrange` with minimum and
  maximum retention times between which to return EIC data points. If
  not specified, the method returns the chromatogram for the full RT
  range.

  For `xcmsSet` objects, it may also be a single number specifying the
  time window around the peak to return EIC data points

- step:

  step (bin) size to use for profile generation. Note that a value of
  `step = 0` is not supported.

- groupidx:

  either character vector with names or integer vector with indicies of
  peak groups for which to get EICs

- sampleidx:

  either character vector with names or integer vector with indicies of
  samples for which to get EICs

- rt:

  `"corrected"` for using corrected retention times, or `"raw"` for
  using raw retention times

## Details

In contrast to the
[`rawEIC`](https://sneumann.github.io/xcms/reference/rawEIC-methods.md)
method, that extracts the actual raw values, this method extracts them
from the object's profile matrix (or if the provided `step` argument
does not match the `profStep` of the `object` the profile matrix is
calculated on the fly and the values returned).

## Value

For `xcmsSet` and `xcmsRaw` objects, an `xcmsEIC` object.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`xcmsEIC-class`](https://sneumann.github.io/xcms/reference/xcmsEIC-class.md),
[`rawEIC`](https://sneumann.github.io/xcms/reference/rawEIC-methods.md)
