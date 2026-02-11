# Integrate areas of missing peaks

For each sample, identify peak groups where that sample is not
represented. For each of those peak groups, integrate the signal in the
region of that peak group and create a new peak.

## Methods

- object = "xcmsSet":

  `fillPeaks.chrom(object, nSlaves=0,expand.mz=1,expand.rt=1, BPPARAM = bpparam())`

## Arguments

- object:

  the `xcmsSet` object

- nSlaves:

  (DEPRECATED): number of slaves/cores to be used for parallel peak
  filling. MPI is used if installed, otherwise the snow package is
  employed for multicore support. If none of the two packages is
  available it uses the parallel package for parallel processing on
  multiple CPUs of the current machine. Users are advised to use the
  `BPPARAM` parameter instead.

- expand.mz:

  Expansion factor for the m/z range used for integration.

- expand.rt:

  Expansion factor for the rentention time range used for integration.

- BPPARAM:

  allows to define a specific parallel processing setup for the current
  task (see
  [`bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html) from
  the `BiocParallel` package help more information). The default uses
  the globally defined parallel setup.

## Details

After peak grouping, there will always be peak groups that do not
include peaks from every sample. This method produces intensity values
for those missing samples by integrating raw data in peak group region.
In a given group, the start and ending retention time points for
integration are defined by the median start and end points of the other
detected peaks. The start and end m/z values are similarly determined.
Intensities can be still be zero, which is a rather unusual intensity
for a peak. This is the case if e.g. the raw data was threshholded, and
the integration area contains no actual raw intensities, or if one
sample is miscalibrated, such thet the raw data points are (just)
outside the integration area.

Importantly, if retention time correction data is available, the
alignment information is used to more precisely integrate the propper
region of the raw data. If the corrected retention time is beyond the
end of the raw data, the value will be not-a-number (NaN).

## Value

A `xcmsSet` objects with filled in peak groups (into and maxo).

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`getPeaks`](https://sneumann.github.io/xcms/reference/getPeaks-methods.md)
[`fillPeaks`](https://sneumann.github.io/xcms/reference/fillPeaks-methods.md)
