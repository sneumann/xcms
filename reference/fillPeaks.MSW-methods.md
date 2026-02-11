# Integrate areas of missing peaks in FTICR-MS data

For each sample, identify peak groups where that sample is not
represented. For each of those peak groups, integrate the signal in the
region of that peak group and create a new peak.

## Methods

- object = "xcmsSet":

  `fillPeaks.MSW(object)`

## Arguments

- object:

  the `xcmsSet` object

## Details

After peak grouping, there will always be peak groups that do not
include peaks from every sample. This method produces intensity values
for those missing samples by integrating raw data in peak group region.
In a given group, the start and ending m/z values for integration are
defined by the median start and end points of the other detected peaks.

## Note

In contrast to the
[`fillPeaks.chrom`](https://sneumann.github.io/xcms/reference/fillPeaks.chrom-methods.md)
method the maximum intensity reported in column `"maxo"` is not the
maximum intensity measured in the expected peak area (defined by columns
`"mzmin"` and `"mzmax"`), but the largest intensity of mz value(s)
closest to the `"mzmed"` of the feature.

## Value

A `xcmsSet` objects with filled in peak groups.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`getPeaks`](https://sneumann.github.io/xcms/reference/getPeaks-methods.md)
[`fillPeaks`](https://sneumann.github.io/xcms/reference/fillPeaks-methods.md)
