# Integrate areas of missing peaks

For each sample, identify peak groups where that sample is not
represented. For each of those peak groups, integrate the signal in the
region of that peak group and create a new peak.

## Methods

- object = "xcmsSet":

  `fillPeaks(object, method="")`

## Arguments

- object:

  the `xcmsSet` object

- method:

  the filling method

## Details

After peak grouping, there will always be peak groups that do not
include peaks from every sample. This method produces intensity values
for those missing samples by integrating raw data in peak group region.
According to the type of raw-data there are 2 different methods
available. for filling gcms/lcms data the method "chrom" integrates
raw-data in the chromatographic domain, whereas "MSW" is used for
peaklists without retention-time information like those from
direct-infusion spectra.

## Value

A `xcmsSet` objects with filled in peak groups.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`getPeaks`](https://sneumann.github.io/xcms/reference/getPeaks-methods.md)
