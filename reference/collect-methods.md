# Collect MS^n peaks into xcmsFragments

Collecting Peaks into
[`xcmsFragments`](https://sneumann.github.io/xcms/reference/xcmsFragments-class.md)s
from several MS-runs using
[`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
and
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md).

## Methods

- object = "xcmsFragments":

  ` collect(object, ...) `

## Arguments

- object:

  (empty)
  [`xcmsFragments-class`](https://sneumann.github.io/xcms/reference/xcmsFragments-class.md)
  object

- xs:

  A
  [`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  object which contains picked ms1-peaks from several experiments

- compMethod:

  ("floor", "round", "none"): compare-method which is used to find the
  parent peak of a MSnpeak through comparing the MZ-values of the
  MS1peaks with the MSnParentPeaks.

- snthresh, mzgap, uniq:

  these are the parameters for the getspec-peakpicker included in
  xcmsRaw.

## Details

After running collect(xFragments,xSet) The peak table of the
xcmsFragments includes the ms1Peaks from all experiments stored in a
xcmsSet-object. Further it contains the relevant msN-peaks from the
xcmsRaw-objects, which were created temporarily with the paths in
xcmsSet.

## Value

A matrix with columns:

- peakID:

  unique identifier of every peak

- MSnParentPeakID:

  PeakID of the parent peak of a msLevel\>1 - peak, it is 0 if the peak
  is msLevel 1.

- msLevel:

  The msLevel of the peak.

- rt:

  retention time of the peak midpoint

- mz:

  the mz-Value of the peak

- intensity:

  the intensity of the peak

- sample:

  the number of the sample from the xcmsSet

- GroupPeakMSn:

  Used for grouped xcmsSet groups

- CollisionEnergy:

  The collision energy of the fragment
