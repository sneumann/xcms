# Class xcmsRaw, a class for handling raw data

This class handles processing and visualization of the raw data from a
single LC/MS or GS/MS run. It includes methods for producing a standard
suite of plots including individual spectra, multi-scan average spectra,
TIC, and EIC. It will also produce a feature list of significant peaks
using matched filtration.

## Objects from the Class

Objects can be created with the
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md)
constructor which reads data from a NetCDF file into a new object.

## Slots

- `acquisitionNum`::

  Numeric representing the acquisition number of the individual
  scans/spectra. Length of `acquisitionNum` is equal to the number of
  spectra/scans in the object and hence equal to the `scantime` slot.
  Note however that this information is only available in mzML files.

- `env`::

  environment with three variables: `mz` - concatenated m/z values for
  all scans, `intensity` - corresponding signal intensity for each m/z
  value, and `profile` - matrix represention of the intensity values
  with columns representing scans and rows representing equally spaced
  m/z values. The profile matrix should be extracted with the
  [`profMat`](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md)
  method.

- `filepath`::

  Path to the raw data file

- `gradient`::

  matrix with first row, `time`, containing the time point for
  interpolation and successive columns representing solvent fractions at
  each point

- `msnAcquisitionNum`::

  for each scan a unique acquisition number as reported via "spectrum
  id" (mzData) or "\<scan num=...\>" and "\<scanOrigin num=...\>"
  (mzXML)

- `msnCollisionEnergy`::

  "CollisionEnergy" (mzData) or "collisionEnergy" (mzXML)

- `msnLevel`::

  for each scan the "msLevel" (both mzData and mzXML)

- `msnPrecursorCharge`::

  "ChargeState" (mzData) and "precursorCharge" (mzXML)

- `msnPrecursorIntensity`::

  "Intensity" (mzData) or "precursorIntensity" (mzXML)

- `msnPrecursorMz`::

  "MassToChargeRatio" (mzData) or "precursorMz" (mzXML)

- `msnPrecursorScan`::

  "spectrumRef" (both mzData and mzXML)

- `msnRt`::

  Retention time of the scan

- `msnScanindex`::

  msnScanindex

- `mzrange`::

  numeric vector of length 2 with minimum and maximum m/z values
  represented in the profile matrix

- `polarity`::

  polarity

- `profmethod`::

  characer value with name of method used for generating the profile
  matrix.

- `profparam`::

  list to store additional profile matrix generation settings. Use the
  [`profinfo`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  method to extract all profile matrix creation relevant information.

- `scanindex`::

  integer vector with starting positions of each scan in the `mz` and
  `intensity` variables (note that index values are based off a 0
  initial position instead of 1).

- `scantime`::

  numeric vector with acquisition time (in seconds) for each scan.

- `tic`::

  numeric vector with total ion count (intensity) for each scan

- `mslevel`::

  Numeric representing the MS level that is present in MS1 slot. This
  slot should be accessed through its getter method `mslevel`.

- `scanrange`::

  Numeric of length 2 specifying the scan range (or `NULL` for the full
  range). This slot should be accessed through its getter method
  `scanrange`. Note that the `scanrange` will always be 1 to the number
  of scans within the `xcmsRaw` object, which does not necessarily have
  to match to the scan index in the original mzML file (e.g. if the
  original data was sub-setted). The `acquisitionNum` information can be
  used to track the original *position* of each scan in the mzML file.

## Methods

- [findPeaks](https://sneumann.github.io/xcms/reference/findPeaks-methods.md):

  `signature(object = "xcmsRaw")`: feature detection using matched
  filtration in the chromatographic time domain

- [getEIC](https://sneumann.github.io/xcms/reference/getEIC-methods.md):

  `signature(object = "xcmsRaw")`: get extracted ion chromatograms in
  specified m/z ranges. This will return the total ion chromatogram
  (TIC) if the m/z range corresponds to the full m/z range (i.e. sum of
  all signals per retention time across all m/z).

- [getPeaks](https://sneumann.github.io/xcms/reference/getPeaks-methods.md):

  `signature(object = "xcmsRaw")`: get data for peaks in specified m/z
  and time ranges

- [getScan](https://sneumann.github.io/xcms/reference/getScan-methods.md):

  `signature(object = "xcmsRaw")`: get m/z and intensity values for a
  single mass scan

- [getSpec](https://sneumann.github.io/xcms/reference/getSpec-methods.md):

  `signature(object = "xcmsRaw")`: get average m/z and intensity values
  for multiple mass scans

- image:

  `signature(x = "xcmsRaw")`: get data for peaks in specified m/z and
  time ranges

- levelplot:

  Create an image of the raw (profile) data m/z against retention time,
  with the intensity color coded.

- mslevel:

  Getter method for the `mslevel` slot.

- [plotChrom](https://sneumann.github.io/xcms/reference/plotChrom-methods.md):

  `signature(object = "xcmsRaw")`: plot a chromatogram from profile data

- [plotRaw](https://sneumann.github.io/xcms/reference/plotRaw-methods.md):

  `signature(object = "xcmsRaw")`: plot locations of raw intensity data
  points

- [plotScan](https://sneumann.github.io/xcms/reference/plotScan-methods.md):

  `signature(object = "xcmsRaw")`: plot a mass spectrum of an individual
  scan from the raw data

- [plotSpec](https://sneumann.github.io/xcms/reference/plotSpec-methods.md):

  `signature(object = "xcmsRaw")`: plot a mass spectrum from profile
  data

- [plotSurf](https://sneumann.github.io/xcms/reference/plotSurf-methods.md):

  `signature(object = "xcmsRaw")`: experimental method for plotting 3D
  surface of profile data with `rgl`.

- [plotTIC](https://sneumann.github.io/xcms/reference/plotTIC-methods.md):

  `signature(object = "xcmsRaw")`: plot total ion count chromatogram

- [profinfo](https://sneumann.github.io/xcms/reference/xcmsSet-class.md):

  `signature(object = "xcmsRaw")`: returns a list containing the profile
  generation method and step (profile m/z step size) and eventual
  additional parameters to the profile function.

- [profMedFilt](https://sneumann.github.io/xcms/reference/profMedFilt-methods.md):

  `signature(object = "xcmsRaw")`: median filter profile data in time
  and m/z dimensions

- profMethod\<-:

  `signature(object = "xcmsRaw")`: change the method of generating the
  `profile` matrix

- [profMethod](https://sneumann.github.io/xcms/reference/profMethod-methods.md):

  `signature(object = "xcmsRaw")`: get the method of generating the
  `profile` matrix

- profMz:

  `signature(object = "xcmsRaw")`: get vector of m/z values for each row
  of the `profile` matrix

- [profRange](https://sneumann.github.io/xcms/reference/profRange-methods.md):

  `signature(object = "xcmsRaw")`: interpret flexible ways of specifying
  subsets of the `profile` matrix

- profStep\<-:

  `signature(object = "xcmsRaw")`: change the m/z step used for
  generating the `profile` matrix

- [profStep](https://sneumann.github.io/xcms/reference/profStep-methods.md):

  `signature(object = "xcmsRaw")`: get the m/z step used for generating
  the `profile` matrix

- revMz:

  `signature(object = "xcmsRaw")`: reverse the order of the data points
  for each scan

- scanrange:

  Getter method for the `scanrange` slot. See slot description above for
  more information.

- sortMz:

  `signature(object = "xcmsRaw")`: sort the data points by increasing
  m/z for each scan

- stitch:

  `signature(object = "xcmsRaw")`: Raw data correction for lock mass
  calibration gaps.

- findmzROI:

  `signature(object = "xcmsRaw")`: internal function to identify regions
  of interest in the raw data as part of the first step of
  centWave-based peak detection.

## Author

Colin A. Smith, <csmith@scripps.edu>, Johannes Rainer
<johannes.rainer@eurac.edu>

## See also

[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md),
[`subset-xcmsRaw`](https://sneumann.github.io/xcms/reference/sub-xcmsRaw-logicalOrNumeric-missing-missing-method.md)
for subsetting by spectra.
