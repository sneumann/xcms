# Feature detection for high resolution LC/MS data

Peak density and wavelet based feature detection for high resolution
LC/MS data in centroid mode

## Methods

- object = "xcmsRaw":

  ` findPeaks.centWave(object, ppm=25, peakwidth=c(20,50), snthresh=10, prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001, fitgauss=FALSE, scanrange= numeric(), noise=0, sleep=0, verbose.columns=FALSE, ROI.list=list()), firstBaselineCheck=TRUE, roiScales=NULL `

## Details

This algorithm is most suitable for high resolution
LC/{TOF,OrbiTrap,FTICR}-MS data in centroid mode. In the first phase of
the method mass traces (characterised as regions with less than `ppm`
m/z deviation in consecutive scans) in the LC/MS map are located. In the
second phase these mass traces are further analysed. Continuous wavelet
transform (CWT) is used to locate chromatographic peaks on different
scales.

## Arguments

- object:

  `xcmsSet` object

- ppm:

  maxmial tolerated m/z deviation in consecutive scans, in ppm (parts
  per million)

- peakwidth:

  Chromatographic peak width, given as range (min,max) in seconds

- snthresh:

  signal to noise ratio cutoff, definition see below.

- prefilter:

  `prefilter=c(k,I)`. Prefilter step for the first phase. Mass traces
  are only retained if they contain at least `k` peaks with intensity
  \>= `I`.

- mzCenterFun:

  Function to calculate the m/z center of the feature: `wMean` intensity
  weighted mean of the feature m/z values, `mean` mean of the feature
  m/z values, `apex` use m/z value at peak apex, `wMeanApex3` intensity
  weighted mean of the m/z value at peak apex and the m/z value left and
  right of it, `meanApex3` mean of the m/z value at peak apex and the
  m/z value left and right of it.

- integrate:

  Integration method. If `=1` peak limits are found through descent on
  the mexican hat filtered data, if `=2` the descent is done on the real
  data. Method 2 is very accurate but prone to noise, while method 1 is
  more robust to noise but less exact.

- mzdiff:

  minimum difference in m/z for peaks with overlapping retention times,
  can be negative to allow overlap

- fitgauss:

  logical, if TRUE a Gaussian is fitted to each peak

- scanrange:

  scan range to process

- noise:

  optional argument which is useful for data that was centroided without
  any intensity threshold, centroids with intensity \< `noise` are
  omitted from ROI detection

- sleep:

  number of seconds to pause between plotting peak finding cycles

- verbose.columns:

  logical, if TRUE additional peak meta data columns are returned

- ROI.list:

  A optional list of ROIs that represents detected mass traces (ROIs).
  If this list is empty (default) then centWave detects the mass trace
  ROIs, otherwise this step is skipped and the supplied ROIs are used in
  the peak detection phase. Each ROI object in the list has the
  following slots: `scmin` start scan index, `scmax` end scan index,
  `mzmin` minimum m/z, `mzmax` maximum m/z, `length` number of scans,
  `intensity` summed intensity.

- firstBaselineCheck:

  logical, if TRUE continuous data within ROI is checked to be above 1st
  baseline

- roiScales:

  numeric, optional vector of scales for each ROI in `ROI.list` to be
  used for the centWave-wavelets

## Value

A matrix with columns:

- mz:

  weighted (by intensity) mean of peak m/z across scans

- mzmin:

  m/z peak minimum

- mzmax:

  m/z peak maximum

- rt:

  retention time of peak midpoint

- rtmin:

  leading edge of peak retention time

- rtmax:

  trailing edge of peak retention time

- into:

  integrated peak intensity

- intb:

  baseline corrected integrated peak intensity

- maxo:

  maximum peak intensity

- sn:

  Signal/Noise ratio, defined as `(maxo - baseline)/sd`, where\
  `maxo` is the maximum peak intensity,\
  `baseline` the estimated baseline value and\
  `sd` the standard deviation of local chromatographic noise.

- egauss:

  RMSE of Gaussian fit

if `verbose.columns` is `TRUE` additionally :

- mu:

  Gaussian parameter mu

- sigma:

  Gaussian parameter sigma

- h:

  Gaussian parameter h

- f:

  Region number of m/z ROI where the peak was localised

- dppm:

  m/z deviation of mass trace across scans in ppm

- scale:

  Scale on which the peak was localised

- scpos:

  Peak position found by wavelet analysis

- scmin:

  Left peak limit found by wavelet analysis (scan number)

- scmax:

  Right peak limit found by wavelet analysis (scan number)

## Author

Ralf Tautenhahn

## References

Ralf Tautenhahn, Christoph Böttcher, and Steffen Neumann "Highly
sensitive feature detection for high resolution LC/MS" BMC
Bioinformatics 2008, 9:504

## See also

[`centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md)
for the new user interface.
[`findPeaks-methods`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
