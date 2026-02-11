# Feature detection based on predicted isotope features for high resolution LC/MS data

Peak density and wavelet based feature detection aiming at isotope peaks
for high resolution LC/MS data in centroid mode

## Methods

- object = "xcmsRaw":

  ` findPeaks.centWave(object, ppm=25, peakwidth=c(20,50), prefilter=c(3,100), mzCenterFun="wMean", integrate=1, mzdiff=-0.001, fitgauss=FALSE, scanrange= numeric(), noise=0, sleep=0, verbose.columns=FALSE, xcmsPeaks, snthresh=6.25, maxcharge=3, maxiso=5, mzIntervalExtension=TRUE) `

## Details

This algorithm is most suitable for high resolution
LC/{TOF,OrbiTrap,FTICR}-MS data in centroid mode. In the first phase of
the method isotope ROIs (regions of interest) in the LC/MS map are
predicted. In the second phase these mass traces are further analysed.
Continuous wavelet transform (CWT) is used to locate chromatographic
peaks on different scales. The resulting peak list and the given peak
list (`xcmsPeaks`) are merged and redundant peaks are removed.

## Arguments

- object:

  `xcmsSet` object

- ppm:

  maxmial tolerated m/z deviation in consecutive scans, in ppm (parts
  per million)

- peakwidth:

  Chromatographic peak width, given as range (min,max) in seconds

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

- xcmsPeaks:

  peak list picked using the `centWave` algorithm with parameter
  `verbose.columns` set to TRUE (columns `scmin` and `scmax` needed)

- snthresh:

  signal to noise ratio cutoff, definition see below.

- maxcharge:

  max. number of the isotope charge.

- maxiso:

  max. number of the isotope peaks to predict for each detected feature.

- mzIntervalExtension:

  logical, if TRUE predicted isotope ROIs (regions of interest) are
  extended in the m/z dimension to increase the detection of low
  intensity and hence noisy peaks.

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
Bioinformatics 2008, 9:504\\ Hendrik Treutler and Steffen Neumann.
"Prediction, detection, and validation of isotope clusters in mass
spectrometry data" Submitted to Metabolites 2016, Special Issue
"Bioinformatics and Data Analysis"

## See also

[`findPeaks.centWave`](https://sneumann.github.io/xcms/reference/findPeaks.centWave-methods.md)
[`findPeaks-methods`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
