# Single-spectrum non-chromatography MS data peak detection

Perform peak detection in mass spectrometry direct injection spectrum
using a wavelet based algorithm.

The `findChromPeaks,OnDiskMSnExp,MSWParam()` method performs peak
detection in single-spectrum non-chromatography MS data using
functionality from the *MassSpecWavelet* package on all samples from an
`OnDiskMSnExp` object. `OnDiskMSnExp` objects encapsule all experiment
specific data and load the spectra data (mz and intensity values) on the
fly from the original files applying also all eventual data
manipulations.

## Usage

``` r
MSWParam(
  snthresh = 3,
  verboseColumns = FALSE,
  scales = c(1, seq(2, 30, 2), seq(32, 64, 4)),
  nearbyPeak = TRUE,
  peakScaleRange = 5,
  ampTh = 0.01,
  minNoiseLevel = ampTh/snthresh,
  ridgeLength = 24,
  peakThr = NULL,
  tuneIn = FALSE,
  ...
)

# S4 method for class 'OnDiskMSnExp,MSWParam'
findChromPeaks(
  object,
  param,
  BPPARAM = bpparam(),
  return.type = "XCMSnExp",
  msLevel = 1L,
  ...
)
```

## Arguments

- snthresh:

  `numeric(1)` defining the signal to noise ratio cutoff.

- verboseColumns:

  `logical(1)` whether additional peak meta data columns should be
  returned.

- scales:

  Numeric defining the scales of the continuous wavelet transform (CWT).

- nearbyPeak:

  `logical(1)` whether to include nearby peaks of major peaks.

- peakScaleRange:

  `numeric(1)` defining the scale range of the peak (larger than 5 by
  default).

- ampTh:

  `numeric(1)` defining the minimum required relative amplitude of the
  peak (ratio of the maximum of CWT coefficients).

- minNoiseLevel:

  `numeric(1)` defining the minimum noise level used in computing the
  SNR.

- ridgeLength:

  `numeric(1)` defining the minimum highest scale of the peak in 2-D CWT
  coefficient matrix.

- peakThr:

  `numeric(1)` with the minimum absolute intensity (above baseline) of
  peaks to be picked. If provided, the smoothing Savitzky-Golay filter
  is used (in the *MassSpecWavelet*) package to estimate the local
  intensity.

- tuneIn:

  `logical(1)` whther to tune in the parameter estimation of the
  detected peaks.

- ...:

  Additional parameters to be passed to the `peakDetectionCWT()` and
  `identifyMajorPeaks()` functions from the *MassSpecWavelet* package.

- object:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  an `OnDiskMSnExp` object containing the MS- and all other
  experiment-relevant data.

      For all other methods: a parameter object.

- param:

  An `MSWParam` object containing all settings for the algorithm.

- BPPARAM:

  A parameter class specifying if and how parallel processing should be
  performed. It defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).
  See documentation of the *BiocParallel* package for more details. If
  parallel processing is enabled, peak detection is performed in
  parallel on several of the input samples.

- return.type:

  Character specifying what type of object the method should return. Can
  be either `"XCMSnExp"` (default), `"list"` or `"xcmsSet"`.

- msLevel:

  `integer(1)` defining the MS level on which the peak detection should
  be performed. Defaults to `msLevel = 1`.

## Value

The `MSWParam()` function returns a `MSWParam` class instance with all
of the settings specified for peak detection by the *MSW* method.

For
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
if `return.type = "XCMSnExp"` an `XCMSnExp` object with the results of
the peak detection. If `return.type = "list"` a list of length equal to
the number of samples with matrices specifying the identified peaks. If
`return.type = "xcmsSet"` an `xcmsSet` object with the results of the
detection.

## Details

This is a wrapper for the peak picker in Bioconductor's
*MassSpecWavelet* package calling `peakDetectionCWT` and
`tuneInPeakInfo` functions. See the *xcmsDirect* vignette for more
information.

Parallel processing (one process per sample) is supported and can be
configured either by the `BPPARAM` parameter or by globally defining the
parallel processing mode using the
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
method from the *BiocParallel* package.

## See also

The
[`do_findPeaks_MSW()`](https://sneumann.github.io/xcms/reference/do_findPeaks_MSW.md)
core API function and
[`findPeaks.MSW()`](https://sneumann.github.io/xcms/reference/findPeaks.MSW-xcmsRaw-method.md)
for the old user interface.

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for the object containing the results of the peak detection.

Other peak detection methods:
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md),
[`findChromPeaks-centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
[`findChromPeaks-centWaveWithPredIsoROIs`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md),
[`findChromPeaks-massifquant`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
[`findChromPeaks-matchedFilter`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)

## Author

Joachim Kutzera, Steffen Neumann, Johannes Rainer

## Examples

``` r

library(MSnbase)
## Create a MSWParam object
mp <- MSWParam(snthresh = 15)
mp
#> Object of class:  MSWParam 
#>  Parameters:
#>  - snthresh: [1] 15
#>  - verboseColumns: [1] FALSE
#>  - scales:  [1]  1  2  4  6  8 10 12 14 16 18 20 22 24 26 28 30 32 36 40 44 48 52 56 60 64
#>  - nearbyPeak: [1] TRUE
#>  - peakScaleRange: [1] 5
#>  - amp.Th: [1] 0.01
#>  - minNoiseLevel: [1] 0.0006666667
#>  - ridgeLength: [1] 24
#>  - peakThr: numeric(0)
#>  - tuneIn: [1] FALSE

## Loading a small subset of direct injection, single spectrum files
library(msdata)
#> 
#>  IMPORTANT: msdata will be deprecated and replaced by the new
#>  MsDataHub package -- see https://bioconductor.org/packages/MsDataHub.
#>  Please open an issue at https://github.com/rformassspectrometry/MsDataHub/issues
#>  to inform us of any data from msdata that would need to be moved to MsDataHub.
fticrf <- list.files(system.file("fticr-mzML", package = "msdata"),
                    recursive = TRUE, full.names = TRUE)
fticr <- readMSData(fticrf[1], msLevel. = 1, mode = "onDisk")

## Perform the MSW peak detection on these:
p <- MSWParam(scales = c(1, 7), peakThr = 80000, ampTh = 0.005,
             SNR.method = "data.mean", winSize.noise = 500)
fticr <- findChromPeaks(fticr, param = p)

head(chromPeaks(fticr))
#>            mz    mzmin    mzmax rt rtmin rtmax    into     maxo        sn intf
#> CP01 403.2367 403.2279 403.2447 -1    -1    -1 4735258 372259.4 22.975342   NA
#> CP02 409.1845 409.1755 409.1936 -1    -1    -1 4158404 310572.1 20.597978   NA
#> CP03 413.2677 413.2585 413.2769 -1    -1    -1 6099006 435462.6 27.217229   NA
#> CP04 414.2719 414.2652 414.2786 -1    -1    -1 1916658 157613.1  8.777433   NA
#> CP05 423.2363 423.2266 423.2459 -1    -1    -1 2708391 174252.7 14.745275   NA
#> CP06 427.2681 427.2583 427.2779 -1    -1    -1 6302089 461385.6 32.468892   NA
#>           maxf sample
#> CP01  814693.1      1
#> CP02  731557.2      1
#> CP03 1018994.8      1
#> CP04  280826.9      1
#> CP05  435858.5      1
#> CP06 1124549.4      1
```
