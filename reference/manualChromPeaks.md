# Manual peak integration and feature definition

The `manualChromPeaks` function allows to *manually* define
chromatographic peaks, integrate the intensities within the specified
peak area and add them to the object's `chromPeaks` matrix. A peak is
not added for a sample if no signal was found in the respective data
file.

Because chromatographic peaks are added to eventually previously
identified peaks, it is suggested to run
[`refineChromPeaks()`](https://sneumann.github.io/xcms/reference/refineChromPeaks.md)
with the
[`MergeNeighboringPeaksParam()`](https://sneumann.github.io/xcms/reference/refineChromPeaks.md)
approach to merge potentially overlapping peaks.

The `manualFeatures` function allows to manually group identified
chromatographic peaks into features by providing their index in the
object's `chromPeaks` matrix.

## Usage

``` r
manualChromPeaks(object, ...)

manualFeatures(object, ...)

# S4 method for class 'MsExperiment'
manualChromPeaks(
  object,
  chromPeaks = matrix(numeric()),
  samples = seq_along(object),
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperiment'
manualChromPeaks(
  object,
  chromPeaks = matrix(numeric()),
  samples = seq_along(object),
  msLevel = 1L,
  chunkSize = 2L,
  BPPARAM = bpparam()
)

# S4 method for class 'XcmsExperiment'
manualFeatures(object, peakIdx = list(), msLevel = 1L)

# S4 method for class 'OnDiskMSnExp'
manualChromPeaks(
  object,
  chromPeaks = matrix(),
  samples = seq_along(fileNames(object)),
  msLevel = 1L,
  BPPARAM = bpparam()
)

# S4 method for class 'XCMSnExp'
manualChromPeaks(
  object,
  chromPeaks = matrix(),
  samples = seq_along(fileNames(object)),
  msLevel = 1L,
  BPPARAM = bpparam()
)

# S4 method for class 'XCMSnExp'
manualFeatures(object, peakIdx = list(), msLevel = 1L)
```

## Arguments

- object:

  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md),
  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  or
  [MSnbase::OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
  object.

- ...:

  ignored.

- chromPeaks:

  For `manualChromPeaks`: `matrix` defining the boundaries of the
  chromatographic peaks with one row per chromatographic peak and
  columns `"mzmin"`, `"mzmax"`, `"rtmin"` and `"rtmax"` defining the m/z
  and retention time region of each peak.

- samples:

  For `manualChromPeaks`: optional `integer` defining individual samples
  in which the peak integration should be performed. Defaults to all
  samples.

- msLevel:

  `integer(1)` defining the MS level in which peak integration should be
  performed. Only a single MS level at a time is supported. Defaults to
  `msLevel = 1L`.

- chunkSize:

  `integer(1)` defining the number of files (samples) that should be
  loaded into memory and processed at the same time. Peak integration is
  then performed in parallel (per sample) on this subset data. This
  setting thus allows to balance between memory demand and speed (due to
  parallel processing). Because parallel processing can only performed
  on the subset of data currently loaded into memory in each iteration,
  the value for `chunkSize` should match the defined parallel setting
  setup. Using a parallel processing setup using 4 CPUs (separate
  processes) but using
  `chunkSize = `1`will not perform any parallel processing, as only the data from one sample is loaded in memory at a time. On the other hand, setting`chunkSize\`
  to the total number of samples in an experiment will load the full MS
  data into memory and will thus in most settings cause an out-of-memory
  error.

- BPPARAM:

  parallel processing settings (see
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details).

- peakIdx:

  For `manualFeatures`: `list` of `integer` vectors with the indices of
  chromatographic peaks in the object's `chromPeaks` matrix that should
  be grouped into features.

## Value

`XcmsExperiment` or `XCMSnExp` with the manually added chromatographic
peaks or features.

## Author

Johannes Rainer

## Examples

``` r

## Read a test dataset.
library(MsDataHub)
fls <- MsDataHub::PestMix1_DDA.mzML()
#> see ?MsDataHub and browseVignettes('MsDataHub') for documentation
#> loading from cache

## Define a data frame with some sample annotations
ann <- data.frame(
    injection_index = 1,
    sample_id = c("Pest_mix"))

## Import the data
library(MsExperiment)
mse <- readMsExperiment(fls)

## Define some arbitrary peak areas
pks <- cbind(
    mzmin = c(512, 234.3), mzmax = c(513, 235),
    rtmin = c(10, 33), rtmax = c(19, 50)
)
pks
#>      mzmin mzmax rtmin rtmax
#> [1,] 512.0   513    10    19
#> [2,] 234.3   235    33    50

res <- manualChromPeaks(mse, pks)
chromPeaks(res)
#>           mz mzmin mzmax     rt rtmin rtmax      into sample      maxo sn
#> CP1 234.9081 234.3   235 45.502    33    50 0.7627435      1 0.4139576 NA

## Peaks were only found in the second file.
```
