# Peak detection in the chromatographic time domain

The *matchedFilter* algorithm identifies peaks in the chromatographic
time domain as described in *Smith 2006*. The intensity values are
binned by cutting The LC/MS data into slices (bins) of a mass unit
(`binSize` m/z) wide. Within each bin the maximal intensity is selected.
The chromatographic peak detection is then performed in each bin by
extending it based on the `steps` parameter to generate slices
comprising bins `current_bin - steps +1` to `current_bin + steps - 1`.
Each of these slices is then filtered with matched filtration using a
second-derative Gaussian as the model peak shape. After filtration peaks
are detected using a signal-to-ratio cut-off. For more details and
illustrations see *Smith 2006*.

The `findChromPeaks,OnDiskMSnExp,MatchedFilterParam()` method performs
peak detection using the *matchedFilter* algorithm on all samples from
an
[`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
object.
[`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
objects encapsule all experiment specific data and load the spectra data
(mz and intensity values) on the fly from the original files applying
also all eventual data manipulations.

## Usage

``` r
MatchedFilterParam(
  binSize = 0.1,
  impute = "none",
  baseValue = numeric(),
  distance = numeric(),
  fwhm = 30,
  sigma = fwhm/2.3548,
  max = 5,
  snthresh = 10,
  steps = 2,
  mzdiff = 0.8 - binSize * steps,
  index = FALSE
)

# S4 method for class 'OnDiskMSnExp,MatchedFilterParam'
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

- binSize:

  `numeric(1)` specifying the width of the bins/slices in m/z dimension.

- impute:

  Character string specifying the method to be used for missing value
  imputation. Allowed values are `"none"` (no linear interpolation),
  `"lin"` (linear interpolation), `"linbase"` (linear interpolation
  within a certain bin-neighborhood) and `"intlin"`. See
  [`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  for more details.

- baseValue:

  The base value to which empty elements should be set. This is only
  considered for `method = "linbase"` and corresponds to the
  [`profBinLinBase()`](https://sneumann.github.io/xcms/reference/profGenerate.md)'s
  `baselevel` argument.

- distance:

  For `method = "linbase"`: number of non-empty neighboring element of
  an empty element that should be considered for linear interpolation.
  See details section for more information.

- fwhm:

  `numeric(1)` specifying the full width at half maximum of matched
  filtration gaussian model peak. Only used to calculate the actual
  sigma, see below.

- sigma:

  `numeric(1)` specifying the standard deviation (width) of the matched
  filtration model peak.

- max:

  `numeric(1)` representing the maximum number of peaks that are
  expected/will be identified per slice.

- snthresh:

  `numeric(1)` defining the signal to noise cutoff to be used in the
  chromatographic peak detection step.

- steps:

  `numeric(1)` defining the number of bins to be merged before
  filtration (i.e. the number of neighboring bins that will be joined to
  the slice in which filtration and peak detection will be performed).

- mzdiff:

  `numeric(1)` defining the minimum difference in m/z for peaks with
  overlapping retention times

- index:

  `logical(1)` specifying whether indicies should be returned instead of
  values for m/z and retention times.

- object:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  an `OnDiskMSnExp` object containing the MS- and all other
  experiment-relevant data.

      For all other methods: a parameter object.

- param:

  An `MatchedFilterParam` object containing all settings for the
  matchedFilter algorithm.

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

- ...:

  ignored.

## Value

The `MatchedFilterParam()` function returns a `MatchedFilterParam` class
instance with all of the settings specified for chromatographic
detection by the *matchedFilter* method.

For
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
if `return.type = "XCMSnExp"` an
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object with the results of the peak detection. If `return.type = "list"`
a list of length equal to the number of samples with matrices specifying
the identified peaks. If `return.type = "xcmsSet"` an `xcmsSet` object
with the results of the peak detection.

## Details

The intensities are binned by the provided m/z values within each
spectrum (scan). Binning is performed such that the bins are centered
around the m/z values (i.e. the first bin includes all m/z values
between `min(mz) - bin_size/2` and `min(mz) + bin_size/2`).

    For more details on binning and missing value imputation see
    [binYonX()] and [imputeLinInterpol()] methods.

Parallel processing (one process per sample) is supported and can be
configured either by the `BPPARAM` parameter or by globally defining the
parallel processing mode using the
[`BiocParallel::register()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
method from the *BiocParallel* package.

## References

Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
*Anal. Chem.* 2006, 78:779-787. doi:
[10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

## See also

The
[`do_findChromPeaks_matchedFilter()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md)
core API function and
[`findPeaks.matchedFilter()`](https://sneumann.github.io/xcms/reference/findPeaks.matchedFilter-xcmsRaw-method.md)
for the old user interface.

[`peaksWithMatchedFilter()`](https://sneumann.github.io/xcms/reference/peaksWithMatchedFilter.md)
for functions to perform matchedFilter peak detection in purely
chromatographic data.

[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for the object containing the results of the chromatographic peak
detection.

Other peak detection methods:
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md),
[`findChromPeaks-centWave`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
[`findChromPeaks-centWaveWithPredIsoROIs`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md),
[`findChromPeaks-massifquant`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
[`findPeaks-MSW`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)

## Author

Colin A Smith, Johannes Rainer

## Examples

``` r

## Create a MatchedFilterParam object. Note that we use a unnecessarily large
## binSize parameter to reduce the run-time of the example.
mfp <- MatchedFilterParam(binSize = 5, snthresh = 15)
mfp
#> Object of class:  MatchedFilterParam 
#>  Parameters:
#>  - binSize: [1] 5
#>  - impute: [1] "none"
#>  - baseValue: numeric(0)
#>  - distance: numeric(0)
#>  - fwhm: [1] 30
#>  - sigma: [1] 12.73994
#>  - max: [1] 5
#>  - snthresh: [1] 15
#>  - steps: [1] 2
#>  - mzdiff: [1] -9.2
#>  - index: [1] FALSE

## Perform the peak detection using matchecFilter on the files from the
## faahKO package. Files are read using the readMSData from the MSnbase
## package
library(faahKO)
library(MSnbase)
fls <- dir(system.file("cdf/KO", package = "faahKO"), recursive = TRUE,
           full.names = TRUE)
raw_data <- readMSData(fls[1], mode = "onDisk")
#> Polarity can not be extracted from netCDF files, please set manually the polarity with the 'polarity' method.
## Perform the chromatographic peak detection using the settings defined
## above. Note that we are also disabling parallel processing in this
## example by registering a "SerialParam"
res <- findChromPeaks(raw_data, param = mfp)
head(chromPeaks(res))
#>             mz mzmin mzmax       rt    rtmin    rtmax      into    intf  maxo
#> CP001 205.0000 205.0 205.0 2784.635 2770.550 2800.284 1778568.9 3580020 84280
#> CP002 205.0000 205.0 205.0 2784.635 2770.550 2800.284 1778568.9 3577971 84280
#> CP003 241.1460 241.1 241.2 3662.574 3646.924 3682.918 1465988.7 2234510 49728
#> CP004 241.1460 241.1 241.2 3662.574 3646.924 3682.918 1465988.7 2234510 49728
#> CP005 244.1000 244.1 244.1 2828.453 2814.369 2842.538  598990.3 1145078 31312
#> CP006 249.1591 249.1 249.2 3659.444 3643.794 3678.223 1435000.7 2367467 49040
#>            maxf i       sn sample
#> CP001 194233.12 1 63.28090      1
#> CP002 194213.46 1 66.00099      1
#> CP003  96022.23 1 25.42409      1
#> CP004  96022.23 1 25.42643      1
#> CP005  64181.64 2 16.99513      1
#> CP006 104291.09 1 36.83500      1
```
