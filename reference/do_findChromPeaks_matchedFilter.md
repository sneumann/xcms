# Core API function for matchedFilter peak detection

This function identifies peaks in the chromatographic time domain as
described in *Smith 2006*. The intensity values are binned by cutting
The LC/MS data into slices (bins) of a mass unit (`binSize` m/z) wide.
Within each bin the maximal intensity is selected. The peak detection is
then performed in each bin by extending it based on the `steps`
parameter to generate slices comprising bins `current_bin - steps +1` to
`current_bin + steps - 1`. Each of these slices is then filtered with
matched filtration using a second-derative Gaussian as the model peak
shape. After filtration peaks are detected using a signal-to-ration
cut-off. For more details and illustrations see *Smith 2006*.

## Usage

``` r
do_findChromPeaks_matchedFilter(
  mz,
  int,
  scantime,
  valsPerSpect,
  binSize = 0.1,
  impute = "none",
  baseValue,
  distance,
  fwhm = 30,
  sigma = fwhm/2.3548,
  max = 5,
  snthresh = 10,
  steps = 2,
  mzdiff = 0.8 - binSize * steps,
  index = FALSE,
  sleep = 0
)
```

## Arguments

- mz:

  Numeric vector with the individual m/z values from all scans/ spectra
  of one file/sample.

- int:

  Numeric vector with the individual intensity values from all
  scans/spectra of one file/sample.

- scantime:

  Numeric vector of length equal to the number of spectra/scans of the
  data representing the retention time of each scan.

- valsPerSpect:

  Numeric vector with the number of values for each spectrum.

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

  `numeric(1)` defining the signal to noise ratio cutoff.

- steps:

  `numeric(1)` defining the number of bins to be merged before
  filtration (i.e. the number of neighboring bins that will be joined to
  the slice in which filtration and peak detection will be performed).

- mzdiff:

  `numeric(1)` representing the minimum difference in m/z dimension
  required for peaks with overlapping retention times; can be negative
  to allow overlap. During peak post-processing, peaks defined to be
  overlapping are reduced to the one peak with the largest signal.

- index:

  `logical(1)` specifying whether indicies should be returned instead of
  values for m/z and retention times.

- sleep:

  `numeric(1)` defining the number of seconds to wait between
  iterations. Defaults to `sleep = 0`. If `> 0` a plot is generated
  visualizing the identified chromatographic peak. Note: this argument
  is for backward compatibility only and will be removed in future.

## Value

A matrix, each row representing an identified chromatographic peak, with
columns:

- `"mz"`: Intensity weighted mean of m/z values of the peak across
  scans.

- `"mzmin"`: Minimum m/z of the peak.

- `"mzmax"`: Maximum m/z of the peak.

- `"rt"`: Retention time of the peak's midpoint.

- `"rtmin"`: Minimum retention time of the peak.

- `"rtmax"`: Maximum retention time of the peak.

- `"into"`: Integrated (original) intensity of the peak.

- `"intf"`: Integrated intensity of the filtered peak.

- `"maxo"`: Maximum intensity of the peak.

- `"maxf"`: Maximum intensity of the filtered peak.

- `"i"`: Rank of peak in merged EIC (`<= max`).

- `"sn"`: Signal to noise ratio of the peak.

## Details

The intensities are binned by the provided m/z values within each
spectrum (scan). Binning is performed such that the bins are centered
around the m/z values (i.e. the first bin includes all m/z values
between `min(mz) - bin_size/2` and `min(mz) + bin_size/2`).

    For more details on binning and missing value imputation see
    [binYonX()] and [imputeLinInterpol()] functions.

## Note

This function exposes core peak detection functionality of the
*matchedFilter* method.

## References

Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
*Anal. Chem.* 2006, 78:779-787. doi:
[10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

## See also

[`binYonX()`](https://sneumann.github.io/xcms/reference/binYonX.md) for
a binning function,
[`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
for the interpolation of missing values.

Other core peak detection functions:
[`do_findChromPeaks_centWave()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWave.md),
[`do_findChromPeaks_centWaveWithPredIsoROIs()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_centWaveWithPredIsoROIs.md),
[`do_findChromPeaks_massifquant()`](https://sneumann.github.io/xcms/reference/do_findChromPeaks_massifquant.md),
[`do_findPeaks_MSW()`](https://sneumann.github.io/xcms/reference/do_findPeaks_MSW.md)

## Author

Colin A Smith, Johannes Rainer

## Examples

``` r

## Load the test file
faahko_sub <- loadXcmsData("faahko_sub")

## Subset to one file and restrict to a certain retention time range
data <- filterRt(filterFile(faahko_sub, 1), c(2500, 3000))

## Get m/z and intensity values
mzs <- mz(data)
ints <- intensity(data)

## Define the values per spectrum:
valsPerSpect <- lengths(mzs)

res <- do_findChromPeaks_matchedFilter(mz = unlist(mzs), int = unlist(ints),
    scantime = rtime(data), valsPerSpect = valsPerSpect)
head(res)
#>            mz mzmin mzmax       rt    rtmin    rtmax      into      intf  maxo
#> [1,] 205.0000 205.0 205.0 2784.635 2770.550 2800.284 1778568.9 3610062.2 84280
#> [2,] 205.9819 205.9 206.0 2786.200 2772.115 2800.284  237993.6  448580.3 10681
#> [3,] 207.0821 207.0 207.1 2712.647 2698.562 2726.731  380873.0  730981.4 18800
#> [4,] 236.0956 236.0 236.1 2518.593 2504.508 2534.242  252282.0  458747.7 12957
#> [5,] 244.1000 244.1 244.1 2828.453 2814.369 2844.103  612169.9 1279308.9 31312
#> [6,] 266.0751 266.0 266.1 2828.453 2815.934 2844.103  113219.0  214886.8  5801
#>           maxf i       sn
#> [1,] 195026.48 1 28.33394
#> [2,]  23860.11 1 16.53987
#> [3,]  40065.74 1 12.87314
#> [4,]  24536.55 1 14.99012
#> [5,]  69898.24 1 24.20989
#> [6,]  11773.56 1 10.83870
```
