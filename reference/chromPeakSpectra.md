# Extract spectra associated with chromatographic peaks

Extract (MS1 or MS2) spectra from an
[XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
or
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object for identified chromatographic peaks. To return spectra for
selected chromatographic peaks, their *peak ID* (i.e., row name in the
`chromPeaks` matrix) can be provided with parameter `peaks`. For
`msLevel = 1L` (only supported for `return.type = "Spectra"` or
`return.type = "List"`) MS1 spectra within the retention time boundaries
(in the file in which the peak was detected) are returned. For
`msLevel = 2L` MS2 spectra are returned for a chromatographic peak if
their precursor m/z is within the retention time and m/z range of the
chromatographic peak. Parameter `method` allows to define whether all or
a single spectrum should be returned:

- `method = "all"`: (default): return all spectra for each
  chromatographic peak.

- `method = "closest_rt"`: return the spectrum with the retention time
  closest to the peak's retention time (at apex).

- `method = "closest_mz"`: return the spectrum with the precursor m/z
  closest to the peaks's m/z (at apex); only supported for
  `msLevel > 1`.

- `method = "largest_tic"`: return the spectrum with the largest total
  signal (sum of peaks intensities).

- `method = "largest_bpi"`: return the spectrum with the largest peak
  intensity (maximal peak intensity).

- `method = "signal"`: only for `object` being a `XCMSnExp`: return the
  spectrum with the sum of intensities most similar to the peak's apex
  signal (`"maxo"`); only supported for `msLevel = 2L`.

Parameter `return.type` allows to specify the *type* of the result
object. With `return.type = "Spectra"` (the default) a
[Spectra::Spectra](https://rdrr.io/pkg/Spectra/man/Spectra.html) object
with all matching spectra is returned. With `return.type = "Spectra"` a
`List` of `Spectra` is returned. The length of the list is equal to the
number of rows of `chromPeaks`. Each element of the list contains thus a
`Spectra` with all spectra for one chromatographic peak (or a `Spectra`
of length 0 if no spectrum was found for the respective chromatographic
peak).

Parameter `chromPeakColumns` allows the user to add specific metadata
columns from the chromatographic peaks (`chromPeaks`) to the returned
spectra object. This can be useful to keep information such as retention
time (`rt`), m/z (`mz`). The columns will be named as they are written
in the `chromPeaks` object with the prefix `"chrom_peak_"`. The *peak
ID* (i.e., the row name of the peak in the `chromPeaks` matrix) is
always added to the spectra object as a metadata column named
`"chrom_peak_id"`.

See also the *LC-MS/MS data analysis* vignette for more details and
examples.

## Usage

``` r
chromPeakSpectra(object, ...)

# S4 method for class 'XcmsExperiment'
chromPeakSpectra(
  object,
  method = c("all", "closest_rt", "closest_mz", "largest_tic", "largest_bpi"),
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  skipFilled = FALSE,
  peaks = character(),
  chromPeakColumns = c("rt", "mz"),
  return.type = c("Spectra", "List"),
  BPPARAM = bpparam()
)

# S4 method for class 'XCMSnExp'
chromPeakSpectra(
  object,
  msLevel = 2L,
  expandRt = 0,
  expandMz = 0,
  ppm = 0,
  method = c("all", "closest_rt", "closest_mz", "signal", "largest_tic", "largest_bpi"),
  skipFilled = FALSE,
  return.type = c("Spectra", "MSpectra", "List", "list"),
  peaks = character()
)
```

## Arguments

- object:

  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  or
  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with identified chromatographic peaks for which spectra should
  be returned.

- ...:

  ignored.

- method:

  `character(1)` specifying which spectra to include in the result.
  Defaults to `method = "all"`. See function description for details.

- msLevel:

  `integer(1)` defining the MS level of the spectra that should be
  returned.

- expandRt:

  `numeric(1)` to expand the retention time range of each peak by a
  constant value on each side.

- expandMz:

  `numeric(1)` to expand the m/z range of each peak by a constant value
  on each side.

- ppm:

  `numeric(1)` to expand the m/z range of each peak (on each side) by a
  value dependent on the peak's m/z.

- skipFilled:

  `logical(1)` whether spectra for filled-in peaks should be reported or
  not. Defaults to `skipFilled = FALSE` thus also spectra for gap-filled
  chromatographic peaks are returned. Set to `skipFilled = TRUE` to get
  only spectra for detected chromatographic peaks.

- peaks:

  `character`, `logical` or `integer` allowing to specify a subset of
  chromatographic peaks in `chromPeaks` for which spectra should be
  returned (providing either their ID, a logical vector same length than
  `nrow(chromPeaks(x))` or their index in `chromPeaks(x)`). Be aware
  that when `peaks` are provided, parameter `skipFilled` is ignored.
  Spectra are returned for any chromatographic peak, detected or
  gap-filled, that are defined with `peaks`.

- chromPeakColumns:

  `character` vector with the names of the columns from `chromPeaks`
  that should be added to the returned spectra object. The columns will
  be named as they are written in the `chromPeaks` object with a prefix
  `"chrom_peak_"`. Defaults to `c("mz", "rt")`.

- return.type:

  `character(1)` defining the type of result object that should be
  returned.

- BPPARAM:

  parallel processing setup. Defaults to
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).

## Value

parameter `return.type` allow to specify the type of the returned
object:

- `return.type = "Spectra"` (default): a `Spectra` object (defined in
  the `Spectra` package). The result contains all spectra for all peaks.
  Metadata column `"peak_id"` provides the ID of the respective peak
  (i.e. its rowname in
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md).

- `return.type = "List"`: `List` of length equal to the number of
  chromatographic peaks is returned, each element being a `Spectra` with
  the spectra for one chromatographic peak.

For backward compatibility options `"MSpectra"` and `"list"` are also
supported but are not suggested.

- `return.type = "MSpectra"` (deprecated): a
  [MSnbase::MSpectra](https://lgatto.github.io/MSnbase/reference/MSpectra.html)
  object with elements being `Spectrum` objects. The result objects
  contains all spectra for all peaks. Metadata column `"peak_id"`
  provides the ID of the respective peak (i.e. its rowname in
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)).

- `return.type = "list"`: `list` of `list`s that are either of length 0
  or contain `Spectrum2` object(s) within the m/z-rt range. The length
  of the list matches the number of peaks.

## Author

Johannes Rainer

## Examples

``` r

## Read a file with DDA LC-MS/MS data
library(MsExperiment)
fl <- system.file("TripleTOF-SWATH/PestMix1_DDA.mzML", package = "msdata")

dda <- readMsExperiment(fl)

## Perform MS1 peak detection
dda <- findChromPeaks(dda, CentWaveParam(peakwidth = c(5, 15),
    prefilter = c(5, 1000)))

## Return all MS2 spectro for each chromatographic peaks as a Spectra object
ms2_sps <- chromPeakSpectra(dda)
ms2_sps
#> MSn data (Spectra) with 93 spectra in a MsBackendMzR backend:
#>       msLevel     rtime scanIndex
#>     <integer> <numeric> <integer>
#> 1           2   237.460      1810
#> 2           2   237.869      1812
#> 3           2   240.759      1842
#> 4           2   241.299      1846
#> 5           2   318.444      2391
#> ...       ...       ...       ...
#> 89          2   552.147      4948
#> 90          2   570.506      5077
#> 91          2   574.435      5108
#> 92          2   574.725      5110
#> 93          2   575.255      5115
#>  ... 37 more variables/columns.
#> 
#> file(s):
#> PestMix1_DDA.mzML
#> Processing:
#>  Filter: select MS level(s) 2 [Wed Feb 11 13:33:41 2026]
#>  Merge 1 Spectra into one [Wed Feb 11 13:33:41 2026] 

## spectra variable *chrom_peak_id* contain the row names of the peaks in the
## chromPeak matrix and allow thus to map chromatographic peaks to the
## returned MS2 spectra
ms2_sps$chrom_peak_id
#>  [1] "CP01" "CP01" "CP01" "CP01" "CP02" "CP03" "CP03" "CP03" "CP03" "CP04"
#> [11] "CP04" "CP04" "CP04" "CP06" "CP06" "CP06" "CP06" "CP06" "CP06" "CP07"
#> [21] "CP07" "CP07" "CP07" "CP07" "CP07" "CP07" "CP08" "CP08" "CP09" "CP09"
#> [31] "CP10" "CP10" "CP11" "CP11" "CP11" "CP11" "CP11" "CP12" "CP12" "CP12"
#> [41] "CP12" "CP12" "CP12" "CP13" "CP13" "CP13" "CP13" "CP13" "CP14" "CP14"
#> [51] "CP14" "CP15" "CP15" "CP15" "CP15" "CP15" "CP15" "CP17" "CP17" "CP17"
#> [61] "CP17" "CP17" "CP18" "CP18" "CP18" "CP18" "CP18" "CP19" "CP19" "CP19"
#> [71] "CP19" "CP19" "CP19" "CP20" "CP20" "CP20" "CP21" "CP21" "CP21" "CP22"
#> [81] "CP22" "CP22" "CP23" "CP23" "CP23" "CP23" "CP25" "CP25" "CP25" "CP27"
#> [91] "CP27" "CP27" "CP27"
chromPeaks(dda)
#>            mz    mzmin    mzmax      rt   rtmin   rtmax        into        intb
#> CP01 219.0961 219.0933 219.0996 240.897 233.595 247.874 15236.49207 15222.54500
#> CP02 219.0947 219.0926 219.0963 319.236 311.979 324.874    76.97947    63.82131
#> CP03 227.9913 227.9902 227.9921 358.349 351.649 367.218  5511.96730  5496.32244
#> CP04 300.0312 300.0289 300.0321 366.568 360.029 375.738 12659.45398 12643.63044
#> CP05 302.0284 302.0273 302.0295 366.568 360.449 375.738  7953.70018  7938.19680
#> CP06 256.1102 256.1095 256.1109 388.167 379.097 397.326  7879.79234  7861.19876
#> CP07 224.0836 224.0831 224.0843 388.587 379.517 397.446  8904.69865  8886.41107
#> CP08 309.9980 309.9962 309.9987 405.715 399.527 415.525  5381.59741  5365.38152
#> CP09 378.0356 378.0341 378.0364 405.715 399.527 415.525  5117.90702  5101.66617
#> CP10 298.2751 298.2739 298.2781 345.711 339.830 354.979 39574.83769 39559.51278
#> CP11 308.0012 307.9983 308.0018 405.715 399.527 415.525  5630.65430  5614.25375
#> CP12 376.0385 376.0359 376.0391 405.715 399.527 415.525  5198.61440  5182.23188
#> CP13 292.1216 292.1204 292.1223 410.555 403.625 420.035 15783.11419 15766.32540
#> CP14 289.1217 289.1202 289.1229 425.434 415.525 433.478 10049.05516 10031.22736
#> CP15 304.1141 304.1119 304.1192 424.614 417.985 434.313 14704.14891 14687.51150
#> CP16 384.9700 384.9684 384.9706 440.533 433.893 449.312  8075.24022  8059.47631
#> CP17 382.9727 382.9715 382.9739 440.533 433.478 449.312  8307.76833  8291.81980
#> CP18 231.0144 231.0134 231.0151 450.992 444.342 460.691  5612.28545  5595.90874
#> CP19 306.1659 306.1624 306.1694 463.786 453.922 471.431 31273.95023 31256.40745
#> CP20 205.0970 205.0956 205.0976 459.421 453.092 466.461  5609.65542  5596.09130
#> CP21 230.9875 230.9862 230.9879 495.109 486.210 502.807  8036.59181  8019.64904
#> CP22 216.1414 216.1398 216.1418 508.678 502.680 516.408  4406.91863  4392.98419
#> CP23 330.2065 330.2008 330.2078 422.514 416.754 431.858   614.91821   599.73617
#> CP24 330.2059 330.2036 330.2081 474.408 460.261 479.107   248.31997   229.20853
#> CP25 330.2067 330.2054 330.2077 552.815 541.355 561.837  5890.30116  5869.52232
#> CP26 330.2057 330.2038 330.2076 505.738 495.519 515.868   315.87230   295.17246
#> CP27 373.0416 373.0406 373.0430 576.468 569.174 580.636  6184.87790  6173.40597
#>            maxo   sn sample
#> CP01 4877.11621 4876      1
#> CP02   19.35541   18      1
#> CP03 1459.38428 1458      1
#> CP04 3156.89258 3156      1
#> CP05 1974.66846 1974      1
#> CP06 1836.12695 1835      1
#> CP07 2268.62915 2268      1
#> CP08 1345.32336 1344      1
#> CP09 1300.64075 1300      1
#> CP10 8060.13037 8059      1
#> CP11 1450.15808 1449      1
#> CP12 1350.63330 1350      1
#> CP13 3203.84229 3203      1
#> CP14 2972.63916 2972      1
#> CP15 5009.70654 5009      1
#> CP16 2224.37622 2223      1
#> CP17 2361.33594 2360      1
#> CP18 1671.33203 1670      1
#> CP19 6300.43652 6299      1
#> CP20 1717.72607 1717      1
#> CP21 2244.02979 2243      1
#> CP22 1576.85681 1576      1
#> CP23  135.81871  135      1
#> CP24   19.64630   19      1
#> CP25 2102.92871 2102      1
#> CP26   27.88095   27      1
#> CP27 2511.13501 2510      1

## Alternatively, return the result as a List of Spectra objects. This list
## is parallel to chromPeaks hence the mapping between chromatographic peaks
## and MS2 spectra is easier.
ms2_sps <- chromPeakSpectra(dda, return.type = "List")
names(ms2_sps)
#>  [1] "CP01" "CP02" "CP03" "CP04" NA     "CP06" "CP07" "CP08" "CP09" "CP10"
#> [11] "CP11" "CP12" "CP13" "CP14" "CP15" NA     "CP17" "CP18" "CP19" "CP20"
#> [21] "CP21" "CP22" "CP23" NA     "CP25" NA     "CP27"
rownames(chromPeaks(dda))
#>  [1] "CP01" "CP02" "CP03" "CP04" "CP05" "CP06" "CP07" "CP08" "CP09" "CP10"
#> [11] "CP11" "CP12" "CP13" "CP14" "CP15" "CP16" "CP17" "CP18" "CP19" "CP20"
#> [21] "CP21" "CP22" "CP23" "CP24" "CP25" "CP26" "CP27"
ms2_sps[[1L]]
#> MSn data (Spectra) with 4 spectra in a MsBackendMzR backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         2   237.460      1810
#> 2         2   237.869      1812
#> 3         2   240.759      1842
#> 4         2   241.299      1846
#>  ... 37 more variables/columns.
#> 
#> file(s):
#> PestMix1_DDA.mzML
#> Processing:
#>  Filter: select MS level(s) 2 [Wed Feb 11 13:33:41 2026]
#>  Merge 1 Spectra into one [Wed Feb 11 13:33:41 2026] 

## Parameter `msLevel` allows to define from which MS level spectra should
## be returned. By default `msLevel = 2L` but with `msLevel = 1L` all
## MS1 spectra with a retention time within the retention time range of
## a chromatographic peak can be returned. Alternatively, selected
## spectra can be returned by specifying the selection criteria/method
## with the `method` parameter. Below we extract for each chromatographic
## peak the MS1 spectra with a retention time closest to the
## chromatographic peak's apex position. Alternatively it would also be
## possible to select the spectrum with the highest total signal or
## highest (maximal) intensity.
ms1_sps <- chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt")
ms1_sps
#> MSn data (Spectra) with 27 spectra in a MsBackendMzR backend:
#>       msLevel     rtime scanIndex
#>     <integer> <numeric> <integer>
#> 1           1   240.897      1843
#> 2           1   319.236      2396
#> 3           1   358.349      2813
#> 4           1   366.568      2924
#> 5           1   366.568      2924
#> ...       ...       ...       ...
#> 23          1   422.514      3571
#> 24          1   474.408      4266
#> 25          1   552.815      4958
#> 26          1   505.738      4549
#> 27          1   576.468      5122
#>  ... 37 more variables/columns.
#> 
#> file(s):
#> PestMix1_DDA.mzML
#> Processing:
#>  Filter: select MS level(s) 1 [Wed Feb 11 13:33:41 2026]
#>  Merge 1 Spectra into one [Wed Feb 11 13:33:41 2026] 

## Parameter peaks would allow to extract spectra for specific peaks only.
## Peaks can be defined with parameter `peaks` which can be either an
## `integer` with the index of the peak in the `chromPeaks` matrix or a
## `character` with its rowname in `chromPeaks`.
chromPeakSpectra(dda, msLevel = 1L, method = "closest_rt", peaks = c(3, 5))
#> MSn data (Spectra) with 2 spectra in a MsBackendMzR backend:
#>     msLevel     rtime scanIndex
#>   <integer> <numeric> <integer>
#> 1         1   358.349      2813
#> 2         1   366.568      2924
#>  ... 37 more variables/columns.
#> 
#> file(s):
#> PestMix1_DDA.mzML
#> Processing:
#>  Filter: select MS level(s) 1 [Wed Feb 11 13:33:41 2026]
#>  Merge 1 Spectra into one [Wed Feb 11 13:33:41 2026] 
```
