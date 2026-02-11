# XCMSnExp filtering and subsetting

The methods listed on this page allow to filter and subset
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
objects. Most of them are inherited from the
[MSnbase::OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
object defined in the *MSnbase* package and have been adapted for
`XCMSnExp` to enable correct subsetting of preprocessing results.

- `[`: subset a `XCMSnExp` object by spectra. Be aware that this removes
  **all** preprocessing results, except adjusted retention times if
  `keepAdjustedRtime = TRUE` is passed to the method.

- `[[`: extracts a single `Spectrum` object (defined in `MSnbase`). The
  reported retention time is the adjusted retention time if alignment
  has been performed.

- `filterChromPeaks`: subset the `chromPeaks` `matrix` in `object`.
  Parameter `method` allows to specify how the chromatographic peaks
  should be filtered. Currently, only `method = "keep"` is supported
  which allows to specify chromatographic peaks to keep with parameter
  `keep` (i.e. provide a `logical`, `integer` or `character` defining
  which chromatographic peaks to keep). Feature definitions (if present)
  are updated correspondingly.

- `filterFeatureDefinitions`: allows to subset the feature definitions
  of an `XCMSnExp` object. Parameter `features` allow to define which
  features to keep. It can be a `logical`, `integer` (index of features
  to keep) or `character` (feature IDs) vector.

- `filterFile`: allows to reduce the `XCMSnExp` to data from only
  selected files. Identified chromatographic peaks for these files are
  retained while correspondence results (feature definitions) are
  removed by default. To force keeping feature definitions use
  `keepFeatures = TRUE`. Adjusted retention times (if present) are
  retained by default if present. Use `keepAdjustedRtime = FALSE` to
  drop them.

- `filterMsLevel`: reduces the `XCMSnExp` object to spectra of the
  specified MS level(s). Chromatographic peaks and identified features
  are also subsetted to the respective MS level. See also the
  `filterMsLevel` documentation in `MSnbase` for details and examples.

- `filterMz`: filters the data set based on the provided m/z value
  range. All chromatographic peaks and features (grouped peaks) with
  their apex falling within the provided mz value range are retained
  (i.e. if `chromPeaks(object)[, "mz"]` is `>= mz[1]` and `<= mz[2]`).
  Adjusted retention times, if present, are kept.

- `filterRt`: filters the data set based on the provided retention time
  range. All chromatographic peaks and features (grouped peaks) within
  the specified retention time window are retained (i.e. if the
  retention time corresponding to the peak's apex is within the
  specified rt range). If retention time correction has been performed,
  the method will by default filter the object by adjusted retention
  times. The argument `adjusted` allows to specify manually whether
  filtering should be performed on raw or adjusted retention times.
  Filtering by retention time does not drop any preprocessing results
  nor does it remove or change alignment results (i.e. adjusted
  retention times). The method returns an empty object if no spectrum or
  feature is within the specified retention time range.

- `split`: splits an `XCMSnExp` object into a `list` of `XCMSnExp`
  objects based on the provided parameter `f`. Note that by default all
  pre-processing results are removed by the splitting, except adjusted
  retention times, if the optional argument `keepAdjustedRtime = TRUE`
  is provided.

## Usage

``` r
# S4 method for class 'XCMSnExp,ANY,ANY,ANY'
x[i, j, ..., drop = TRUE]

# S4 method for class 'XCMSnExp,ANY,ANY'
x[[i, j, drop = FALSE]]

# S4 method for class 'XCMSnExp'
filterMsLevel(object, msLevel., keepAdjustedRtime = hasAdjustedRtime(object))

# S4 method for class 'XCMSnExp'
filterFile(
  object,
  file,
  keepAdjustedRtime = hasAdjustedRtime(object),
  keepFeatures = FALSE
)

# S4 method for class 'XCMSnExp'
filterMz(object, mz, msLevel., ...)

# S4 method for class 'XCMSnExp'
filterRt(object, rt, msLevel., adjusted = hasAdjustedRtime(object))

# S4 method for class 'XCMSnExp,ANY'
split(x, f, drop = FALSE, ...)

# S4 method for class 'XCMSnExp'
filterChromPeaks(
  object,
  keep = rep(TRUE, nrow(chromPeaks(object))),
  method = "keep",
  ...
)

# S4 method for class 'XCMSnExp'
filterFeatureDefinitions(object, features = integer())
```

## Arguments

- x:

  For `[` and `[[`: an `XCMSnExp` object.

- i:

  For `[`: `numeric` or `logical` vector specifying to which spectra the
  data set should be reduced. For `[[`: a single integer or character.

- j:

  For `[` and `[[`: not supported.

- ...:

  Optional additional arguments.

- drop:

  For `[` and `[[`: not supported.

- object:

  A
  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object.

- msLevel.:

  For `filterMz`, `filterRt`: `numeric` defining the MS level(s) to
  which operations should be applied or to which the object should be
  subsetted.

- keepAdjustedRtime:

  For `filterFile`, `filterMsLevel`, `[`, `split`: `logical(1)` defining
  whether the adjusted retention times should be kept, even if e.g.
  features are being removed (and the retention time correction was
  performed on these features).

- file:

  For `filterFile`: `integer` defining the file index within the object
  to subset the object by file or `character` specifying the file names
  to sub set. The indices are expected to be increasingly ordered, if
  not they are ordered internally.

- keepFeatures:

  For `filterFile`: `logical(1)` whether correspondence results (feature
  definitions) should be kept or dropped. Defaults to
  `keepFeatures = FALSE` hence feature definitions are removed from the
  returned object by default.

- mz:

  For `filterMz`: `numeric(2)` defining the lower and upper mz value for
  the filtering.

- rt:

  For `filterRt`: `numeric(2)` defining the retention time window (lower
  and upper bound) for the filtering.

- adjusted:

  For `filterRt`: `logical` indicating whether the object should be
  filtered by original (`adjusted = FALSE`) or adjusted retention times
  (`adjusted = TRUE`). For `spectra`: whether the retention times in the
  individual `Spectrum` objects should be the adjusted or raw retention
  times.

- f:

  For `split` a vector of length equal to the length of x defining how
  `x` should be splitted. It is converted internally to a `factor`.

- keep:

  For `filterChromPeaks`: `logical`, `integer` or `character` defining
  which chromatographic peaks should be retained.

- method:

  For `filterChromPeaks`: `character(1)` allowing to specify the method
  by which chromatographic peaks should be filtered. Currently only
  `method = "keep"` is supported (i.e. specify with parameter `keep`
  which chromatographic peaks should be retained).

- features:

  For `filterFeatureDefinitions`: either a `integer` specifying the
  indices of the features (rows) to keep, a `logical` with a length
  matching the number of rows of `featureDefinitions` or a `character`
  with the feature (row) names.

## Value

All methods return an
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object.

## Details

All subsetting methods try to ensure that the returned data is
consistent. Correspondence results for example are removed by default if
the data set is sub-setted by file, since the correspondence results are
dependent on the files on which correspondence was performed. This can
be changed by setting `keepFeatures = TRUE`. For adjusted retention
times, most subsetting methods support the argument `keepAdjustedRtime`
(even the `[` method) that forces the adjusted retention times to be
retained even if the default would be to drop them.

## Note

The `filterFile` method removes also process history steps not related
to the files to which the object should be sub-setted and updates the
`fileIndex` attribute accordingly. Also, the method does not allow
arbitrary ordering of the files or re-ordering of the files within the
object.

Note also that most of the filtering methods, and also the subsetting
operations `[` drop all or selected preprocessing results. To
consolidate the alignment results, i.e. ensure that adjusted retention
times are always preserved, use the
[`applyAdjustedRtime()`](https://sneumann.github.io/xcms/reference/applyAdjustedRtime.md)
function on the object that contains the alignment results. This
replaces the raw retention times with the adjusted ones.

## See also

[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for base class documentation.

[`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
for similar filter functions on `XChromatograms` objects.

## Author

Johannes Rainer

## Examples

``` r

## Loading a test data set with identified chromatographic peaks
library(MSnbase)
data(faahko_sub)
## Update the path to the files for the local system
dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")

## Disable parallel processing for this example
register(SerialParam())

## Subset the dataset to the first and third file.
xod_sub <- filterFile(faahko_sub, file = c(1, 3))

## The number of chromatographic peaks per file for the full object
table(chromPeaks(faahko_sub)[, "sample"])
#> 
#>   1   2   3 
#>  87 100  61 

## The number of chromatographic peaks per file for the subset
table(chromPeaks(xod_sub)[, "sample"])
#> 
#>  1  2 
#> 87 61 

basename(fileNames(faahko_sub))
#> [1] "ko15.CDF" "ko16.CDF" "ko18.CDF"
basename(fileNames(xod_sub))
#> [1] "ko15.CDF" "ko18.CDF"

## Filter on mz values; chromatographic peaks and features within the
## mz range are retained (as well as adjusted retention times).
xod_sub <- filterMz(faahko_sub, mz = c(300, 400))
head(chromPeaks(xod_sub))
#>        mz mzmin mzmax       rt    rtmin    rtmax       into       intb   maxo
#> CP002 302   302   302 2617.185 2595.275 2640.659   687146.6   671297.8  30552
#> CP003 344   344   344 2679.783 2646.919 2709.517  5210015.9  5135916.9 152320
#> CP005 366   366   366 2679.783 2642.224 2718.907  3365174.0  3279468.3  79928
#> CP006 343   343   343 2678.218 2637.529 2712.647 24147443.2 23703761.7 672064
#> CP007 365   365   365 2679.783 2634.399 2717.342 14975760.8 14525699.1 357632
#> CP009 313   313   313 2783.070 2768.985 2806.544  1744615.5  1698756.0  74080
#>        sn sample
#> CP002  46      1
#> CP003  68      1
#> CP005  49      1
#> CP006  87      1
#> CP007 111      1
#> CP009  44      1
nrow(chromPeaks(xod_sub))
#> [1] 89
nrow(chromPeaks(faahko_sub))
#> [1] 248

## Filter on rt values. All chromatographic peaks and features within the
## retention time range are retained. Filtering is performed by default on
## adjusted retention times, if present.
xod_sub <- filterRt(faahko_sub, rt = c(2700, 2900))

range(rtime(xod_sub))
#> [1] 2700.127 2898.877
head(chromPeaks(xod_sub))
#>          mz mzmin mzmax       rt    rtmin    rtmax       into       intb   maxo
#> CP008 579.1 579.1 579.1 2786.200 2765.855 2806.544  1669436.8  1642951.6  87656
#> CP009 313.0 313.0 313.0 2783.070 2768.985 2806.544  1744615.5  1698756.0  74080
#> CP010 301.0 301.0 301.0 2786.200 2765.855 2809.674  3051847.8  2958809.7 118512
#> CP011 279.0 279.0 279.0 2787.765 2764.290 2814.369 17140627.0 16526792.4 805248
#> CP012 279.1 279.1 279.1 2825.323 2820.629 2828.453   106769.4   106763.2  17760
#> CP018 453.2 453.2 453.2 2720.472 2678.218 2786.200  2987076.3  2986966.8  34208
#>          sn sample
#> CP008    57      1
#> CP009    44      1
#> CP010    62      1
#> CP011    45      1
#> CP012 17759      1
#> CP018 33983      1
range(chromPeaks(xod_sub)[, "rt"])
#> [1] 2701.693 2895.746

nrow(chromPeaks(faahko_sub))
#> [1] 248
nrow(chromPeaks(xod_sub))
#> [1] 21

## Extract a single Spectrum
faahko_sub[[4]]
#> Object of class "Spectrum1"
#>  Retention time: 41:46 
#>  MSn level: 1 
#>  Total ion count: 427 
#>  Polarity: -1 

## Subsetting using [ removes all preprocessing results - using
## keepAdjustedRtime = TRUE would keep adjusted retention times, if present.
xod_sub <- faahko_sub[fromFile(faahko_sub) == 1]
#> Warning: Removed preprocessing results
xod_sub
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 0.48 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 1278 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#>  MSnbase version: 2.33.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: ko15.CDF
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   ko15.CDF 
#> protocolData: none
#> featureData
#>   featureNames: F1.S0001 F1.S0002 ... F1.S1278 (1278 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (35 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> - - - xcms preprocessing - - -

## Using split does also remove preprocessing results, but it supports the
## optional parameter keepAdjustedRtime.
## Split the object into a list of XCMSnExp objects, one per file
xod_list <- split(faahko_sub, f = fromFile(faahko_sub))
xod_list
#> $`1`
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 0.48 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 1278 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#>  MSnbase version: 2.33.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: ko15.CDF
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   ko15.CDF 
#> protocolData: none
#> featureData
#>   featureNames: F1.S0001 F1.S0002 ... F1.S1278 (1278 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (35 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> - - - xcms preprocessing - - -
#> 
#> $`2`
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 0.48 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 1278 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#>  MSnbase version: 2.33.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: ko16.CDF
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   ko16.CDF 
#> protocolData: none
#> featureData
#>   featureNames: F2.S0001 F2.S0002 ... F2.S1278 (1278 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (35 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> - - - xcms preprocessing - - -
#> 
#> $`3`
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 0.48 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 1278 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#>  MSnbase version: 2.33.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: ko18.CDF
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   ko18.CDF 
#> protocolData: none
#> featureData
#>   featureNames: F3.S0001 F3.S0002 ... F3.S1278 (1278 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (35 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> - - - xcms preprocessing - - -
#> 
```
