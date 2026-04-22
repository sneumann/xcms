# Data container storing xcms preprocessing results

The `XCMSnExp` object is a container for the results of a G/LC-MS data
preprocessing that comprises chromatographic peak detection, alignment
and correspondence. These results can be accessed with the
`chromPeaks()`, `adjustedRtime()` and `featureDefinitions()` functions;
see below (after the Usage, Arguments, Value and Slots sections) for
more details). Along with the results, the object contains the
processing history that allows to track each processing step along with
the used settings. This can be extracted with the `processHistory()`
function. `XCMSnExp` objects, by directly extending the
[MSnbase::OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
object from the *MSnbase* package, inherit all of its functionality and
allows thus an easy access to the full raw data at any stage of an
analysis. To support interaction with packages requiring the *old*
objects, `XCMSnExp` objects can be coerced into `xcmsSet` objects using
the `as()` method (see examples below). All preprocessing results will
be passed along to the resulting object.

General functions for `XCMSnExp` objects are (see further below for
specific function to handle chromatographic peak data, alignment and
correspondence results):

`processHistoryTypes()` returns the available *types* of process
histories. These can be passed with argument `type` to the
`processHistory` method to extract specific process step(s).

`hasFilledChromPeaks()`: whether filled-in peaks are present or not.

[`profMat()`](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md):
creates a *profile matrix*, which is a n x m matrix, n (rows)
representing equally spaced m/z values (bins) and m (columns) the
retention time of the corresponding scans. Each cell contains the
maximum intensity measured for the specific scan and m/z values. See
[`profMat()`](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md)
for more details and description of the various binning methods.

`hasAdjustedRtime()`: whether the object provides adjusted retention
times.

`hasFeatures()`: whether the object contains correspondence results
(i.e. features).

`hasChromPeaks()`: whether the object contains peak detection results.

`hasFilledChromPeaks()`: whether the object contains any filled-in
chromatographic peaks.

`adjustedRtime()`,`adjustedRtime<-`: extract/set adjusted retention
times. `adjustedRtime<-` should not be called manually, it is called
internally by the
[`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
methods. For `XCMSnExp` objects, `adjustedRtime<-` does also apply
retention time adjustments to eventually present chromatographic peaks.
The `bySample` parameter allows to specify whether the adjusted
retention time should be grouped by sample (file).

`featureDefinitions()`, `featureDefinitions<-`: extract or set the
correspondence results, i.e. the mz-rt features (peak groups). Similar
to the `chromPeaks()` it is possible to extract features for specified
m/z and/or rt ranges. The function supports also the parameter `type`
that allows to specify which features to be returned if any of `rt` or
`mz` is specified. For details see help of `chromPeaks()`. See also
[`featureSummary()`](https://sneumann.github.io/xcms/reference/featureSummary.md)
for a function to calculate simple feature summaries.

`chromPeaks()`, `chromPeaks<-`: extract or set the matrix containing the
information on identified chromatographic peaks. Rownames of the matrix
represent unique IDs of the respective peaks within the experiment.
Parameter `bySample` allows to specify whether peaks should be returned
ungrouped (default `bySample = FALSE`) or grouped by sample
(`bySample = TRUE`). The `chromPeaks<-` method for `XCMSnExp` objects
removes also all correspondence (peak grouping) and retention time
correction (alignment) results. The optional arguments `rt`, `mz`, `ppm`
and `type` allow to extract only chromatographic peaks overlapping the
defined retention time and/or m/z ranges. Argument `type` allows to
define how *overlapping* is determined: for `type == "any"` (the
default), all peaks that are even partially overlapping the region are
returned (i.e. for which either `"mzmin"` or `"mzmax"` of the
`chromPeaks` or `featureDefinitions` matrix are within the provided m/z
range), for `type == "within"` the full peak has to be within the region
(i.e. both `"mzmin"` and `"mzmax"` have to be within the m/z range) and
for `type == "apex_within"` the peak's apex position (highest signal of
the peak) has to be within the region (i.e. the peak's or features m/z
has to be within the m/z range). See description of the return value for
details on the returned matrix. Users usually don't have to use the
`chromPeaks<-` method directly as detected chromatographic peaks are
added to the object by the
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
method. Also, `chromPeaks<-` will replace any existing `chromPeakData`.

`chromPeakData()` and `chromPeakData<-` allow to get or set arbitrary
chromatographic peak annotations. These are returned or ar returned as a
`DataFrame`. Note that the number of rows and the rownames of the
`DataFrame` have to match those of `chromPeaks`. Parameter `columns`
allows to extract only selected columns from the `chromPeakData`. By
default (`columns = character()`) all columns are returned.

[`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
extracts the retention time for each scan. The `bySample` parameter
allows to return the values grouped by sample/file and `adjusted`
whether adjusted or raw retention times should be returned. By default
the method returns adjusted retention times, if they are available (i.e.
if retention times were adjusted using the
[`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
method).

[`mz()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
extracts the mz values from each scan of all files within an `XCMSnExp`
object. These values are extracted from the original data files and
eventual processing steps are applied *on the fly*. Using the `bySample`
parameter it is possible to switch from the default grouping of mz
values by spectrum/scan to a grouping by sample/file.

[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
extracts the intensity values from each scan of all files within an
`XCMSnExp` object. These values are extracted from the original data
files and eventual processing steps are applied *on the fly*. Using the
`bySample` parameter it is possible to switch from the default grouping
of intensity values by spectrum/scan to a grouping by sample/file.

[`spectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
extracts the `Spectrum` objects containing all data from `object`. The
values are extracted from the original data files and eventual
processing steps are applied *on the fly*. By setting `bySample = TRUE`,
the spectra are returned grouped by sample/file. If the `XCMSnExp`
object contains adjusted retention times, these are returned by default
in the `Spectrum` objects (can be overwritten by setting
`adjusted = FALSE`).

`processHistory()`: returns a `list` of
[`ProcessHistory()`](https://sneumann.github.io/xcms/reference/ProcessHistory-class.md)
objects (or objects inheriting from this base class) representing the
individual processing steps that have been performed, eventually along
with their settings (`Param` parameter class). Optional arguments
`fileIndex`, `type` and `msLevel` allow to restrict to process steps of
a certain type or performed on a certain file or MS level.

`dropChromPeaks()`: drops any identified chromatographic peaks and
returns the object without that information. Note that for `XCMSnExp`
objects the method drops by default also results from a correspondence
(peak grouping) analysis. Adjusted retention times are removed if the
alignment has been performed *after* peak detection. This can be
overruled with `keepAdjustedRtime = TRUE`.

`dropFeatureDefinitions()`: drops the results from a correspondence
(peak grouping) analysis, i.e. the definition of the mz-rt features and
returns the object without that information. Note that for `XCMSnExp`
objects the method will also by default drop retention time adjustment
results, if these were performed after the last peak grouping (i.e.
which base on the results from the peak grouping that are going to be
removed). All related process history steps are removed too as well as
eventually filled in peaks (by
[`fillChromPeaks()`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)).
The parameter `keepAdjustedRtime` can be used to avoid removal of
adjusted retention times.

`dropAdjustedRtime()`: drops any retention time adjustment information
and returns the object without adjusted retention time. For `XCMSnExp`
objects, this also reverts the retention times reported for the
chromatographic peaks in the peak matrix to the original, raw, ones
(after chromatographic peak detection). Note that for `XCMSnExp` objects
the method drops also all peak grouping results if these were performed
*after* the retention time adjustment. All related process history steps
are removed too.

[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
performs chromatographic peak detection on the provided `XCMSnExp`
objects. For more details see the method for `XCMSnExp()`. Note that by
default (with parameter `add = FALSE`) previous peak detection results
are removed. Use `add = TRUE` to perform a second round of peak
detection and add the newly identified peaks to the previous peak
detection results. Correspondence results (features) are always removed
prior to peak detection. Previous alignment (retention time adjustment)
results are kept, i.e. chromatographic peak detection is performed using
adjusted retention times if the data was first aligned using e.g.
obiwarp
([`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)).

`dropFilledChromPeaks()`: drops any filled-in chromatographic peaks
(filled in by the
[`fillChromPeaks()`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)
method) and all related process history steps.

[`spectrapply()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
applies the provided function to each `Spectrum` in the object and
returns its results. If no function is specified the function simply
returns the `list` of `Spectrum` objects.

`XCMSnExp` objects can be combined with the
[`c()`](https://rdrr.io/r/base/c.html) function. This combines
identified chromatographic peaks and the objects' pheno data but
discards alignment results or feature definitions.

[`plot()`](https://rdrr.io/r/graphics/plot.default.html) plots the
spectrum data (see
[`MSnbase::plot()`](https://lgatto.github.io/MSnbase/reference/plot-methods.html)
for `MSnExp` objects in the *MSnbase* package for more details. For
`type = "XIC"`, identified chromatographic peaks will be indicated as
rectangles with border color `peakCol`.

## Usage

``` r
processHistoryTypes()

# S4 method for class 'XCMSnExp'
hasFilledChromPeaks(object)

# S4 method for class 'OnDiskMSnExp'
profMat(
  object,
  method = "bin",
  step = 0.1,
  baselevel = NULL,
  basespace = NULL,
  mzrange. = NULL,
  fileIndex,
  ...
)

# S4 method for class 'XCMSnExp'
hasAdjustedRtime(object)

# S4 method for class 'XCMSnExp'
hasFeatures(object, msLevel = integer())

# S4 method for class 'XCMSnExp'
hasChromPeaks(object, msLevel = integer())

# S4 method for class 'XCMSnExp'
hasFilledChromPeaks(object)

# S4 method for class 'XCMSnExp'
adjustedRtime(object, bySample = FALSE)

# S4 method for class 'XCMSnExp'
adjustedRtime(object) <- value

# S4 method for class 'XCMSnExp'
featureDefinitions(
  object,
  mz = numeric(),
  rt = numeric(),
  ppm = 0,
  type = c("any", "within", "apex_within"),
  msLevel = integer()
)

# S4 method for class 'XCMSnExp'
featureDefinitions(object) <- value

# S4 method for class 'XCMSnExp'
chromPeaks(
  object,
  bySample = FALSE,
  rt = numeric(),
  mz = numeric(),
  ppm = 0,
  msLevel = integer(),
  type = c("any", "within", "apex_within"),
  isFilledColumn = FALSE
)

# S4 method for class 'XCMSnExp'
chromPeaks(object) <- value

# S4 method for class 'XCMSnExp'
rtime(object, bySample = FALSE, adjusted = hasAdjustedRtime(object))

# S4 method for class 'XCMSnExp'
mz(object, bySample = FALSE, BPPARAM = bpparam())

# S4 method for class 'XCMSnExp'
intensity(object, bySample = FALSE, BPPARAM = bpparam())

# S4 method for class 'XCMSnExp'
spectra(
  object,
  bySample = FALSE,
  adjusted = hasAdjustedRtime(object),
  BPPARAM = bpparam()
)

# S4 method for class 'XCMSnExp'
processHistory(object, fileIndex, type, msLevel)

# S4 method for class 'XCMSnExp'
dropChromPeaks(object, keepAdjustedRtime = FALSE)

# S4 method for class 'XCMSnExp'
dropFeatureDefinitions(object, keepAdjustedRtime = FALSE, dropLastN = -1)

# S4 method for class 'XCMSnExp'
dropAdjustedRtime(object)

# S4 method for class 'XCMSnExp'
profMat(
  object,
  method = "bin",
  step = 0.1,
  baselevel = NULL,
  basespace = NULL,
  mzrange. = NULL,
  fileIndex,
  ...
)

# S4 method for class 'XCMSnExp,Param'
findChromPeaks(
  object,
  param,
  BPPARAM = bpparam(),
  return.type = "XCMSnExp",
  msLevel = 1L,
  add = FALSE
)

# S4 method for class 'XCMSnExp'
dropFilledChromPeaks(object)

# S4 method for class 'XCMSnExp'
spectrapply(object, FUN = NULL, BPPARAM = bpparam(), ...)

# S3 method for class 'XCMSnExp'
c(...)

# S4 method for class 'XCMSnExp'
chromPeakData(object, columns = character(), ...)

# S4 method for class 'XCMSnExp'
chromPeakData(object) <- value

# S4 method for class 'XCMSnExp,missing'
plot(x, y, type = c("spectra", "XIC"), peakCol = "#ff000060", ...)
```

## Arguments

- object:

  either a `MsFeatureData` or a `XCMSnExp` object.

- method:

  `character(1)` defining the profile matrix generation method. Allowed
  are `"bin"`, `"binlin"`, `"binlinbase"` and `"intlin"`. See details
  section for more information.

- step:

  `numeric(1)` representing the m/z bin size.

- baselevel:

  `numeric(1)` representing the base value to which empty elements (i.e.
  m/z bins without a measured intensity) should be set. Only considered
  if `method = "binlinbase"`. See `baseValue` parameter of
  [`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  for more details.

- basespace:

  `numeric(1)` representing the m/z length after which the signal will
  drop to the base level. Linear interpolation will be used between
  consecutive data points falling within `2 * basespace` to each other.
  Only considered if `method = "binlinbase"`. If not specified, it
  defaults to `0.075`. Internally this parameter is translated into the
  `distance` parameter of the
  [`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  function by `distance = floor(basespace / step)`. See `distance`
  parameter of
  [`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  for more details.

- mzrange.:

  Optional `numeric(2)` manually specifying the mz value range to be
  used for binnind. If not provided, the whole m/z value range is used.

- fileIndex:

  For `processHistory()`: optional `integer` specifying the index of the
  files/samples for which the
  [`ProcessHistory()`](https://sneumann.github.io/xcms/reference/ProcessHistory-class.md)
  objects should be retrieved.

- ...:

  Additional parameters.

- msLevel:

  `integer` specifying the MS level(s) for which identified
  chromatographic peaks should be returned.

- bySample:

  `logical(1)` specifying whether results should be grouped by sample.

- value:

  For `adjustedRtime<-`: a `list` (length equal to the number of
  samples) with numeric vectors representing the adjusted retention
  times per scan.

      For `featureDefinitions<-`: a `DataFrame` with peak
      grouping information. See return value for the `featureDefinitions`
      method for the expected format.

      For `chromPeaks<-`: a `matrix` with information on
      detected peaks. See return value for the `chromPeaks` method for the
      expected format.

- mz:

  optional `numeric(2)` defining the mz range for which chromatographic
  peaks should be returned.

- rt:

  optional `numeric(2)` defining the retention time range for which
  chromatographic peaks should be returned.

- ppm:

  optional `numeric(1)` specifying the ppm by which the `mz` range
  should be extended. For a value of `ppm = 10`, all peaks within
  `mz[1] - ppm / 1e6` and `mz[2] + ppm / 1e6` are returned.

- type:

  For `processHistory()`: restrict returned
  [`ProcessHistory()`](https://sneumann.github.io/xcms/reference/ProcessHistory-class.md)
  objects to analysis steps of a certain type. Use the
  `processHistoryTypes()` to list all supported values. For
  `chromPeaks()]: `character`specifying which peaks to return if`rt`or`mz`are defined. For`type
  =
  "any"`all chromatographic peaks partially overlapping the range defined by`mz`and/or`rt`are returned,`type
  = "within"`returns only peaks completely within the region and`type =
  "apex_within"\` peaks for which the peak's apex is within the region.

- isFilledColumn:

  `logical(1)` whether a column `"is_filled"` is included in the
  returned
  "matrix"`providing the information if a peak was filled in. Alternatively, this information would be provided by the`chromPeakData()\`
  data frame.

- adjusted:

  `logical(1)` whether adjusted or raw (i.e. the original retention
  times reported in the files) should be returned.

- BPPARAM:

  Parameter class for parallel processing. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html).

- keepAdjustedRtime:

  For `dropFeatureDefinitions,XCMSnExp()`: `logical(1)` defining whether
  eventually present retention time adjustment should not be dropped. By
  default dropping feature definitions drops retention time adjustment
  results too.

- dropLastN:

  For `dropFeatureDefinitions,XCMSnExp()`: `numeric(1)` defining the
  number of peak grouping related process history steps to remove. By
  default `dropLastN = -1`, dropping the chromatographic peaks removes
  all process history steps related to peak grouping. Setting e.g.
  `dropLastN = 1` will only remove the most recent peak grouping related
  process history step.

- param:

  A
  [`CentWaveParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md),
  [`MatchedFilterParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md),
  [`MassifquantParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-massifquant.md),
  [`MSWParam()`](https://sneumann.github.io/xcms/reference/findPeaks-MSW.md)
  or
  [`CentWavePredIsoParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWaveWithPredIsoROIs.md)
  object with the settings for the chromatographic peak detection
  algorithm.

- return.type:

  Character specifying what type of object the method should return. Can
  be either `"XCMSnExp"` (default), `"list"` or `"xcmsSet"`.

- add:

  For
  [`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md):
  if newly identified chromatographic peaks should be added to the peak
  matrix with the already identified chromatographic peaks. By default
  (`add = FALSE`) previous peak detection results will be removed.

- FUN:

  For `spectrapply`: a function that should be applied to each spectrum
  in the object.

- columns:

  For `chromPeakData()`: optional `character` with the names of the
  columns to include in the returned data frame. By default
  (`columns = character()`) all columns are reported.

- x:

  For [`plot()`](https://rdrr.io/r/base/plot.html): `XCMSnExp` object.

- y:

  For [`plot()`](https://rdrr.io/r/base/plot.html): not used.

- peakCol:

  For [`plot()`](https://rdrr.io/r/base/plot.html): the color that
  should be used to indicate identified chromatographic peaks (only in
  combination with `type = "XIC"` and if chromatographic peaks are
  present).

## Value

For
[`profMat()`](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md):
a `list` with a the profile matrix `matrix` (or matrices if `fileIndex`
was not specified or if \`length(fileIndex) \> 1). See
[profile-matrix](https://sneumann.github.io/xcms/reference/profMat-xcmsSet.md)
for general help and information about the profile matrix.

For `adjustedRtime()`: if `bySample = FALSE` a `numeric` vector with the
adjusted retention for each spectrum of all files/samples within the
object. If `bySample = TRUE` a `list` (length equal to the number of
samples) with adjusted retention times grouped by sample. Returns `NULL`
if no adjusted retention times are present.

For `featureDefinitions()`: a `DataFrame` with peak grouping
information, each row corresponding to one mz-rt feature (grouped peaks
within and across samples) and columns `"mzmed"` (median mz value),
`"mzmin"` (minimal mz value), `"mzmax"` (maximum mz value), `"rtmed"`
(median retention time), `"rtmin"` (minimal retention time), `"rtmax"`
(maximal retention time) and `"peakidx"`. Column `"peakidx"` contains a
`list` with indices of chromatographic peaks (rows) in the matrix
returned by the `chromPeaks()` method that belong to that feature group.
The method returns `NULL` if no feature definitions are present.
`featureDefinitions()` supports also parameters `mz`, `rt`, `ppm` and
`type` to return only features within certain ranges (see description of
`chromPeaks()` for details).

For `chromPeaks`: if `bySample = FALSE` a `matrix` (each row being a
chromatographic peak, rownames representing unique IDs of the peaks)
with at least the following columns:

- `"mz"` (intensity-weighted mean of mz values of the peak across
  scans/retention times),

- `"mzmin"` (minimal mz value),

- `"mzmax"` (maximal mz value),

- `"rt"` (retention time of the peak apex),

- `"rtmin"` (minimal retention time),

- `"rtmax"` (maximal retention time),

- `"into"` (integrated, original, intensity of the peak),

- `"maxo"` (maximum intentity of the peak),

- `"sample"` (sample index in which the peak was identified) and

Depending on the employed peak detection algorithm and the
`verboseColumns` parameter of it, additional columns might be returned.
If parameter `isFilledColumn` was set to `TRUE` a column named
`"is_filled"` is also returned. For `bySample = TRUE` the
chromatographic peaks are returned as a `list` of matrices, each
containing the chromatographic peaks of a specific sample. For samples
in which no peaks were detected a matrix with 0 rows is returned.

For [`rtime()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
if `bySample = FALSE` a numeric vector with the retention times of each
scan, if `bySample = TRUE` a `list` of numeric vectors with the
retention times per sample.

For [`mz()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html): if
`bySample = FALSE` a `list` with the mz values (numeric vectors) of each
scan. If `bySample = TRUE` a `list` with the mz values per sample.

For
[`intensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
if `bySample = FALSE` a `list` with the intensity values (numeric
vectors) of each scan. If `bySample = TRUE` a `list` with the intensity
values per sample.

For
[`spectra()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
if `bySample = FALSE` a `list` with `Spectrum` objects. If
`bySample = TRUE` the result is grouped by sample, i.e. as a `list` of
lists, each element in the *outer* `list` being the `list` of spectra of
the specific file.

For `processHistory()`: a `list` of
[`ProcessHistory()`](https://sneumann.github.io/xcms/reference/ProcessHistory-class.md)
objects providing the details of the individual data processing steps
that have been performed.

## Slots

- `.processHistory`:

  `list` with `XProcessHistory` objects tracking all individual analysis
  steps that have been performed.

- `msFeatureData`:

  `MsFeatureData` class extending `environment` and containing the
  results from a chromatographic peak detection (element
  `"chromPeaks"`), peak grouping (element `"featureDefinitions"`) and
  retention time correction (element `"adjustedRtime"`) steps. This
  object should not be manipulated directly.

## Note

The `"chromPeaks"` element in the `msFeatureData` slot is equivalent to
the `@peaks` slot of the `xcmsSet` object, the `"featureDefinitions"`
contains information from the `@groups` and `@groupidx` slots from an
`xcmsSet` object.

## Chromatographic peak data

Chromatographic peak data is added to an `XCMSnExp` object by the
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
function. Functions to access chromatographic peak data are:

- `hasChromPeaks()` whether chromatographic peak data is available, see
  below for help of the function.

- `chromPeaks()` access chromatographic peaks (see below for help).

- `dropChromPeaks()` remove chromatographic peaks (see below for help).

- `dropFilledChromPeaks()` remove filled-in peaks (see below for help).

- `[fillChromPeaks()]` fill-in missing peaks (see respective help page).

- `[plotChromPeaks()]` plot identified peaks for a file (see respective
  help page).

- `[plotChromPeakImage()]` plot distribution of peaks along the
  retention time axis (see respective help page).

## Adjusted retention times

Adjusted retention times are stored in an `XCMSnExp` object besides the
original, raw, retention times, allowing to switch between raw and
adjusted times. It is also possible to replace the raw retention times
with the adjusted ones with the
[`applyAdjustedRtime()`](https://sneumann.github.io/xcms/reference/applyAdjustedRtime.md)
function. The adjusted retention times are added to an `XCMSnExp` by the
[`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
function. All functions related to the access of adjusted retention
times are:

- `hasAdjustedRtime()` whether adjusted retention times are available
  (see below for help).

- `dropAdjustedRtime()` remove adjusted retention times (see below for
  help).

- [`applyAdjustedRtime()`](https://sneumann.github.io/xcms/reference/applyAdjustedRtime.md)
  replace the raw retention times with the adjusted ones (see respective
  help page).

- [`plotAdjustedRtime()`](https://sneumann.github.io/xcms/reference/plotAdjustedRtime.md)
  plot differences between adjusted and raw retention times (see
  respective help page).

## Correspondence results, features

The correspondence analysis
[`groupChromPeaks()`](https://sneumann.github.io/xcms/reference/groupChromPeaks.md)
adds the definition of LC-MS features to an `XCMSnExp` object. All
functions related to these are listed below:

- `hasFeatures()` whether correspondence results are available (see
  below for help).

- `featureDefinitions()` access the definitions of the features (see
  below for help).

- `dropFeatureDefinitions()` remove correspondence results (see below
  for help).

- [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md)
  access values for features (see respective help page).

- [`featureSummary()`](https://sneumann.github.io/xcms/reference/featureSummary.md)
  perform a simple summary of the defined features (see respective help
  page).

- [`overlappingFeatures()`](https://sneumann.github.io/xcms/reference/overlappingFeatures.md)
  identify features that are overlapping or close in the m/z - rt space
  (see respective help page).

- [`quantify()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  extract feature intensities and put them, along with feature
  definitions and phenodata information, into a
  [`SummarizedExperiment::SummarizedExperiment()`](https://rdrr.io/pkg/SummarizedExperiment/man/SummarizedExperiment-class.html).
  See help page for details.

## See also

[MSnbase::OnDiskMSnExp](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html),
and
[MSnbase::pSet](https://lgatto.github.io/MSnbase/reference/pSet-class.html)
for a complete list of inherited methods.

    [findChromPeaks()] for available peak detection methods
    returning a `XCMSnExp` object as a result.

    [groupChromPeaks()] for available peak grouping
    methods and `featureDefinitions` for the method to extract
    the feature definitions representing the peak grouping results.
    [adjustRtime()] for retention time adjustment methods.

    [chromatogram()] to extract chromatographic MS data.

    [featureChromatograms()] to extract chromatograms for each
    feature.

    [chromPeakSpectra()] to extract MS1 or MS2 spectra for each
    chromatographic peak.

    [featureSpectra()] to extract MS1 or MS2 spectra for features.

[`fillChromPeaks()`](https://sneumann.github.io/xcms/reference/fillChromPeaks.md)
for the method to fill-in eventually missing chromatographic peaks for a
feature in some samples.

## Author

Johannes Rainer

## Examples

``` r

## Load a test data set with detected peaks
library(MSnbase)
#> Loading required package: BiocGenerics
#> Loading required package: generics
#> 
#> Attaching package: ‘generics’
#> The following objects are masked from ‘package:base’:
#> 
#>     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
#>     setequal, union
#> 
#> Attaching package: ‘BiocGenerics’
#> The following objects are masked from ‘package:stats’:
#> 
#>     IQR, mad, sd, var, xtabs
#> The following objects are masked from ‘package:base’:
#> 
#>     Filter, Find, Map, Position, Reduce, anyDuplicated, aperm, append,
#>     as.data.frame, basename, cbind, colnames, dirname, do.call,
#>     duplicated, eval, evalq, get, grep, grepl, is.unsorted, lapply,
#>     mapply, match, mget, order, paste, pmax, pmax.int, pmin, pmin.int,
#>     rank, rbind, rownames, sapply, saveRDS, table, tapply, unique,
#>     unsplit, which.max, which.min
#> Loading required package: Biobase
#> Welcome to Bioconductor
#> 
#>     Vignettes contain introductory material; view with
#>     'browseVignettes()'. To cite Bioconductor, see
#>     'citation("Biobase")', and for packages 'citation("pkgname")'.
#> Loading required package: mzR
#> Loading required package: Rcpp
#> Loading required package: S4Vectors
#> Loading required package: stats4
#> 
#> Attaching package: ‘S4Vectors’
#> The following object is masked from ‘package:utils’:
#> 
#>     findMatches
#> The following objects are masked from ‘package:base’:
#> 
#>     I, expand.grid, unname
#> 
#> This is MSnbase version 2.37.0 
#>   Visit https://lgatto.github.io/MSnbase/ to get started.
#>  Consider switching to the 'R for Mass Spectrometry'
#>  packages - see https://RforMassSpectrometry.org for details.
#> 
#> Attaching package: ‘MSnbase’
#> The following object is masked from ‘package:base’:
#> 
#>     trimws
data(faahko_sub)
## Update the path to the files for the local system
dirname(faahko_sub) <- system.file("cdf/KO", package = "faahKO")

## Disable parallel processing for this example
register(SerialParam())

## The results from the peak detection are now stored in the XCMSnExp
## object
faahko_sub
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 1.22 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 3834 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#>  MSnbase version: 2.33.3 
#> - - - Meta data  - - -
#> phenoData
#>   rowNames: ko15.CDF ko16.CDF ko18.CDF
#>   varLabels: sampleNames
#>   varMetadata: labelDescription
#> Loaded from:
#>   [1] ko15.CDF...  [3] ko18.CDF
#>   Use 'fileNames(.)' to see all files.
#> protocolData: none
#> featureData
#>   featureNames: F1.S0001 F1.S0002 ... F3.S1278 (3834 total)
#>   fvarLabels: fileIdx spIdx ... spectrum (35 total)
#>   fvarMetadata: labelDescription
#> experimentData: use 'experimentData(object)'
#> - - - xcms preprocessing - - -
#> Chromatographic peak detection:
#>  method: centWave 
#>  248 peaks identified in 3 samples.
#>  On average 82.7 chromatographic peaks per sample.

## The detected peaks can be accessed with the chromPeaks method.
head(chromPeaks(faahko_sub))
#>          mz mzmin mzmax       rt    rtmin    rtmax       into       intb   maxo
#> CP001 453.2 453.2 453.2 2506.073 2501.378 2527.982  1007409.0  1007380.8  38152
#> CP002 302.0 302.0 302.0 2617.185 2595.275 2640.659   687146.6   671297.8  30552
#> CP003 344.0 344.0 344.0 2679.783 2646.919 2709.517  5210015.9  5135916.9 152320
#> CP004 430.1 430.1 430.1 2681.348 2639.094 2712.647  2395840.3  2299899.6  65752
#> CP005 366.0 366.0 366.0 2679.783 2642.224 2718.907  3365174.0  3279468.3  79928
#> CP006 343.0 343.0 343.0 2678.218 2637.529 2712.647 24147443.2 23703761.7 672064
#>          sn sample
#> CP001 38151      1
#> CP002    46      1
#> CP003    68      1
#> CP004    42      1
#> CP005    49      1
#> CP006    87      1

## The settings of the chromatographic peak detection can be accessed with
## the processHistory method
processHistory(faahko_sub)
#> [[1]]
#> Object of class "XProcessHistory"
#>  type: Peak detection 
#>  date: Wed Mar 12 08:36:21 2025 
#>  info:  
#>  fileIndex: 1,2,3 
#>  Parameter class: CentWaveParam 
#>  MS level(s) 1 
#> 

## Also the parameter class for the peak detection can be accessed
processParam(processHistory(faahko_sub)[[1]])
#> Object of class:  CentWaveParam 
#>  Parameters:
#>  - ppm: [1] 25
#>  - peakwidth: [1] 20 50
#>  - snthresh: [1] 40
#>  - prefilter: [1]     3 10000
#>  - mzCenterFun: [1] "wMean"
#>  - integrate: [1] 1
#>  - mzdiff: [1] -0.001
#>  - fitgauss: [1] FALSE
#>  - noise: [1] 10000
#>  - verboseColumns: [1] FALSE
#>  - roiList: list()
#>  - firstBaselineCheck: [1] TRUE
#>  - roiScales: numeric(0)
#>  - extendLengthMSW: [1] FALSE
#>  - verboseBetaColumns: [1] FALSE

## The XCMSnExp inherits all methods from the pSet and OnDiskMSnExp classes
## defined in Bioconductor's MSnbase package. To access the (raw) retention
## time for each spectrum we can use the rtime method. Setting bySample = TRUE
## would cause the retention times to be grouped by sample
head(rtime(faahko_sub))
#> F1.S0001 F1.S0002 F1.S0003 F1.S0004 F1.S0005 F1.S0006 
#> 2501.378 2502.943 2504.508 2506.073 2507.638 2509.203 

## Similarly it is possible to extract the mz values or the intensity values
## using the mz and intensity method, respectively, also with the option to
## return the results grouped by sample instead of the default, which is
## grouped by spectrum. Finally, to extract all of the data we can use the
## spectra method which returns Spectrum objects containing all raw data.
## Note that all these methods read the information from the original input
## files and subsequently apply eventual data processing steps to them.
mzs <- mz(faahko_sub, bySample = TRUE)
length(mzs)
#> [1] 3
lengths(mzs)
#>      1      2      3 
#> 445769 440794 442516 

## The full data could also be read using the spectra data, which returns
## a list of Spectrum object containing the mz, intensity and rt values.
## spctr <- spectra(faahko_sub)
## To get all spectra of the first file we can split them by file
## head(split(spctr, fromFile(faahko_sub))[[1]])

############
## Filtering
##
## XCMSnExp objects can be filtered by file, retention time, mz values or
## MS level. For some of these filter preprocessing results (mostly
## retention time correction and peak grouping results) will be dropped.
## Below we filter the XCMSnExp object by file to extract the results for
## only the second file.
xod_2 <- filterFile(faahko_sub, file = 2)
xod_2
#> MSn experiment data ("XCMSnExp")
#> Object size in memory: 0.48 Mb
#> - - - Spectra data - - -
#>  MS level(s): 1 
#>  Number of spectra: 1278 
#>  MSn retention times: 41:41 - 74:60 minutes
#> - - - Processing information - - -
#> Data loaded [Wed Mar 12 08:36:21 2025] 
#> Filter: select file(s) 2. [Wed Apr 22 11:02:32 2026] 
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
#> Chromatographic peak detection:
#>  method: centWave 
#>  100 peaks identified in 1 samples.
#>  On average 100 chromatographic peaks per sample.

## Now the objects contains only the idenfified peaks for the second file
head(chromPeaks(xod_2))
#>          mz mzmin mzmax       rt    rtmin    rtmax       into        intb
#> CP088 384.1 384.1 384.1 2603.100 2592.145 2614.054   274906.3   274885.96
#> CP089 533.3 533.3 533.3 2686.042 2681.347 2689.172    70456.3    70450.04
#> CP090 344.0 344.0 344.0 2686.042 2648.483 2714.211  5700652.0  5574744.87
#> CP091 366.0 366.0 366.0 2686.042 2639.094 2725.166  3448376.1  3338483.07
#> CP092 343.0 343.0 343.0 2686.042 2645.353 2717.341 26229546.1 25615117.17
#> CP093 365.0 365.0 365.0 2686.042 2639.094 2726.731 15565868.0 14884380.84
#>         maxo    sn sample
#> CP088  15187 15186      1
#> CP089  11819 11818      1
#> CP090 151744    87      1
#> CP091  77640    42      1
#> CP092 692352    68      1
#> CP093 345536    48      1

##########
## Coercing to an xcmsSet object
##
## We can also coerce the XCMSnExp object into an xcmsSet object:
xs <- as(faahko_sub, "xcmsSet")
#> Note: you might want to set/adjust the 'sampclass' of the returned xcmSet object before proceeding with the analysis.
head(peaks(xs))
#>          mz mzmin mzmax       rt    rtmin    rtmax       into       intb   maxo
#> CP001 453.2 453.2 453.2 2506.073 2501.378 2527.982  1007409.0  1007380.8  38152
#> CP002 302.0 302.0 302.0 2617.185 2595.275 2640.659   687146.6   671297.8  30552
#> CP003 344.0 344.0 344.0 2679.783 2646.919 2709.517  5210015.9  5135916.9 152320
#> CP004 430.1 430.1 430.1 2681.348 2639.094 2712.647  2395840.3  2299899.6  65752
#> CP005 366.0 366.0 366.0 2679.783 2642.224 2718.907  3365174.0  3279468.3  79928
#> CP006 343.0 343.0 343.0 2678.218 2637.529 2712.647 24147443.2 23703761.7 672064
#>          sn sample
#> CP001 38151      1
#> CP002    46      1
#> CP003    68      1
#> CP004    42      1
#> CP005    49      1
#> CP006    87      1
```
