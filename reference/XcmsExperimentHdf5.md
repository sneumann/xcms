# xcms result object for very large data sets

The *xcms* result objects
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
and
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
keep all preprocessing results in memory and can thus (depending on the
size of the data set) require a large amount of memory. In contrast, the
`XcmsExperimentHdf5` class, by using an on-disk data storage mechanism,
has a much lower memory footprint allowing also to analyze very large
data sets on regular computer systems such as desktop or laptop
computers. With some exceptions, including additional parameters, the
functionality and usability of this object is identical to the default
`XcmsExperiment` object.

This help page lists functions that have additional or different
parameters or properties than the respective methods for
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
objects. For all other functions not listed here the usability is
identical to those for the
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
object (see the respective help page for information).

## Usage

``` r
toXcmsExperimentHdf5(object, hdf5File = tempfile())

toXcmsExperiment(object, ...)

# S4 method for class 'XcmsExperimentHdf5'
chromPeakData(
  object,
  msLevel = integer(),
  peaks = character(),
  columns = character(),
  return.type = c("DataFrame", "data.frame"),
  bySample = FALSE
)

# S4 method for class 'XcmsExperimentHdf5'
filterChromPeaks(
  object,
  keep = rep(TRUE, nrow(chromPeaks(object))),
  method = "keep",
  ...
)

# S4 method for class 'XcmsExperimentHdf5,PeakGroupsParam'
adjustRtimePeakGroups(object, param = PeakGroupsParam(), msLevel = 1L)

# S4 method for class 'XcmsExperimentHdf5'
filterFeatureDefinitions(object, features = integer())
```

## Arguments

- object:

  `XcmsExperimentHdf5` object.

- hdf5File:

  For `toXcmsExperimentHdf5()`: `character(1)` with the path and name of
  the (not yet existing) file where the preprocessing results should be
  stored to.

- ...:

  additional parameters eventually passed to downstream functions.

- msLevel:

  For
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  and
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md):
  optional `integer` with the MS level(s) from which the data should be
  returned. By default `msLevel = integer()` results from all MS levels
  are returned (if present). For
  [`refineChromPeaks()`](https://sneumann.github.io/xcms/reference/refineChromPeaks.md):
  `integer(1)` with the MS level from which chromatographic peaks should
  be refined.

- peaks:

  For
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md):
  optional `character` with the ID of chromatographic peaks (row name in
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md))
  for which the data should be returned. By default
  (`peaks = character()`) the data for all chromatographic peaks is
  returned.

- columns:

  For
  `chromPeakData()~: optional `character` allowing to define a subset of columns that should be included in the returned data frame. By default (`columns
  = character()\`) the full data is returned.

- return.type:

  For
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md):
  `character(1)` specifying the type of object that should be returned.
  Can be either `return.type = "DataFrame"` (the default) to return a
  `DataFrame`, or `return.type = "data.frame"` to return the results as
  a `data.frame`.

- bySample:

  For
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  and
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md):
  `logical(1)` whether the data should be returned *by sample*, i.e. as
  a `list` of `matrix` or `data.frame` objects, one for each sample.

- keep:

  For
  [`filterChromPeaks()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  defining the chromatographic peaks to keep: either a `logical` with
  the same length than the number of chromatographic peaks, an `integer`
  with the indices or a `character` with the IDs of the chromatographic
  peaks to keep.

- method:

  For
  [`filterChromPeaks()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  `character(1)`; currently only `method = "keep"` is supported.

- param:

  *parameter* object defining and configuring the algorithm to be used.

- features:

  For
  [`filterFeatureDefinitions()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  defining the features to keep: either a `logical` with the same length
  than the number of features, an `integer` with the indices or a
  `character` with the ID of the features to keep.

## Value

See description of the individual methods for information.

## Details

The `XcmsExperimentHdf5` object stores all preprocessing results (except
adjusted retention times, which are stored as an additional spectra
variable in the object's
[`Spectra::Spectra()`](https://rdrr.io/pkg/Spectra/man/Spectra.html)
object), in a file in HDF5 format.

`XcmsExperimentHdf5` uses a different naming scheme for chromatographic
peaks: for efficiency reasons, chromatographic peak data is organized by
sample and MS level. The chrom peak IDs are hence in the format
`"CP"<MS level>"S"<sample id><chrom peak index>` with `<MS level>` being
the MS level in which the chromatographic peaks were detected and
`<sample id>` the ID of the sample (usually related to the index in the
original `MsExperiment` object) and the `<chrom peak index>` the index
of the chromatographic peak in the chrom peak matrix **of that sample**
and MS level.

HDF5 files do not support parallel processing, thus preprocessing
results need to be stored or loaded sequentially.

All functionality for `XcmsExperimentHdf5` objects is optimized to
reduce memory demand at the cost of eventually lower performance.

## Conversion between `XcmsExperiment` and `XcmsExperimentHdf5`

To use the `XcmsExperimentHdf5` class for preprocessing results, the
`hdf5File` parameter of the
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
function needs to be defined, specifying the path and name of the HDF5
file to store the results. In addition it is possible to convert a
`XcmsExperiment` object to a `XcmsExperimentHdf5` object with the
`toXcmsExperimentHdf5()` function. All present preprocessing results
will be stored to the specified HDF5 file. To load all preprocessing
results into memory and hence change from a `XcmsExperimentHdf5` to a
`XcmsExperiment` object, the `toXcmsExperiment()` function can be used.

## Using the HDF5 file-based on-disk data storage

Calling
[`findChromPeaks()`](https://sneumann.github.io/xcms/reference/findChromPeaks.md)
on an `MsExperiment` using the parameter `hdf5File` will return an
instance of the `XcmsExperimentHdf5` class and hence use the on-disk
data storage mode described on this page. The results are stored in the
file specified with parameter `hdf5File`.

## Subset

- `[`: subset the `XcmsExperimentHdf5` object to the specified samples.
  Parameters `keepChromPeaks` (default `TRUE`), `keepAdjustedRtime`
  (default `TRUE`) and `keepFeatures` (default `FALSE`) allow to
  configure whether present chromatographic peaks, alignment or
  correspondence results should be retained. This will only change
  information in the object (i.e., the reference to the respective
  entries in the HDF5 file), but will **not** change the content of the
  HDF5 file. Thus, *reverting* the retention times of detected
  chromatographic peaks is **not** supported and `keepChromPeaks = TRUE`
  with `keepAdjustedRtime = FALSE` will throw an error. Note that with
  `keepChromPeaks = FALSE` also `keepFeatures` is set to `FALSE`.

- [`filterChromPeaks()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  and
  [`filterFeatureDefinitions()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  to filter the chromatographic peak and correspondence results,
  respectively. See documentation below for details. Subset using
  unsorted or duplicated indices is not supported.

## Functionality related to chromatographic peaks

- [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  gains parameter `bySample = FALSE` that, if set to `TRUE` returns a
  `list` of `chromPeaks` matrices, one for each sample. Due to the way
  data is organized in `XcmsExperimentHdf5` objects this is more
  efficient than `bySample = FALSE`. Thus, in cases where chrom peak
  data is subsequently evaluated or processed by sample, it is suggested
  to use `bySample = TRUE`.

- [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  gains a new parameter `peaks = character()` which allows to specify
  from which chromatographic peaks data should be returned. For these
  chromatographic peaks the ID (row name in
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md))
  should be provided with the `peaks` parameter. This can reduce the
  memory requirement for cases in which only data of some selected
  chromatographic peaks needs to be extracted. Also,
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  supports the `bySample` parameter described for
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  above. All other parameters present also for
  [`chromPeakData()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  of `XcmsExperiment` objects, such as `columns` are supported.

- [`filterChromPeaks()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  allows to filter the chromatographic peaks specifying which should be
  retainend using the `keep` parameter. This can be either a `logical`,
  `character` or `integer` vector. Duplicated or unsorted indices are
  **not** supported. Eventually present feature definitions will be
  updated as well. The function returns the object with the filtered
  chromatographic peaks.

## Retention time alignment

- [`adjustRtimePeakGroups()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
  and
  [`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
  with `PeakGroupsParam`: parameter `extraPeaks` of `PeakGroupsParam` is
  **ignored**. Anchor peaks are thus only defined using the
  `minFraction` and the optional `subset` parameter.

## Correspondence analysis results

- [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md):
  similarly to
  [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  for
  [XcmsExperiment](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  objects, this method returns a `data.frame` with the characteristics
  for the defined LC-MS features. The function for `XcmsExperimentHdf5`
  does however **not** return the `"peakidx"` column with the indices of
  the chromatographic peaks per feature. Also, the columns are returned
  in alphabetic order.

- [`featureValues()`](https://sneumann.github.io/xcms/reference/XCMSnExp-peak-grouping-results.md):
  for parameter `value`, the option `value = "index"` (i.e. returning
  the index of the chromatographic peaks within the
  [`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  matrix per feature) is **not** supported.

- [`filterFeatureDefinitions()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md):
  filter the feature definitions keeping only the specified features.
  Parameter `features` can be used to define the features to retain. It
  supports a `logical`, `integer` indices or `character` with the IDs of
  the features (i.e., their row names in
  [`featureDefinitions()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)).
  The function returns the input `XcmsExperimentHdf5` with the filtered
  content.

## Author

Johannes Rainerr, Philippine Louail

## Examples

``` r

## Create a MsExperiment object representing the data from an LC-MS
## experiment.
library(MsExperiment)
library(faahKO)

## Define the raw data files
fls <- c(system.file('cdf/KO/ko15.CDF', package = "faahKO"),
         system.file('cdf/KO/ko16.CDF', package = "faahKO"),
         system.file('cdf/KO/ko18.CDF', package = "faahKO"))

## Define a data frame with the sample characterization
df <- data.frame(mzML_file = basename(fls),
                sample = c("ko15", "ko16", "ko18"))
## Importe the data. This will initialize a `Spectra` object representing
## the raw data and assign these to the individual samples.
mse <- readMsExperiment(spectraFiles = fls, sampleData = df)

## Perform chromatographic peak detection storing the data in an HDF5 file
## Parameter `hdf5File` has to be provided and needs to be the path and
## name of a (not yet existing) file to which results are going to be
## stored. For the example below we use a temporary file.
xmse <- findChromPeaks(mse, param = CentWaveParam(prefilter = c(4, 100000)),
    hdf5File = tempfile())
xmse
#> Object of class XcmsExperimentHdf5 
#>  Spectra: MS1 (3834) 
#>  Experiment data: 3 sample(s)
#>  Sample data links:
#>   - spectra: 3 sample(s) to 3834 element(s).
#>  xcms results:
#>   - chromatographic peaks in MS level(s): 1 
#>  results storage file:
#>    /tmp/RtmpJElXXi/file7a3b78245254

## Extract selected columnds from the chromatographic peak detection
## results
chromPeaks(xmse, columns = c("rt", "mz", "into")) |> head()
#>                   rt  mz     into sample
#> CP1S1000001 2682.913 360  5641322      1
#> CP1S1000002 2679.783 344  5210016      1
#> CP1S1000003 2678.218 343 24147443      1
#> CP1S1000004 2679.783 365 14975761      1
#> CP1S1000005 2659.438 365  3520591      1
#> CP1S1000006 2784.635 280  2537599      1

## Extract the results per sample
res <- chromPeaks(xmse, columns = c("rt", "mz", "into"), bySample = TRUE)

## The chromatographic peaks of the second sample:
res[[2]] |> head()
#>                   rt  mz     into
#> CP1S2000001 2686.042 360 10248211
#> CP1S2000002 2686.042 344  5700652
#> CP1S2000003 2686.042 343 26229546
#> CP1S2000004 2596.840 365  2358688
#> CP1S2000005 2686.042 365 15565868
#> CP1S2000006 2797.154 279 10916521

## Convert the result object to the in-memory representation:
xmse_mem <- toXcmsExperiment(xmse)
xmse_mem
#> Object of class XcmsExperiment 
#>  Spectra: MS1 (3834) 
#>  Experiment data: 3 sample(s)
#>  Sample data links:
#>   - spectra: 3 sample(s) to 3834 element(s).
#>  xcms results:
#>   - chromatographic peaks: 181 in MS level(s): 1 
```
