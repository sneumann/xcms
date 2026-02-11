# Class xcmsSet, a class for preprocessing peak data

This class transforms a set of peaks from multiple LC/MS or GC/MS
samples into a matrix of preprocessed data. It groups the peaks and does
nonlinear retention time correction without internal standards. It fills
in missing peak values from raw data. Lastly, it generates extracted ion
chromatograms for ions of interest.

## Objects from the Class

Objects can be created with the
[`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet.md)
constructor which gathers peaks from a set NetCDF files. Objects can
also be created by calls of the form `new("xcmsSet", ...)`.

## Slots

- peaks:

  `matrix` containing peak data.

- filled:

  A vector with peak indices of peaks which have been added by a
  [`fillPeaks`](https://sneumann.github.io/xcms/reference/fillPeaks-methods.md)
  method.

- groups:

  Matrix containing statistics about peak groups.

- groupidx:

  List containing indices of peaks in each group.

- phenoData:

  A `data.frame` containing the experimental design factors.

- rt:

  `list` containing two lists, `raw` and `corrected`, each containing
  retention times for every scan of every sample.

- filepaths:

  Character vector with absolute path name of each NetCDF file.

- profinfo:

  `list` containing the values `method` - profile generation method, and
  `step` - profile m/z step size and eventual additional parameters to
  the profile function.

- dataCorrection:

  `logical` vector filled if the waters Lock mass correction parameter
  is used.

- polarity:

  A string ("positive" or "negative" or NULL) describing whether only
  positive or negative scans have been used reading the raw data.

- progressInfo:

  Progress informations for some xcms functions (for GUI).

- progressCallback:

  Function to be called, when progressInfo changes (for GUI).

- mslevel:

  Numeric representing the MS level on which the peak picking was
  performed (by default on MS1). This slot should be accessed through
  its getter method `mslevel`.

- scanrange:

  Numeric of length 2 specifying the scan range (or `NULL` for the full
  range). This slot should be accessed through its getter method
  `scanrange`. The scan range provided in this slot represents the scans
  to which the whole raw data is subsetted.

- .processHistory:

  Internal slot to be used to keep track of performed processing steps.
  This slot should not be directly accessed by the user.

## Methods

- [c](https://sneumann.github.io/xcms/reference/c.xcmsSet.md):

  `signature("xcmsSet")`: combine objects together

- filepaths\<-:

  `signature(object = "xcmsSet")`: set `filepaths` slot

- filepaths:

  `signature(object = "xcmsSet")`: get `filepaths` slot

- [diffreport](https://sneumann.github.io/xcms/reference/diffreport-methods.md):

  `signature(object = "xcmsSet")`: create report of differentially
  regulated ions including EICs

- [fillPeaks](https://sneumann.github.io/xcms/reference/fillPeaks-methods.md):

  `signature(object = "xcmsSet")`: fill in peak data for groups with
  missing peaks

- [getEIC](https://sneumann.github.io/xcms/reference/getEIC-methods.md):

  `signature(object = "xcmsSet")`: get list of EICs for each sample in
  the set

- [getXcmsRaw](https://sneumann.github.io/xcms/reference/getXcmsRaw-methods.md):

  `signature(object = "xcmsSet", sampleidx = 1, profmethod = profMethod(object), profstep = profStep(object), profparam=profinfo(object), mslevel = NULL, scanrange = NULL, rt=c("corrected", "raw"), BPPARAM = bpparam())`:
  read the raw data for one or more files in the `xcmsSet` and return
  it. The default parameters will apply all settings used in the
  original
  [`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet.md) call
  to generate the `xcmsSet` object to be applied also to the raw data.
  Parameter `sampleidx` allows to specify which raw file(s) should be
  loaded. Argument `BPPARAM` allows to setup parallel processing.

- groupidx\<-:

  `signature(object = "xcmsSet")`: set `groupidx` slot

- groupidx:

  `signature(object = "xcmsSet")`: get `groupidx` slot

- [groupnames](https://sneumann.github.io/xcms/reference/groupnames-methods.md):

  `signature(object = "xcmsSet")`: get textual names for peak groups

- groups\<-:

  `signature(object = "xcmsSet")`: set `groups` slot

- groups:

  `signature(object = "xcmsSet")`: get `groups` slot

- [groupval](https://sneumann.github.io/xcms/reference/groupval-methods.md):

  `signature(object = "xcmsSet")`: get matrix of values from peak data
  with a row for each peak group

- [group](https://sneumann.github.io/xcms/reference/group-methods.md):

  `signature(object = "xcmsSet")`: find groups of peaks across samples
  that share similar m/z and retention times

- mslevel:

  Getter method for the `mslevel` slot.

- peaks\<-:

  `signature(object = "xcmsSet")`: set `peaks` slot

- peaks:

  `signature(object = "xcmsSet")`: get `peaks` slot

- [plotrt](https://sneumann.github.io/xcms/reference/plotrt-methods.md):

  `signature(object = "xcmsSet")`: plot retention time deviation
  profiles

- profinfo\<-:

  `signature(object = "xcmsSet")`: set `profinfo` slot

- profinfo:

  `signature(object = "xcmsSet")`: get `profinfo` slot

- profMethod:

  `signature(object = "xcmsSet")`: extract the method used to generate
  the profile matrix.

- profStep:

  `signature(object = "xcmsSet")`: extract the profile step used for the
  generation of the profile matrix.

- [retcor](https://sneumann.github.io/xcms/reference/retcor-methods.md):

  `signature(object = "xcmsSet")`: use initial grouping of peaks to do
  nonlinear loess retention time correction

- sampclass\<-:

  `signature(object = "xcmsSet")`: Replaces the column “class” in the
  `phenoData` slot. See details for more information.

- sampclass:

  `signature(object = "xcmsSet")`: Returns the content of the column
  “class” from the `phenoData` slot or, if not present, the interaction
  of the experimental design factors (i.e. of the `phenoData`
  `data.frame`). See details for more information.

- phenoData\<-:

  `signature(object = "xcmsSet")`: set the `phenoData` slot

- phenoData:

  `signature(object = "xcmsSet")`: get the `phenoData` slot

- progressCallback\<-:

  `signature(object = "xcmsSet")`: set the `progressCallback` slot

- progressCallback:

  `signature(object = "xcmsSet")`: get the `progressCallback` slot

- scanrange:

  Getter method for the `scanrange` slot. See scanrange slot description
  above for more details.

- sampnames\<-:

  `signature(object = "xcmsSet")`: set rownames in the `phenoData` slot

- [sampnames](https://sneumann.github.io/xcms/reference/sampnames-methods.md):

  `signature(object = "xcmsSet")`: get rownames in the `phenoData` slot

- [split](https://sneumann.github.io/xcms/reference/split.xcmsSet.md):

  `signature("xcmsSet")`: divide the xcmsSet into a list of xcmsSet
  objects depending on the provided factor. Note that only peak data
  will be preserved, i.e. eventual peak grouping information will be
  lost.

- `object$name`, `object$name<-value`:

  Access and set `name` column in `phenoData`

- `object[, i]`:

  Conducts subsetting of a `xcmsSet` instance. Only subsetting on
  columns, i.e. samples, is supported. Subsetting is performed on all
  slots, also on `groups` and `groupidx`. Parameter `i` can be an
  integer vector, a logical vector or a character vector of sample names
  (matching `sampnames`).

## Details

The `phenoData` slot (and `phenoData` parameter in the
[`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet.md)
function) is intended to contain a `data.frame` describing all
experimental factors, i.e. the samples along with their properties. If
this `data.frame` contains a column named “class”, this will be returned
by the `sampclass` method and will thus be used by all methods to
determine the sample grouping/class assignment (e.g. to define the
colors in various plots or for the
[`group`](https://sneumann.github.io/xcms/reference/group-methods.md)
method).

The `sampclass<-` method adds or replaces the “class” column in the
`phenoData` slot. If a `data.frame` is submitted to this method, the
interaction of its columns will be stored into the “class” column.

Also, similar to other classes in Bioconductor, the `$` method can be
used to directly access all columns in the `phenoData` slot (e.g. use
`xset$name` on a `xcmsSet` object called “xset” to extract the values
from a column named “name” in the `phenoData` slot).

## Author

Colin A. Smith, <csmith@scripps.edu>, Johannes Rainer
<johannes.rainer@eurac.edu>

## See also

[`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet.md)
