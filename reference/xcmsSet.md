# Constructor for xcmsSet objects which finds peaks in NetCDF/mzXML files

This function handles the construction of xcmsSet objects. It finds
peaks in batch mode and pre-sorts files from subdirectories into
different classes suitable for grouping.

## Usage

``` r
xcmsSet(files = NULL, snames = NULL, sclass = NULL, phenoData = NULL,
        profmethod = "bin", profparam = list(),
        polarity = NULL, lockMassFreq=FALSE,
  mslevel=NULL, nSlaves=0, progressCallback=NULL,
        scanrange = NULL, BPPARAM = bpparam(),
        stopOnError = TRUE, ...)
```

## Arguments

- files:

  path names of the NetCDF/mzXML files to read

- snames:

  sample names. By default the file name without extension is used.

- sclass:

  sample classes.

- phenoData:

  `data.frame` or `AnnotatedDataFrame` defining the sample names and
  classes and other sample related properties. If not provided, the
  argument `sclass` or the subdirectories in which the samples are
  stored will be used to specify sample grouping.

- profmethod:

  Method to use for profile generation. Supported values are `"bin"`,
  `"binlin"`, `"binlinbase"` and `"intlin"` (for methods
  [`profBin`](https://sneumann.github.io/xcms/reference/profGenerate.md),
  [`profBinLin`](https://sneumann.github.io/xcms/reference/profGenerate.md),
  [`profBinLinBase`](https://sneumann.github.io/xcms/reference/profGenerate.md)
  and
  [`profIntLin`](https://sneumann.github.io/xcms/reference/profGenerate.md),
  respectively). See help on
  [`profBin`](https://sneumann.github.io/xcms/reference/profGenerate.md)
  for a complete list of available methods and their supported
  parameters.

- profparam:

  parameters to use for profile generation.

- polarity:

  filter raw data for positive/negative scans

- lockMassFreq:

  Performs correction for Waters LockMass function

- mslevel:

  perform peak picking on data of given mslevel

- nSlaves:

  *DEPRECATED*, use `BPPARAM` argument instead.

- progressCallback:

  function to be called, when progressInfo changes (useful for GUIs)

- scanrange:

  scan range to read

- BPPARAM:

  a `BiocParallel` parameter object to control how and if parallel
  processing should be performed. Such objects can be created by the
  `SerialParam`, `MulticoreParam` or `SnowParam` functions.

- stopOnError:

  Logical specifying whether the feature detection call should stop on
  the first encountered error (the default), or whether feature
  detection is performed in all files regardless eventual failures for
  individual files in which case all errors are reported as `warnings`.

- ...:

  further arguments to the `findPeaks` method of the `xcmsRaw` class

## Details

The default values of the `files`, `snames`, `sclass`, and `phenoData`
arguments cause the function to recursively search for readable files.
The filename without extention is used for the sample name. The
subdirectory path is used for the sample class. If the files contain
both positive and negative spectra, the polarity can be selected
explicitly. The default (NULL) is to read all scans.

If `phenoData` is provided, it is stored to the `phenoData` slot of the
returned `xcmsSet` class. If that `data.frame` contains a column named
“class”, its content will be returned by the
[`sampclass`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
method and thus be used for the group/class assignment of the individual
files (e.g. for peak grouping etc.). For more details see the help of
the
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md).

The step size (in m/z) to use for profile generation can be submitted
either using the `profparam` argument (e.g. `profparam=list(step=0.1)`)
or by submitting `step=0.1`. By specifying a value of `0` the profile
matrix generation can be skipped.

The feature/peak detection algorithm can be specified with the `method`
argument which defaults to the `"matchFilter"` method
([`findPeaks.matchedFilter`](https://sneumann.github.io/xcms/reference/findPeaks.matchedFilter-xcmsRaw-method.md)).
Possible values are returned by
`getOption("BioC")$xcms$findPeaks.methods`.

The lock mass correction allows for the lock mass scan to be added back
in with the last working scan. This correction gives better
reproducibility between sample sets.

## Value

A `xcmsSet` object.

## Author

Colin A. Smith, <csmith@scripps.edu>

## Note

The arguments `profmethod` and `profparam` have no influence on the
feature/peak detection. The step size parameter `step` for the profile
generation in the
[`findPeaks.matchedFilter`](https://sneumann.github.io/xcms/reference/findPeaks.matchedFilter-xcmsRaw-method.md)
peak detection algorithm can be passed using the `...`.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md),
[`profStep`](https://sneumann.github.io/xcms/reference/profStep-methods.md),
[`profMethod`](https://sneumann.github.io/xcms/reference/profMethod-methods.md),
[`profBin`](https://sneumann.github.io/xcms/reference/profGenerate.md)
