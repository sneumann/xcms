# The profile matrix

The *profile* matrix is an n x m matrix, n (rows) representing equally
spaced m/z values (bins) and m (columns) the retention time of the
corresponding scans. Each cell contains the maximum intensity measured
for the specific scan and m/z values falling within the m/z bin.

    The `profMat` method creates a new profile matrix or returns the
    profile matrix within the object's `@env` slot, if available.
    Settings for the profile matrix generation, such as `step` (the bin
    size), `method` or additional settings are extracted from the
    respective slots of the `xcmsRaw` object.
    Alternatively it is possible to specify all of the settings as
    additional parameters.

    For [MsExperiment()] or [XcmsExperiment()] objects, the method returns
    a `list` of profile matrices, one for each sample in `object`. Using
    parameter `fileIndex` it is also possible to create a profile matrix only
    for selected samples (files).

## Usage

``` r
# S4 method for class 'MsExperiment'
profMat(
  object,
  method = "bin",
  step = 0.1,
  baselevel = NULL,
  basespace = NULL,
  mzrange. = NULL,
  fileIndex = seq_along(object),
  chunkSize = 1L,
  msLevel = 1L,
  BPPARAM = bpparam(),
  ...
)

# S4 method for class 'xcmsRaw'
profMat(object, method, step, baselevel, basespace, mzrange.)
```

## Arguments

- object:

  An `xcmsRaw`, `OnDiskMSnExp`, `XCMSnExp`, `MsExperiment` or
  `XcmsExperiment` object.

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

  For `MsExperiment` or `XcmsExperiment`: `integer` defining the idex
  (or indices) of the sample(s) from which the profile matrix should be
  created.

- chunkSize:

  For `MsExperiment` or `XcmsExperiment`: `integer(1)` defining the
  number of files from which data should be loaded and processed in one
  iteration. By default one file at a time is processed `chunkSize = 1L`
  which requires less memory. For parallel processing, the `chunkSize`
  should be `>=` than the number of parallel processes that should be
  used.

- msLevel:

  For `MsExperiment` or `XcmsExperiment`: `integer(1)` defining the MS
  level from which the profile matrix should be generated.

- BPPARAM:

  For `MsExperiment` or `XcmsExperiment`: parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for more details. Defaults to `BPPARAM = bpparam()`.

- ...:

  ignored.

## Value

`profMat` returns the profile matrix (rows representing scans, columns
equally spaced m/z values). For `object` being a `MsExperiment` or
`XcmsExperiment`, the method returns a `list` of profile matrices, one
for each file (sample).

## Details

Profile matrix generation methods:

- `"bin"`: The default profile matrix generation method that does a
  simple binning, i.e. aggregating of intensity values falling within an
  m/z bin.

- `"binlin"`: Binning followed by linear interpolation to impute missing
  values. The value for m/z bins without a measured intensity are
  inferred by a linear interpolation between neighboring bins with a
  measured intensity.

- `"binlinbase"`: Binning followed by a linear interpolation to impute
  values for empty elements (m/z bins) within a user-definable proximity
  to non-empty elements while stetting the element's value to the
  `baselevel` otherwise. See `impute = "linbase"` parameter of
  [`imputeLinInterpol()`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
  for more details.

- `"intlin"`: Set the elements' values to the integral of the linearly
  interpolated data from plus to minus half the step size.

## Author

Johannes Rainer

## Examples

``` r
file <- system.file('cdf/KO/ko15.CDF', package = "faahKO")
## Load the data without generating the profile matrix (profstep = 0)
xraw <- xcmsRaw(file, profstep = 0)
## Extract the profile matrix
profmat <- profMat(xraw, step = 0.3)
#> Create profile matrix with method 'bin' and step 0.3 ... 
#> OK
dim(profmat)
#> [1] 1335 1278
## If not otherwise specified, the settings from the xraw object are used:
profinfo(xraw)
#> $method
#> [1] "bin"
#> 
#> $step
#> [1] 0
#> 
## To extract a profile matrix with linear interpolation use
profmat <- profMat(xraw, step = 0.3, method = "binlin")
#> Create profile matrix with method 'binlin' and step 0.3 ... 
#> OK
## Alternatively, the profMethod of the xraw objects could be changed
profMethod(xraw) <- "binlin"
profmat_2 <- profMat(xraw, step = 0.3)
#> Create profile matrix with method 'binlin' and step 0.3 ... 
#> OK
all.equal(profmat, profmat_2)
#> [1] TRUE
```
