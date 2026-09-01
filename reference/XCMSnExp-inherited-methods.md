# XCMSnExp data manipulation methods inherited from MSnbase

The methods listed on this page are
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
methods inherited from its parent, the
[`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
class from the *MSnbase* package, that alter the raw data or are related
to data subsetting. Thus calling any of these methods causes all *xcms*
pre-processing results to be removed from the
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object to ensure its data integrity.

[`bin()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
allows to *bin* spectra. See
[`MSnbase::bin()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
documentation in the *MSnbase* package for more details and examples.

[`clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.html):
removes unused `0` intensity data points. See
[`MSnbase::clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.html)
documentation in the *MSnbase* package for details and examples.

[`filterAcquisitionNum()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
filters the
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object keeping only spectra with the provided acquisition numbers. See
[`MSnbase::filterAcquisitionNum()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
for details and examples.

The [`normalize()`](https://rdrr.io/pkg/BiocGenerics/man/normalize.html)
method performs basic normalization of spectra intensities. See
[`MSnbase::normalize()`](https://rdrr.io/pkg/BiocGenerics/man/normalize.html)
documentation in the *MSnbase* package for details and examples.

The
[`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html)
method performs peak picking. See documentation for that function in the
*MSnbase* package for details and examples.

The
[`removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.html)
method removes mass peaks (intensities) lower than a threshold. Note
that these peaks refer to *mass* peaks, which are different to the
chromatographic peaks detected and analyzed in a metabolomics
experiment! See
[`MSnbase::removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.html)
documentation for details and examples.

The [`smooth()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
method smooths spectra. See
[`MSnbase::smooth()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
documentation in *MSnbase* for details and examples.

## Usage

``` r
# S4 method for class 'XCMSnExp'
bin(x, binSize = 1L, msLevel.)

# S4 method for class 'XCMSnExp'
clean(object, all = FALSE, verbose = FALSE, msLevel.)

# S4 method for class 'XCMSnExp'
filterAcquisitionNum(object, n, file)

# S4 method for class 'XCMSnExp'
normalize(object, method = c("max", "sum"), ...)

# S4 method for class 'XCMSnExp'
pickPeaks(
  object,
  halfWindowSize = 3L,
  method = c("MAD", "SuperSmoother"),
  SNR = 0L,
  ...
)

# S4 method for class 'XCMSnExp'
removePeaks(object, t = "min", verbose = FALSE, msLevel.)

# S4 method for class 'XCMSnExp'
smooth(
  x,
  method = c("SavitzkyGolay", "MovingAverage"),
  halfWindowSize = 2L,
  verbose = FALSE,
  ...
)
```

## Arguments

- x:

  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  or
  [`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
  object.

- binSize:

  `numeric(1)` defining the size of a bin (in Dalton).

- msLevel.:

  For [`bin()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  [`clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.html),
  [`filterMsLevel()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html),
  [`removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.html):
  `integer(1)` defining the MS level(s) to which operations should be
  applied or to which the object should be subsetted.

- object:

  `XCMSnExp` or `OnDiskMSnExp` object.

- all:

  For
  [`clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.html):
  `logical(1)`, if `TRUE` all zeros are removed.

- verbose:

  `logical(1)` whether progress information should be displayed.

- n:

  For
  [`filterAcquisitionNum()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `integer` defining the acquisition numbers of the spectra to which the
  data set should be sub-setted.

- file:

  For
  [`filterAcquisitionNum()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `integer` defining the file index within the object to subset the
  object by file.

- method:

  For
  [`normalize()`](https://rdrr.io/pkg/BiocGenerics/man/normalize.html):
  `character(1)` specifying the normalization method. See
  [`MSnbase::normalize()`](https://rdrr.io/pkg/BiocGenerics/man/normalize.html)
  in the *MSnbase* package for details. For
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html):
  `character(1)` defining the method. See help for
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html)
  in the *MSnbase* package for options. For
  [`smooth()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `character(1)` defining the method. See
  [`MSnbase::smooth()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  in the *MSnbase* package for options and details.

- ...:

  Optional additional arguments.

- halfWindowSize:

  For
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html)
  and
  [`smooth()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html):
  `integer(1)` defining the window size for the peak picking. See help
  for `pickPeaks` and \[MSnbase::smooth()\` in the *MSnbase* package for
  details and options.

- SNR:

  For
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html):
  `numeric(1)` defining the signal to noise ratio to be considered. See
  the documentation for
  [`pickPeaks()`](https://lgatto.github.io/MSnbase/reference/pickPeaks-method.html)
  in the *MSnbase* package for details.

- t:

  For
  [`removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.html):
  either a `numeric(1)` or `"min"` defining the threshold (method) to be
  used. See
  [`MSnbase::removePeaks()`](https://lgatto.github.io/MSnbase/reference/removePeaks-methods.html)
  for details.

## Value

For all methods: a `XCMSnExp` object.

## See also

[XCMSnExp-filter](https://sneumann.github.io/xcms/reference/XCMSnExp-filter-methods.md)
for methods to filter and subset `XCMSnExp` objects.
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for base class documentation.
[`MSnbase::OnDiskMSnExp()`](https://lgatto.github.io/MSnbase/reference/OnDiskMSnExp-class.html)
for the documentation of the parent class.

## Author

Johannes Rainer
