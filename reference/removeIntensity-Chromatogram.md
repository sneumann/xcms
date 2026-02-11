# Remove intensities from chromatographic data

`removeIntensities` allows to remove intensities from chromatographic
data matching certain conditions (depending on parameter `which`). The
intensities are actually not *removed* but replaced with `NA_real_`. To
actually **remove** the intensities (and the associated retention times)
use
[`MSnbase::clean()`](https://lgatto.github.io/MSnbase/reference/clean-methods.html)
afterwards.

Parameter `which` allows to specify which intensities should be replaced
by `NA_real_`. By default (`which = "below_threshod"` intensities below
`threshold` are removed. If `x` is a `XChromatogram` or `XChromatograms`
object (and hence provides also chromatographic peak definitions within
the object) `which = "outside_chromPeak"` can be selected which removes
any intensity which is outside the boundaries of identified
chromatographic peak(s) in the chromatographic data.

Note that
[`filterIntensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
might be a better approach to subset/filter chromatographic data.

## Usage

``` r
# S4 method for class 'Chromatogram'
removeIntensity(object, which = "below_threshold", threshold = 0)

# S4 method for class 'MChromatograms'
removeIntensity(object, which = "below_threshold", threshold = 0)

# S4 method for class 'XChromatogram'
removeIntensity(
  object,
  which = c("below_threshold", "outside_chromPeak"),
  threshold = 0
)
```

## Arguments

- object:

  an object representing chromatographic data. Can be a
  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html),
  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html),
  [`XChromatogram()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  or
  [`XChromatograms()`](https://sneumann.github.io/xcms/reference/XChromatogram.md)
  object.

- which:

  `character(1)` defining the condition to remove intensities. See
  description for details and options.

- threshold:

  `numeric(1)` defining the threshold below which intensities are
  removed (if `which = "below_threshold"`).

## Value

the input object with matching intensities being replaced by `NA`.

## Author

Johannes Rainer

## Examples

``` r

library(MSnbase)
chr <- Chromatogram(rtime = 1:10 + rnorm(n = 10, sd = 0.3),
    intensity = c(5, 29, 50, NA, 100, 12, 3, 4, 1, 3))

## Remove all intensities below 20
res <- removeIntensity(chr, threshold = 20)
intensity(res)
#>  [1]  NA  29  50  NA 100  NA  NA  NA  NA  NA
```
