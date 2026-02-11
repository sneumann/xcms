# matchedFilter-based peak detection in purely chromatographic data

`findChromPeaks` on a
[`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
or
[`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
object with a
[MatchedFilterParam](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
parameter object performs matchedFilter-based peak detection on purely
chromatographic data. See
[matchedFilter](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
for details on the method and
[MatchedFilterParam](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
for details on the parameter class. Note that not all settings from the
`MatchedFilterParam` will be used. See
[`peaksWithMatchedFilter()`](https://sneumann.github.io/xcms/reference/peaksWithMatchedFilter.md)
for the arguments used for peak detection on purely chromatographic
data.

## Usage

``` r
# S4 method for class 'Chromatogram,MatchedFilterParam'
findChromPeaks(object, param, ...)
```

## Arguments

- object:

  a
  [`MSnbase::Chromatogram()`](https://lgatto.github.io/MSnbase/reference/Chromatogram-class.html)
  or
  [`MSnbase::MChromatograms()`](https://lgatto.github.io/MSnbase/reference/MChromatograms-class.html)
  object.

- param:

  a
  [MatchedFilterParam](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
  object specifying the settings for the peak detection. See
  [`peaksWithMatchedFilter()`](https://sneumann.github.io/xcms/reference/peaksWithMatchedFilter.md)
  for the description of arguments used for peak detection.

- ...:

  currently ignored.

## Value

If called on a `Chromatogram` object, the method returns a `matrix` with
the identified peaks. Columns `"mz"`, `"mzmin"` and `"mzmax"` in the
[`chromPeaks()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
peak matrix provide the mean m/z and the maximum and minimum m/z value
of the `Chromatogram` object. See
[`peaksWithMatchedFilter()`](https://sneumann.github.io/xcms/reference/peaksWithMatchedFilter.md)
for details on the remaining columns.

## See also

[`peaksWithMatchedFilter()`](https://sneumann.github.io/xcms/reference/peaksWithMatchedFilter.md)
for the downstream function and
[matchedFilter](https://sneumann.github.io/xcms/reference/findChromPeaks-matchedFilter.md)
for details on the method.

## Author

Johannes Rainer

## Examples

``` r

## Loading a test data set with identified chromatographic peaks
faahko_sub <- loadXcmsData("faahko_sub2")
faahko_sub <- filterRt(faahko_sub, c(2500, 3700))
#> Filter spectra

##
od <- as(filterFile(faahko_sub, 1L), "MsExperiment")

## Extract chromatographic data for a small m/z range
chr <- chromatogram(od, mz = c(272.1, 272.3))[1, 1]

## Identify peaks with default settings
xchr <- findChromPeaks(chr, MatchedFilterParam())

## Plot the identified peaks
plot(xchr)
```
