# Estimate precursor intensity for MS level 2 spectra

[`estimatePrecursorIntensity()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
determines the precursor intensity for a MS 2 spectrum based on the
intensity of the respective signal from the neighboring MS 1 spectra
(i.e. based on the peak with the m/z matching the precursor m/z of the
MS 2 spectrum). Based on parameter `method` either the intensity of the
peak from the previous MS 1 scan is used (`method = "previous"`) or an
interpolation between the intensity from the previous and subsequent MS1
scan is used (`method = "interpolation"`, which considers also the
retention times of the two MS1 scans and the retention time of the MS2
spectrum).

## Usage

``` r
# S4 method for class 'MsExperiment'
estimatePrecursorIntensity(
  object,
  ppm = 10,
  tolerance = 0,
  method = c("previous", "interpolation"),
  BPPARAM = bpparam()
)

# S4 method for class 'OnDiskMSnExp'
estimatePrecursorIntensity(
  object,
  ppm = 10,
  tolerance = 0,
  method = c("previous", "interpolation"),
  BPPARAM = bpparam()
)
```

## Arguments

- object:

  `MsExperiment`, `XcmsExperiment`, `OnDiskMSnExp` or `XCMSnExp` object.

- ppm:

  `numeric(1)` defining the maximal acceptable difference (in ppm) of
  the precursor m/z and the m/z of the corresponding peak in the MS 1
  scan.

- tolerance:

  `numeric(1)` with the maximal allowed difference of m/z values between
  the precursor m/z of a spectrum and the m/z of the respective ion on
  the MS1 scan.

- method:

  `character(1)` defining the method how the precursor intensity should
  be determined (see description above for details). Defaults to
  `method = "previous"`.

- BPPARAM:

  parallel processing setup. See
  [`BiocParallel::bpparam()`](https://rdrr.io/pkg/BiocParallel/man/register.html)
  for details.

## Value

`numeric` with length equal to the number of spectra in `x`. `NA` is
returned for MS 1 spectra or if no matching peak in a MS 1 scan can be
found for an MS 2 spectrum

## Author

Johannes Rainer with feedback and suggestions from Corey Broeckling
