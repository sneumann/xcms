# Copy MSn data in an xcmsRaw to the MS slots

The MS2 and MSn data is stored in separate slots, and can not directly
be used by e.g. findPeaks(). `msn2xcmsRaw()` will copy the MSn spectra
into the "normal" `xcmsRaw` slots.

## Usage

``` r
msn2xcmsRaw(xmsn)
```

## Arguments

- xmsn:

  an object of class `xcmsRaw` that contains spectra read with
  includeMSn=TRUE

## Details

The default gap value is determined from the 90th percentile of the
pair-wise differences between adjacent mass values.

## Value

An xcmsRaw object

## Author

Steffen Neumann <sneumann@ipb-halle.de>

## See also

[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md),

## Examples

``` r
 msnfile <- system.file("microtofq/MSMSpos20_6.mzML", package = "msdata")
 xrmsn <- xcmsRaw(msnfile, includeMSn=TRUE)
#> Create profile matrix with method 'bin' and step 1 ... 
#> OK
 xr <- msn2xcmsRaw(xrmsn)
 p <- findPeaks(xr, method="centWave")
#> Detecting mass traces at 25 ppm ... 
#> OK
#> Detecting chromatographic peaks in 10 regions of interest ...
#>  OK: 10 found.
```
