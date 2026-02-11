# Get and set m/z step for generating profile data

These methods get and set the m/z step for generating profile (matrix)
data from raw mass spectral data. Smaller steps yield more precision at
the cost of greater memory usage.

## Methods

- object = "xcmsRaw":

  `profStep(object)`

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`profMethod`](https://sneumann.github.io/xcms/reference/profMethod-methods.md)

## Examples

``` r
  if (FALSE) { # \dontrun{
    library(faahKO)
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xset <- xcmsRaw(cdffiles[1])

    xset
    plotSurf(xset, mass=c(200,500))

    profStep(xset)<-0.1 ## decrease the bin size to get better resolution
    plotSurf(xset, mass=c(200, 500))
    ##works nicer on high resolution data.
  } # }
```
