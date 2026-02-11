# Automatic parameter for Lock mass fixing `AutoLockMass` \~~

`AutoLockMass` - This function decides where the lock mass scans are in
the xcmsRaw object. This is done by using the scan time differences.

## Methods

- object = "xcmsRaw":

  `signature(object = "xcmsRaw")`

## Arguments

- object:

  An
  [`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
  object

## Value

`AutoLockMass` A numeric vector of scan locations corresponding to lock
Mass scans

## Author

Paul Benton, <hpaul.benton08@imperial.ac.uk>

## Examples

``` r
if (FALSE) library(xcms)
    library(faahKO)
    ## These files do not have this problem
    ## to correct for but just for an example
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xr<-xcmsRaw(cdffiles[1])
#> Create profile matrix with method 'bin' and step 1 ... 
#> OK
    xr
#> An "xcmsRaw" object with 1278 mass spectra
#> 
#> Time range: 2501.4-4499.8 seconds (41.7-75 minutes)
#> Mass range: 200-600 m/z
#> Intensity range: 70-1373180 
#> 
#> MSn data on  0  mass(es)
#>  with  0  MSn spectra
#> Profile method: bin 
#> Profile step: 1 m/z (401 grid points from 200 to 600 m/z)
#> 
#> Memory usage: 10.8 MB
    ##Lets assume that the lockmass starts at 1 and is every 100 scans
    lockMass<-xcms:::makeacqNum(xr, freq=100, start=1)
    ## these are equalvent
    lockmass2<-AutoLockMass(xr)
#> Warning: 
#> Lock mass frequency wasn't detected
    all((lockmass == lockmass2) == TRUE)
#> Error: object 'lockmass' not found

    ob<-stitch(xr, lockMass)
#> Warning: no non-missing arguments to min; returning Inf
#> Warning: no non-missing arguments to max; returning -Inf
 # \dontrun{}
```
