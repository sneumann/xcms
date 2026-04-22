# Group Peaks via High Resolution Alignment

Runs high resolution alignment on single spectra samples stored in a
given xcmsSet.

## Methods

- object = "xcmsSet":

  ` group(object, method="mzClust", mzppm = 20, mzabs = 0, minsamp = 1, minfrac=0) `

## Arguments

- object:

  a xcmsSet with peaks

- mzppm:

  the relative error used for clustering/grouping in ppm (parts per
  million)

- mzabs:

  the absolute error used for clustering/grouping

- minsamp:

  set the minimum number of samples in one bin

- minfrac:

  set the minimum fraction of each class in one bin

## Value

Returns a xcmsSet with slots groups and groupindex set.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),

## References

Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill and David
F. Grant\
*Alignment of high resolution mass spectra: development of a heuristic
approach for metabolomics*.\
Metabolomics, Vol. 2, No. 2, 75-83 (2006)

## Examples

``` r
if (FALSE) { # \dontrun{
library(MsDataHub)
mzMLfiles <- c(MsDataHub::HAM004_641fE_14.11.07..Exp1.extracted.mzML(),
               MsDataHub::HAM004_641fE_14.11.07..Exp2.extracted.mzML(),
               MsDataHub::HAM005_641fE_14.11.07..Exp1.extracted.mzML(),
               MsDataHub::HAM005_641fE_14.11.07..Exp2.extracted.mzML())

xs <- xcmsSet(method="MSW", files=mzMLfiles, scales=c(1,7),
              SNR.method='data.mean' , winSize.noise=500,
              peakThr=80000,  amp.Th=0.005)

xsg <- group(xs, method="mzClust")
} # }
```
