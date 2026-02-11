# Group peaks from different samples together

Group peaks together across samples by creating a master peak list and
assigning corresponding peaks from all samples. It is inspired by the
alignment algorithm of mzMine. For further details check
<http://mzmine.sourceforge.net/> and

Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for processing and
visualization of mass spectrometry based molecular profile data.
Bioinformatics (Oxford, England) 2006, 22:634?636.

Currently, there is no equivalent to `minfrac` or `minsamp`.

## Methods

- object = "xcmsSet":

  ` group(object, mzVsRTbalance=10, mzCheck=0.2, rtCheck=15, kNN=10) `

## Arguments

- object:

  the `xcmsSet` object

- mzVsRTbalance:

  Multiplicator for mz value before calculating the (euclidean) distance
  between two peaks.

- mzCheck:

  Maximum tolerated distance for mz.

- rtCheck:

  Maximum tolerated distance for RT.

- kNN:

  Number of nearest Neighbours to check

## Value

An `xcmsSet` object with peak group assignments and statistics.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`group.density`](https://sneumann.github.io/xcms/reference/group.density.md)
and
[`group.mzClust`](https://sneumann.github.io/xcms/reference/group.mzClust.md)

## Examples

``` r
if (FALSE) library(xcms)
    library(faahKO)
    ## These files do not have this problem to correct for
    ## but just for an example
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)

    xset<-xcmsSet(cdffiles)

    gxset<-group(xset, method="nearest")
#> Processing sample number 8 ... 
#> OK
#> Processing sample number 1 ... 
#> OK
#> Processing sample number 7 ... 
#> OK
#> Processing sample number 3 ... 
#> OK
#> Processing sample number 12 ... 
#> OK
#> Processing sample number 9 ... 
#> OK
#> Processing sample number 6 ... 
#> OK
#> Processing sample number 5 ... 
#> OK
#> Processing sample number 4 ... 
#> OK
#> Processing sample number 11 ... 
#> OK
#> Processing sample number 10 ... 
#> OK
    nrow(gxset@groups) == 1096 ## the number of features before minFrac
#> [1] FALSE

    post.minFrac<-function(object, minFrac=0.5){
        ix.minFrac<-sapply(1:length(unique(sampclass(object))),
                           function(x, object, mf){
                               meta<-groups(object)
                               minFrac.idx<-numeric(length=nrow(meta))
                               idx<-which(
                                   meta[,levels(sampclass(object))[x]] >=
                                   mf*length(which(levels(sampclass(object))[x]
                                                   == sampclass(object)) ))
                               minFrac.idx[idx]<-1
                               return(minFrac.idx)
                           }, object, minFrac)
        ix.minFrac<-as.logical(apply(ix.minFrac, 1, sum))
        ix<-which(ix.minFrac == TRUE)
        return(ix)
    }

    ## using the above function we can get a post processing minFrac
    idx<-post.minFrac(gxset)

    gxset.post<-gxset ## copy the xcmsSet object
    gxset.post@groupidx<-gxset@groupidx[idx]
    gxset.post@groups<-gxset@groups[idx,]

    nrow(gxset.post@groups) == 465 ## number of features after minFrac
#> [1] FALSE
 # \dontrun{}
```
