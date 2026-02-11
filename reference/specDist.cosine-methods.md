# a Distance function based on matching peaks

This method calculates the distance of two sets of peaks using the
cosine-distance.

## Usage

``` r
specDist.cosine(peakTable1, peakTable2, mzabs=0.001, mzppm=10, mzExp=0.6,
                intExp=3, nPdiff=2, nPmin=8, symmetric=FALSE)
```

## Methods

- peakTable1 = "matrix", peakTable2 = "matrix":

  ` specDist.cosine(peakTable1, peakTable2, mzabs = 0.001, mzppm = 10, mzExp = 0.6, intExp = 3, nPdiff = 2, nPmin = 8, symmetric = FALSE) `

## Details

The result is the cosine-distance of the product from weighted factors
of mz and intensity from matching peaks in the two peaktables. The
factors are calculated as wFact = mz^mzExp \* int^intExp. if no distance
is calculated (for example because no matching peaks were found) the
return-value is NA.

## Arguments

- peakTable1:

  a Matrix containing at least m/z-values, row must be called "mz"

- peakTable2:

  the matrix for the other mz-values

- mzabs:

  maximum absolute deviation for two matching peaks

- mzppm:

  relative deviations in ppm for two matching peaks

- symmetric:

  use symmetric pairwise m/z-matches only, or each match

- mzExp:

  the exponent used for mz

- intExp:

  the exponent used for intensity

- nPdiff:

  the maximum nrow-difference of the two peaktables

- nPmin:

  the minimum absolute sum of peaks from both praktables

## Author

Joachim Kutzera, <jkutzer@ipb-halle.de>
