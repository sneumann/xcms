# a Distance function based on matching peaks

This method calculates the distance of two sets of peaks.

## Usage

``` r
specDist.meanMZmatch(peakTable1, peakTable2, matchdist=1, matchrate=1,
                     mzabs=0.001, mzppm=10, symmetric=TRUE)
```

## Methods

- peakTable1 = "matrix", peakTable2 = "matrix":

  ` specDist.meanMZmatch(peakTable1, peakTable2, matchdist=1, matchrate=1, mzabs=0.001, mzppm=10, symmetric=TRUE) `

## Details

The result of the calculation is a weighted sum of two values. Value one
is the mean absolute difference of the matching peaks, value two is the
relation of matching peaks and non matching peaks. if no distance is
calculated (for example because no matching peaks were found) the
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

- matchdist:

  the weight for value one (see details)

- matchrate:

  the weight for value two

## Author

Joachim Kutzera, <jkutzer@ipb-halle.de>
