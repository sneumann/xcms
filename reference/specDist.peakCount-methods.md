# a Distance function based on matching peaks

This method calculates the distance of two sets of peaks by just
returning the number of matching peaks (m/z-values).

## Usage

``` r
specDist.peakCount(peakTable1, peakTable2, mzabs=0.001, mzppm=10, symmetric=FALSE)
```

## Methods

- peakTable1 = "matrix", peakTable2 = "matrix":

  ` specDist.peakCount(peakTable1, peakTable2, mzppm=10,symmetric=FALSE ) `

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

## Author

Joachim Kutzera, <jkutzer@ipb-halle.de>
