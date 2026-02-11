# Allocate double, integer, or logical matricies

Allocate double, integer, or logical matricies in one step without
copying memory around.

## Usage

``` r
doubleMatrix(nrow = 0, ncol = 0)
integerMatrix(nrow = 0, ncol = 0)
logicalMatrix(nrow = 0, ncol = 0)
```

## Arguments

- nrow:

  number of matrix rows

- ncol:

  number of matrix columns

## Value

Matrix of double, integer, or logical values. Memory is not zeroed.

## Author

Colin A. Smith, <csmith@scripps.edu>
