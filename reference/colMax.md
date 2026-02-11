# Find row and column maximum values

Find row and column maximum values for numeric arrays.

## Usage

``` r
colMax(x, na.rm = FALSE, dims = 1)
rowMax(x, na.rm = FALSE, dims = 1)
which.colMax(x, na.rm = FALSE, dims = 1)
which.rowMax(x, na.rm = FALSE, dims = 1)
```

## Arguments

- x:

  an array of two or more dimensions, containing numeric values

- na.rm:

  logical. Should missing values (including 'NaN') be omitted from the
  calculations? (not currently implemented)

- dims:

  Which dimensions are regarded as "rows" or "columns" to maximize. For
  `rowMax`, the maximum is over dimensions `dims+1, ...`; for `colMax`
  it is over dimensions `1:dims`.

## Details

These functions are designed to act like the `colSums` series of
functions except that they only currently handle real arrays and will
not remove `NA` values.

## Value

A numeric array of suitable size, or a vector if the result is
one-dimensional. The `dimnames` (or `names` for a vector result) are
taken from the original array.

For the `which.*` functions, an integer array of suitable size, or a
vector if the result is one-dimensional. The indecies returned are for
accessing `x` one-dimensionally (i.e. `x[index]`). For `which.colMax()`,
the actual row indecies my be determined using
`(which.colMax(x)-1) %% nrow(x) + 1`. For `which.rowMax()`, the actual
column indecies may be determined using `ceiling(rowMax(x)/nrow(x))`.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`colSums`](https://rdrr.io/r/base/colSums.html)
