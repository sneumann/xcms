# Apply a median filter to a matrix

For each element in a matix, replace it with the median of the values
around it.

## Usage

``` r
medianFilter(x, mrad, nrad)
```

## Arguments

- x:

  numeric matrix to median filter

- mrad:

  number of rows on either side of the value to use for median
  calculation

- nrad:

  number of rows on either side of the value to use for median
  calculation

## Value

A matrix whose values have been median filtered

## Author

Colin A. Smith, <csmith@scripps.edu>

## Examples

``` r
mat <- matrix(1:25, nrow=5)
mat
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]    1    6   11   16   21
#> [2,]    2    7   12   17   22
#> [3,]    3    8   13   18   23
#> [4,]    4    9   14   19   24
#> [5,]    5   10   15   20   25
medianFilter(mat, 1, 1)
#>      [,1] [,2] [,3] [,4] [,5]
#> [1,]  4.0  6.5 11.5 16.5 19.0
#> [2,]  4.5  7.0 12.0 17.0 19.5
#> [3,]  5.5  8.0 13.0 18.0 20.5
#> [4,]  6.5  9.0 14.0 19.0 21.5
#> [5,]  7.0  9.5 14.5 19.5 22.0
```
