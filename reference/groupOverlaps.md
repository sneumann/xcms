# Group overlapping ranges

`groupOverlaps` identifies overlapping ranges in the input data and
groups them by returning their indices in `xmin` `xmax`.

## Usage

``` r
groupOverlaps(xmin, xmax)
```

## Arguments

- xmin:

  `numeric` (same length than `xmax`) with the lower boundary of the
  range.

- xmax:

  `numeric` (same length than `xmin`) with the upper boundary of the
  range.

## Value

`list` with the indices of grouped elements.

## Author

Johannes Rainer

## Examples

``` r

x <- c(2, 12, 34.2, 12.4)
y <- c(3, 16, 35, 36)

groupOverlaps(x, y)
#> [[1]]
#> [1] 1
#> 
#> [[2]]
#> [1] 2 3 4
#> 
```
