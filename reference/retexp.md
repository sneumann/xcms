# Set retention time window to a specified width

Expands (or contracts) the retention time window in each row of a matrix
as defined by the `retmin` and `retmax` columns.

## Usage

``` r
retexp(peakrange, width = 200)
```

## Arguments

- peakrange:

  maxtrix with columns `retmin` and `retmax`

- width:

  new width for the window

## Value

The altered matrix.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`getEIC`](https://sneumann.github.io/xcms/reference/getEIC-methods.md)
