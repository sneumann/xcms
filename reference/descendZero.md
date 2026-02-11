# Find start and end points of a peak

Decends down the sides of a data peak and finds either the points
greater than or equal to the zero intercept, the intercept with a given
value, or the bottom of the first valley on each side.

## Usage

``` r
descendZero(y, istart = which.max(y))
descendValue(y, value, istart = which.max(y))
descendMin(y, istart = which.max(y))
```

## Arguments

- y:

  numeric vector with values

- istart:

  starting point for descent

- value:

  numeric value to descend to

## Value

An integer vector of length 2 with the starting and ending indicies of
the peak start and end points.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

`descendValue`

## Examples

``` r
normdist <- dnorm(seq(-4, 4, .1)) - .1
xcms:::descendZero(normdist)
#> ilower iupper 
#>     25     57 
normdist[xcms:::descendZero(normdist)]
#> [1] 0.01092083 0.01092083
xcms:::descendValue(normdist, .15)
#> ilower iupper 
#>     32     50 
normdist[xcms:::descendValue(normdist, .15)]
#> [1] 0.1660852 0.1660852
xcms:::descendMin(normdist)
#> ilower iupper 
#>      1     81 
```
