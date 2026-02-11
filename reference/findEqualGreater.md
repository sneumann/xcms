# Find values in sorted vectors

Find values in sorted vectors.

## Usage

``` r
findEqualGreater(x, value)
findEqualLess(x, value)
findEqualGreaterM(x, values)
findRange(x, values, NAOK = FALSE)
```

## Arguments

- x:

  numeric vector sorted in increasing order

- value:

  value to find in `x`

- values:

  numeric values to find in `x`

- NAOK:

  don't check for NA values in `x`

## Details

`findEqualGreater` finds the index of the first value in `x` that is
equal or greater than `value`. `findEqualLess` does same except that it
finds equal or less. `findEqualGreaterM` creates an index of a vector by
finding specified values. `findRange` locates the start and stop
indicides of a range of two `x` values.

The only things that save time at this point are `findeEqualGreaterM`
(when the length of values approaches the lenght of x) and `findRange`
(when `NAOK` is set to TRUE). They run in log(N) and N time,
respectively.

## Value

An integer vector with the position(s) of the values(s).

## Author

Colin A. Smith, <csmith@scripps.edu>
