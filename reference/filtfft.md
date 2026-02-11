# Apply an convolution filter using an FFT

Expands a vector to the length of the filter and then convolutes it
using two successive FFTs.

## Usage

``` r
filtfft(y, filt)
```

## Arguments

- y:

  numeric vector of data to be filtered

- filt:

  filter with length `nextn(length(y))`

## Value

A numeric vector the same length as y.

## Author

Colin A. Smith, <csmith@scripps.edu>
