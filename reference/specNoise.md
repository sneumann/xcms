# Calculate noise for a sparse continuum mass spectrum

Given a sparse continuum mass spectrum, determine regions where no
signal is present, substituting half of the minimum intensity for those
regions. Calculate the noise level as the weighted mean of the regions
with signal and the regions without signal. If there is only one raw
peak, return zero.

## Usage

``` r
specNoise(spec, gap = quantile(diff(spec[, "mz"]), 0.9))
```

## Arguments

- spec:

  matrix with named columns `mz` and `intensity`

- gap:

  threshold above which to data points are considerd to be separated by
  a blank region and not bridged by an interpolating line

## Details

The default gap value is determined from the 90th percentile of the
pair-wise differences between adjacent mass values.

## Value

A numeric noise level

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`getSpec`](https://sneumann.github.io/xcms/reference/getSpec-methods.md),
[`specPeaks`](https://sneumann.github.io/xcms/reference/specPeaks.md)
