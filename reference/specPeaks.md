# Identify peaks in a sparse continuum mode spectrum

Given a spectrum, identify and list significant peaks as determined by
several criteria.

## Usage

``` r
specPeaks(spec, sn = 20, mzgap = 0.2)
```

## Arguments

- spec:

  matrix with named columns `mz` and `intensity`

- sn:

  minimum signal to noise ratio

- mzgap:

  minimal distance between adjacent peaks, with smaller peaks being
  excluded

## Details

Peaks must meet two criteria to be considered peaks: 1) Their s/n ratio
must exceed a certain threshold. 2) They must not be within a given
distance of any greater intensity peaks.

## Value

A matrix with columns:

- mz:

  m/z at maximum peak intensity

- intensity:

  maximum intensity of the peak

- fwhm:

  full width at half max of the peak

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`getSpec`](https://sneumann.github.io/xcms/reference/getSpec-methods.md),
[`specNoise`](https://sneumann.github.io/xcms/reference/specNoise.md)
