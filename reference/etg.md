# Empirically Transformed Gaussian function

A general function for asymmetric chromatographic peaks.

## Usage

``` r
etg(x, H, t1, tt, k1, kt, lambda1, lambdat, alpha, beta)
```

## Arguments

- x:

  times to evaluate function at

- H:

  peak height

- t1:

  time of leading edge inflection point

- tt:

  time of trailing edge inflection point

- k1:

  leading edge parameter

- kt:

  trailing edge parameter

- lambda1:

  leading edge parameter

- lambdat:

  trailing edge parameter

- alpha:

  leading edge parameter

- beta:

  trailing edge parameter

## Value

The function evaluated at times `x`.

## References

Jianwei Li. Development and Evaluation of Flexible Empirical Peak
Functions for Processing Chromatographic Peaks. Anal. Chem., 69 (21),
4452-4462, 1997. <http://dx.doi.org/10.1021/ac970481d>

## Author

Colin A. Smith, <csmith@scripps.edu>
