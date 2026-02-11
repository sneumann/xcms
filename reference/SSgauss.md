# Gaussian Model

This `selfStart` model evalueates the Gaussian model and its gradient.
It has an `initial` attribute that will evalueate the inital estimates
of the parameters `mu`, `sigma`, and `h`.

## Usage

``` r
SSgauss(x, mu, sigma, h)
```

## Arguments

- x:

  a numeric vector of values at which to evaluate the model

- mu:

  mean of the distribution function

- sigma:

  standard deviation of the distribution fuction

- h:

  height of the distribution function

## Details

Initial values for `mu` and `h` are chosen from the maximal value of
`x`. The initial value for `sigma` is determined from the area under `x`
divided by `h*sqrt(2*pi)`.

## Value

A numeric vector of the same length as `x`. It is the value of the
expression `h*exp(-(x-mu)^2/(2*sigma^2)`, which is a modified gaussian
function where the maximum height is treated as a separate parameter not
dependent on `sigma`. If arguments `mu`, `sigma`, and `h` are names of
objects, the gradient matrix with respect to these names is attached as
an attribute named `gradient`.

## Author

Colin A. Smith, <csmith@scripps.edu>

## See also

[`nls`](https://rdrr.io/r/stats/nls.html),
[`selfStart`](https://rdrr.io/r/stats/selfStart.html)
