# Generate p-values for a vector of t-statistics

Generate p-values for a vector of Welch's two-sample t-statistics based
on the t distribution.

## Usage

``` r
pval(X, classlabel, teststat)
```

## Arguments

- X:

  original data matrix

- classlabel:

  integer vector with classlabel

- teststat:

  numeric vector with Welch's two-sample t-statistics

## Value

A numeric vector of p-values.

## Author

Colin A. Smith, <csmith@scripps.edu>
