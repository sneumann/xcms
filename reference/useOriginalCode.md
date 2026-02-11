# Enable usage of old xcms code

This function allows to enable the usage of old, partially deprecated
code from xcms by setting a corresponding global option. See details for
functions affected.

## Usage

``` r
useOriginalCode(x)
```

## Arguments

- x:

  `logical(1)` to specify whether or not original old code should be
  used in corresponding functions. If not provided the function simply
  returns the value of the global option.

## Value

`logical(1)` indicating whether old code is being used.

## Details

The functions/methods that are affected by this option are:

- [do_findChromPeaks_matchedFilter](https://sneumann.github.io/xcms/reference/do_findChromPeaks_matchedFilter.md):
  use the original code that iteratively creates a subset of the binned
  (profile) matrix. This is helpful for computers with limited memory or
  matchedFilter settings with a very small bin size.

- [getPeaks](https://sneumann.github.io/xcms/reference/getPeaks-methods.md)

## Note

For parallel processing using the SOCKS method (e.g. by
[`BiocParallel::SnowParam()`](https://rdrr.io/pkg/BiocParallel/man/SnowParam-class.html)
on Windows computers) this option might not be passed to the individual
R processes performing the calculations. In such cases it is suggested
to specify the option manually and system-wide by adding the line
`options(XCMSuseOriginalCode = TRUE)` in a file called *.Rprofile* in
the folder in which new R processes are started (usually the user's home
directory; to ensure that the option is correctly read add a new line to
the file too). See also [Startup](https://rdrr.io/r/base/Startup.html)
from the base R documentation on how to specify system-wide options for
R.

Usage of old code is strongly dicouraged. This function is thought to be
used mainly in the transition phase from xcms to xcms version 3.

## Author

Johannes Rainer
