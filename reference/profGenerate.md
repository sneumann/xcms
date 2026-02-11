# Generation of profile data

Generates profile (binned) data in a given range from an indexed pair of
vectors.

## Usage

``` r
profBin(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
         param = list())
profBinLin(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinLinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
profBinLinBase(x, y, num, xstart = min(x), xend = max(x), param = list())
profBinLinBaseM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
                param = list())
profIntLin(x, y, num, xstart = min(x), xend = max(x), param = list())
profIntLinM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
profMaxIdx(x, y, num, xstart = min(x), xend = max(x), param = list())
profMaxIdxM(x, y, zidx, num, xstart = min(x), xend = max(x), NAOK = FALSE,
            param = list())
```

## Arguments

- x:

  numeric vector of value positions

- y:

  numeric vector of values to bin

- zidx:

  starting position of each new segment

- num:

  number of equally spaced x bins

- xstart:

  starting x value

- xend:

  ending x value

- NAOK:

  allow NA values (faster)

- param:

  parameters for profile generation

## Details

These functions take a vector of unequally spaced `y` values and
transform them into either a vector or matrix, depending on whether
there is an index or not. Each point in the vector or matrix represents
the data for the point centered at its corresponding `x` value, plus or
minus half the `x` step size (`xend-xstart/(num-1)`).

The `Bin` functions set each matrix or vector value to the maximal point
that gets binned into it.

The `BinLin` functions do the same except that they linearly interpolate
values into which nothing was binned.

The `BinLinBase` functions do the same except that they populate empty
parts of spectra with a base value. They take to two parameters: 1)
`baselevel`, the intensity level to fill in for empty parts of the
spectra. It defaluts to half of the minimum intensity. 2) `basespace`,
the m/z length after which the signal will drop to the base level.
Linear interpolation will be used between consecuitive data points
falling within `2*basespace` of eachother. It defaluts to 0.075.

The `IntLin` functions set each matrix or vector value to the integral
of the linearly interpolated data from plus to minus half the step size.

The `MaxIdx` functions work similarly to the `Bin` functions execpt that
the return the integer index of which x,y pair would be placed in a
particular cell.

## Note

There are some issues with the `profBinLin` method, see
<https://github.com/sneumann/xcms/issues/46> and
<https://github.com/sneumann/xcms/issues/49>. Thus it is suggested to
use the functions
[`binYonX`](https://sneumann.github.io/xcms/reference/binYonX.md) in
combination with
[`imputeLinInterpol`](https://sneumann.github.io/xcms/reference/imputeLinInterpol.md)
instead.

## Value

For `prof*`, a numeric vector of length `num`.

For `prof*M`, a matrix with dimensions `num` by `length(zidx)`.

For `MaxIdx`, the data type is integer, for all others it is double.

## Author

Colin A. Smith, <csmith@scripps.edu>

## Examples

``` r
if (FALSE) { # \dontrun{
    library(faahKO)
    cdfpath <- system.file("cdf", package = "faahKO")
    cdffiles <- list.files(cdfpath, recursive = TRUE, full.names = TRUE)
    xraw <- xcmsRaw(cdffiles[1])

    image(xraw) ## not how with intLin the intensity's blur
    profMethod(xraw) <- "bin"
    image(xraw) ## now with 'bin' there is no blurring good for centroid data
    ##try binlinbase for profile data
} # }
```
