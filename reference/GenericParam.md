# Generic parameter class

The `GenericParam` class allows to store generic parameter information
such as the name of the function that was/has to be called (slot `fun`)
and its arguments (slot `args`). This object is used to track the
process history of the data processings of an
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object. This is in contrast to e.g. the
[`CentWaveParam()`](https://sneumann.github.io/xcms/reference/findChromPeaks-centWave.md)
object that is passed to the actual processing method.

## Usage

``` r
GenericParam(fun = character(), args = list())
```

## Arguments

- fun:

  `character` representing the name of the function.

- args:

  `list` (ideally named) with the arguments to the function.

## Value

The `GenericParam()` function returns a `GenericParam` object.

## Slots

- `fun`:

  `character` specifying the function name.

- `args`:

  `list` (ideally named) with the arguments to the function.

## See also

[`processHistory()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for how to access the process history of an
[XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object.

## Author

Johannes Rainer

## Examples

``` r
prm <- GenericParam(fun = "mean")

prm <- GenericParam(fun = "mean", args = list(na.rm = TRUE))
```
