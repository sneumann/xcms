# Base class for loading raw data from a file

Data sources which read data from a file should inherit from this class.
The `xcms` package provides classes to read from `netCDF`, `mzData`,
`mzXML`, and `mzML` files using `xcmsFileSource`.

This class should be considered virtual and will not work if passed to
[`loadRaw-methods`](https://sneumann.github.io/xcms/reference/loadRaw-methods.md).
The reason it is not explicitly virtual is that there does not appear to
be a way for a class to be both virtual and have a data part (which lets
functions treat objects as if they were character strings).

This class validates that a file exists at the path given.

## Objects from the Class

`xcmsFileSource` objects should not be instantiated directly. Instead,
create subclasses and instantiate those.

## Slots

- `.Data`::

  Object of class `"character"`. File path of a file from which to read
  raw data as the object's data part

## Extends

Class `"`[`character`](https://rdrr.io/r/methods/BasicClasses.html)`"`,
from data part. Class
`"`[`xcmsSource`](https://sneumann.github.io/xcms/reference/xcmsSource-class.md)`"`,
directly.

## Methods

- `xcmsSource`:

  `signature(object = "character")`: Create an `xcmsFileSource` object
  referencing the given file name.

## Author

Daniel Hackney <dan@haxney.org>

## See also

[`xcmsSource`](https://sneumann.github.io/xcms/reference/xcmsSource-class.md)
