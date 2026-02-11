# Virtual class for raw data sources

This virtual class provides an implementation-independent way to load
mass spectrometer data from various sources for use in an
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
object. Subclasses can be defined to enable data to be loaded from
user-specified sources. The virtual class
[`xcmsFileSource`](https://sneumann.github.io/xcms/reference/xcmsFileSource-class.md)
is included out of the box which contains a file name as a character
string.

When implementing child classes of `xcmsSource`, a corresponding
[`loadRaw-methods`](https://sneumann.github.io/xcms/reference/loadRaw-methods.md)
method must be provided which accepts the `xcmsSource` child class and
returns a list in the format described in
[`loadRaw-methods`](https://sneumann.github.io/xcms/reference/loadRaw-methods.md).

## Objects from the Class

A virtual Class: No objects may be created from it.

## Author

Daniel Hackney, <dan@haxney.org>

## See also

[`xcmsSource-methods`](https://sneumann.github.io/xcms/reference/xcmsSource-methods.md)
for creating `xcmsSource` objects in various ways.
