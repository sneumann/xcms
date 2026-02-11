# Class xcmsEIC, a class for multi-sample extracted ion chromatograms

This class is used to store and plot parallel extracted ion
chromatograms from multiple sample files. It integrates with the
`xcmsSet` class to display peak area integrated during peak
identification or fill-in.

## Objects from the Class

Objects can be created with the
[`getEIC`](https://sneumann.github.io/xcms/reference/getEIC-methods.md)
method of the `xcmsSet` class. Objects can also be created by calls of
the form `new("xcmsEIC", ...)`.

## Slots

- `eic`::

  list containing named entries for every sample. for each entry, a list
  of two column EIC matricies with retention time and intensity

- `mzrange`::

  two column matrix containing starting and ending m/z for each EIC

- `rtrange`::

  two column matrix containing starting and ending time for each EIC

- `rt`::

  either `"raw"` or `"corrected"` to specify retention times contained
  in the object

- `groupnames`::

  group names from `xcmsSet` object used to generate EICs

## Methods

- [groupnames](https://sneumann.github.io/xcms/reference/groupnames-methods.md):

  `signature(object = "xcmsEIC")`: get `groupnames` slot

- mzrange:

  `signature(object = "xcmsEIC")`: get `mzrange` slot

- [plot](https://sneumann.github.io/xcms/reference/plot.xcmsEIC.md):

  `signature(x = "xcmsEIC")`: plot the extracted ion chromatograms

- rtrange:

  `signature(object = "xcmsEIC")`: get `rtrange` slot

- [sampnames](https://sneumann.github.io/xcms/reference/sampnames-methods.md):

  `signature(object = "xcmsEIC")`: get sample names

## Author

Colin A. Smith, <csmith@scripps.edu>

## Note

No notes yet.

## See also

[`getEIC`](https://sneumann.github.io/xcms/reference/getEIC-methods.md)
