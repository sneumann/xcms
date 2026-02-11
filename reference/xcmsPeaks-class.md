# A matrix of peaks

A matrix of peak information. The actual columns depend on how it is
generated (i.e. the
[`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
method).

## Objects from the Class

Objects can be created by calls of the form `new("xcmsPeaks", ...)`.

## Slots

- `.Data`::

  The matrix holding the peak information

## Extends

Class `"`[`matrix`](https://rdrr.io/r/methods/StructureClasses.html)`"`,
from data part. Class
`"`[`array`](https://rdrr.io/r/methods/StructureClasses.html)`"`, by
class "matrix", distance 2. Class
`"`[`structure`](https://rdrr.io/r/methods/StructureClasses.html)`"`, by
class "matrix", distance 3. Class
`"`[`vector`](https://rdrr.io/r/methods/BasicClasses.html)`"`, by class
"matrix", distance 4, with explicit coerce.

## Methods

None yet. Some utilities for working with peak data would be nice.

## Author

Michael Lawrence

## See also

[`findPeaks`](https://sneumann.github.io/xcms/reference/findPeaks-methods.md)
for detecting peaks in an
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md).
