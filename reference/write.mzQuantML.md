# Save an xcmsSet object to an PSI mzQuantML file

Export in XML data formats: Write the processed data in an xcmsSet to
mzQuantML.

## Methods

- object = "xcmsSet":

  `write.mzQuantML(object, filename)`

## Arguments

- object:

  the `xcmsRaw` or `xcmsSet` object

- filename:

  filename (may include full path) for the output file. Pipes or URLs
  are not allowed.

## Details

The write.mzQuantML() function will write a (grouped) xcmsSet into the
PSI standard format mzQuantML, see <http://www.psidev.info/mzquantml>

## Value

None.

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`xcmsSet`](https://sneumann.github.io/xcms/reference/xcmsSet.md),
[`verify.mzQuantML`](https://sneumann.github.io/xcms/reference/verify.mzQuantML.md),
