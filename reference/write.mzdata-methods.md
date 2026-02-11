# Save an xcmsRaw object to a file

Write the raw data to a (simple) mzData file.

## Methods

- object = "xcmsRaw":

  `write.mzdata(object, filename)`

## Arguments

- object:

  the `xcmsRaw` object

- filename:

  filename (may include full path) for the mzData file. Pipes or URLs
  are not allowed.

## Details

This function will export a given xcmsRaw object to an mzData file. The
mzData file will contain a \<spectrumList\> containing the \<spectrum\>
with mass and intensity values in 32 bit precision. Other formats are
currently not supported. Any header information (e.g. additional
\<software\> information or \<cvParams\>) will be lost. Currently, also
any MSn information will not be stored.

## Value

None.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md),
