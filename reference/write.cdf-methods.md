# Save an xcmsRaw object to file

Write the raw data to a (simple) CDF file.

## Methods

- object = "xcmsRaw":

  `write.cdf(object, filename)`

## Arguments

- object:

  the `xcmsRaw` object

- filename:

  filename (may include full path) for the CDF file. Pipes or URLs are
  not allowed.

## Details

Currently the only application known to read the resulting file is XCMS.
Others, especially those which build on the AndiMS library, will refuse
to load the output.

## Value

None.

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md),
