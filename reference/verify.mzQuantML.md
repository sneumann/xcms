# Verify an mzQuantML file

Export in XML data formats: verify the written data

## Usage

``` r
verify.mzQuantML(filename, xsdfilename)
```

## Arguments

- filename:

  filename (may include full path) for the output file. Pipes or URLs
  are not allowed.

- xsdfilename:

  Filename of the XSD to verify against (may include full path)

## Details

The verify.mzQuantML() function will verify an PSI standard format
mzQuantML document against the XSD schemda, see
<http://www.psidev.info/mzquantml>

## Value

None.

## See also

[`write.mzQuantML`](https://sneumann.github.io/xcms/reference/write.mzQuantML.md)
