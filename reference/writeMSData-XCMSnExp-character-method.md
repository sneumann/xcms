# Export MS data to mzML/mzXML files

`writeMSData` exports mass spectrometry data in mzML or mzXML format. If
adjusted retention times are present, these are used as retention time
of the exported spectra.

## Usage

``` r
# S4 method for class 'XCMSnExp,character'
writeMSData(
  object,
  file,
  outformat = c("mzml", "mzxml"),
  copy = FALSE,
  software_processing = NULL,
  ...
)
```

## Arguments

- object:

  [XCMSnExp](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
  object with the mass spectrometry data.

- file:

  `character` with the file name(s). The length of this parameter has to
  match the number of files/samples of `object`.

- outformat:

  `character(1)` defining the format of the output files ( either
  `"mzml"` or `"mzxml"`).

- copy:

  `logical(1)` if metadata (data processing, software used, original
  file names etc) should be copied from the original files.

- software_processing:

  optionally provide specific data processing steps. See documentation
  of the `software_processing` parameter of
  [`MSnbase::writeMSData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html).

- ...:

  Additional parameters to pass down to the
  [`MSnbase::writeMSData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
  function in the `MSnbase` package, such as `outformat` to specify the
  output format (`"mzml"` or `"mzxml"`) or `copy` to specify whether
  general information from the original MS data files (such as data
  processing, software etc) should be copied to the new files.

## See also

[`MSnbase::writeMSData()`](https://rdrr.io/pkg/ProtGenerics/man/protgenerics.html)
function in the `MSnbase` package.

## Author

Johannes Rainer
