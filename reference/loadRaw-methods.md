# Read binary data from a source

This function extracts the raw data which will be used an
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md)
object. Further processing of data is done in the
[`xcmsRaw`](https://sneumann.github.io/xcms/reference/xcmsRaw.md)
constructor.

## Methods

- `signature(object = "xcmsSource")`:

  Uses `loadRaw,xcmsSource-method` to extract raw data. Subclasses of
  [`xcmsSource`](https://sneumann.github.io/xcms/reference/xcmsSource-class.md)
  can provide different ways of fetching data.

## Arguments

- object:

  Specification of a data source (such as a file name or database query)

## Details

The implementing methods decide how to gather the data.

## Value

A list containing elements describing the data source. The `rt`,
`scanindex`, `tic`, and `acquisitionNum` components each have one entry
per scan. They are *parallel* in the sense that `rt[1]`, `scanindex[1]`,
and `acquisitionNum[1]` all refer to the same scan. The list containst
the following components:

- `rt`:

  Numeric vector with acquisition time (in seconds) for each scan

- `tic`:

  Numeric vector with Total Ion Count for each scan

- `scanindex`:

  Integer vector with starting positions of each scan in the `mz` and
  `intensity` components. It is an exclusive offset, so `scanindex[i]`
  is the offset in `mz` and `intensity` *before* the beginning of scan
  `i`. This means that the `mz` (respectively `intensity`) values for
  scan `i` would be from `scanindex[i] + 1` to `scanindex[i + 1]`

- `mz`:

  Concatenated vector of m/z values for all scans

- `intensity`:

  Concatenated vector of intensity values for all scans

## Author

Daniel Hackney, <dan@haxney.org>

## See also

[`xcmsRaw-class`](https://sneumann.github.io/xcms/reference/xcmsRaw-class.md),
[`xcmsSource`](https://sneumann.github.io/xcms/reference/xcmsSource-class.md)
