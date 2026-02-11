# Create report of aligned peak intensities

Create a report showing all aligned peaks.

## Methods

- object = "xcmsSet":

  `peakTable(object, filebase = character(), ...)`

## Arguments

- object:

  the `xcmsSet` object

- filebase:

  base file name to save report, `.tsv` file and `_eic` will be appended
  to this name for the tabular report and EIC directory, respectively.
  if blank nothing will be saved

- ...:

  arguments passed down to
  [`groupval`](https://sneumann.github.io/xcms/reference/groupval-methods.md),
  which provides the actual intensities.

## Details

This method handles creation of summary reports similar to
[`diffreport`](https://sneumann.github.io/xcms/reference/diffreport-methods.md).
It returns a summary report that can optionally be written out to a
tab-separated file.

If a base file name is provided, the report (see Value section) will be
saved to a tab separated file.

## Value

A data frame with the following columns:

- mz:

  median m/z of peaks in the group

- mzmin:

  minimum m/z of peaks in the group

- mzmax:

  maximum m/z of peaks in the group

- rt:

  median retention time of peaks in the group

- rtmin:

  minimum retention time of peaks in the group

- rtmax:

  maximum retention time of peaks in the group

- npeaks:

  number of peaks assigned to the group

- Sample Classes:

  number samples from each sample class represented in the group

- ...:

  one column for every sample class

- Sample Names:

  integrated intensity value for every sample

- ...:

  one column for every sample

## See also

[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),

## Examples
