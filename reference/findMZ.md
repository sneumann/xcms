# Find fragment ions in xcmsFragment objects

This is a method to find a fragment mass with a ppm window in a
xcmsFragment object

## Usage

``` r
findMZ(object, find, ppmE=25, print=TRUE)
```

## Arguments

- object:

  xcmsFragment object type

- find:

  The fragment ion to be found

- ppmE:

  the ppm error window for searching

- print:

  If we should print a nice little report

## Details

The method simply searches for a given fragment ion in an xcmsFragment
object type given a certain ppm error window

## Value

A data frame with the following columns:

- PrecursorMz:

  The precursor m/z of the fragment

- MSnParentPeakID:

  An index ID of the location of the precursor peak in the xcmsFragment
  object

- msLevel:

  The level of the found fragment ion

- rt:

  the Retention time of the found ion

- mz:

  the actual m/z of the found fragment ion

- intensity:

  The intensity of the fragment ion

- sample:

  Which sample the fragment ion came from

- GroupPeakMSn:

  an ID if the peaks were grouped by an xcmsSet grouping

- CollisionEnergy:

  The collision energy of the precursor scan

## References

H. Paul Benton, D.M. Wong, S.A.Strauger, G. Siuzdak "XC\\MS^2\\"
Analytical Chemistry 2008

## See also

[`findneutral`](https://sneumann.github.io/xcms/reference/findneutral.md),

## Author

H. Paul Benton, <hpaul.beonton08@imperial.ac.uk>
