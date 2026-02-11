# Constructor for xcmsFragments objects which holds Tandem MS peaks

EXPERIMANTAL FEATURE

xcmsFragments is an object similar to xcmsSet, which holds peaks picked
(or collected) from one or several xcmsRaw objects.

There are still discussions going on about the exact API for MS\$^n\$
data, so this is likely to change in the future. The code is not yet
pipeline-ified.

## Usage

``` r
xcmsFragments(xs, ...)
```

## Arguments

- xs:

  A
  [`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md)
  object which contains picked ms1-peaks from one or several experiments

- ...:

  further arguments to the `collect` method

## Details

After running collect(xFragments,xSet) The peaktable of the
xcmsFragments includes the ms1Peaks from all experinemts stored in a
xcmsSet-object. Further it contains the relevant MSn-peaks from the
xcmsRaw-objects, which were created temporarily with the paths in
xcmsSet.

## Value

An `xcmsFragments` object.

## Author

Joachim Kutzera, Steffen Neumann, <sneumann@ipb-halle.de>

## See also

[`xcmsFragments-class`](https://sneumann.github.io/xcms/reference/xcmsFragments-class.md),
[`collect`](https://sneumann.github.io/xcms/reference/collect-methods.md)
