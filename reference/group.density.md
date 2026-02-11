# Group peaks from different samples together

Group peaks together across samples using overlapping m/z bins and
calculation of smoothed peak distributions in chromatographic time.

## Methods

- object = "xcmsSet":

  ` group(object, bw = 30, minfrac = 0.5, minsamp = 1, mzwid = 0.25, max = 50, sleep = 0) `

## Arguments

- object:

  the `xcmsSet` object

- minfrac:

  minimum fraction of samples necessary in at least one of the sample
  groups for it to be a valid group

- minsamp:

  minimum number of samples necessary in at least one of the sample
  groups for it to be a valid group

- bw:

  bandwidth (standard deviation or half width at half maximum) of
  gaussian smoothing kernel to apply to the peak density chromatogram

- mzwid:

  width of overlapping m/z slices to use for creating peak density
  chromatograms and grouping peaks across samples

- max:

  maximum number of groups to identify in a single m/z slice

- sleep:

  seconds to pause between plotting successive steps of the peak
  grouping algorithm. peaks are plotted as points showing relative
  intensity. identified groups are flanked by dotted vertical lines.

## Value

An `xcmsSet` object with peak group assignments and statistics.

## See also

[`do_groupChromPeaks_density`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.md)
for the core API function performing the analysis.
[`xcmsSet-class`](https://sneumann.github.io/xcms/reference/xcmsSet-class.md),
[`density`](https://rdrr.io/r/stats/density.html)
