# Core API function for peak grouping using mzClust

The `do_groupPeaks_mzClust` function performs high resolution
correspondence on single spectra samples.

## Usage

``` r
do_groupPeaks_mzClust(
  peaks,
  sampleGroups,
  ppm = 20,
  absMz = 0,
  minFraction = 0.5,
  minSamples = 1
)
```

## Arguments

- peaks:

  A `matrix` or `data.frame` with the mz values and retention times of
  the identified chromatographic peaks in all samples of an experiment.
  Required columns are `"mz"`, `"rt"` and `"sample"`. The latter should
  contain `numeric` values representing the index of the sample in which
  the peak was found.

- sampleGroups:

  For `PeakDensityParam`: A vector of the same length than samples
  defining the sample group assignments (i.e. which samples belong to
  which sample group). This parameter is mandatory for
  `PeakDensityParam` and has to be defined also if there is no sample
  grouping in the experiment (in which case all samples should be
  assigned to the same group). Samples for which a `NA` is provided will
  not be considered in the feature definitions step. Providing `NA` for
  all blanks in an experiment will for example avoid features to be
  defined for signals (chrom peaks) present only in blank samples.

- ppm:

  For `MzClustParam`: `numeric(1)` representing the relative m/z error
  for the clustering/grouping (in parts per million). For
  `PeakDensityParam`: `numeric(1)` to define m/z-dependent, increasing
  m/z bin sizes. If `ppm = 0` (the default) m/z bins are defined by the
  sequence of values from the smallest to the larges m/z value with a
  constant bin size of `binSize`. For `ppm` \> 0 the size of each bin is
  increased in addition by the `ppm` of the (upper) m/z boundary of the
  bin. The maximal bin size (used for the largest m/z values) would then
  be `binSize` plus `ppm` parts-per-million of the largest m/z value of
  all peaks in the data set.

- absMz:

  For `NearestPeaksParam` and `MzClustParam`: `numeric(1)` maximum
  tolerated distance for m/z values.

- minFraction:

  For `PeakDensityParam`: `numeric(1)` defining the minimum fraction of
  samples in at least one sample group in which the peaks have to be
  present to be considered as a peak group (feature).

- minSamples:

  For `PeakDensityParam`: `numeric(1)` with the minimum number of
  samples in at least one sample group in which the peaks have to be
  detected to be considered a peak group (feature).

## Value

A `list` with elements `"featureDefinitions"` and `"peakIndex"`.
`"featureDefinitions"` is a `matrix`, each row representing an (mz-rt)
feature (i.e. peak group) with columns:

- `"mzmed"`: median of the peaks' apex mz values.

- `"mzmin"`: smallest mz value of all peaks' apex within the feature.

- `"mzmax"`: largest mz value of all peaks' apex within the feature.

- `"rtmed"`: always `-1`.

- `"rtmin"`: always `-1`.

- `"rtmax"`: always `-1`.

- `"npeaks"`: the total number of peaks assigned to the feature. Note
  that this number can be larger than the total number of samples, since
  multiple peaks from the same sample could be assigned to a group.

`"peakIndex"` is a `list` with the indices of all peaks in a peak group
in the `peaks` input matrix.

## References

Saira A. Kazmi, Samiran Ghosh, Dong-Guk Shin, Dennis W. Hill and David
F. Grant\
*Alignment of high resolution mass spectra: development of a heuristic
approach for metabolomics*.\
Metabolomics, Vol. 2, No. 2, 75-83 (2006)

## See also

Other core peak grouping algorithms:
[`do_groupChromPeaks_density()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.md),
[`do_groupChromPeaks_nearest()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_nearest.md)
