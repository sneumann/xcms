# Core API function for chromatic peak grouping using a nearest neighbor approach

The `do_groupChromPeaks_nearest` function groups peaks across samples by
creating a master peak list and assigning corresponding peaks from all
samples to each peak group (i.e. feature). The method is inspired by the
correspondence algorithm of mzMine (Katajamaa 2006).

## Usage

``` r
do_groupChromPeaks_nearest(
  peaks,
  sampleGroups,
  mzVsRtBalance = 10,
  absMz = 0.2,
  absRt = 15,
  kNN = 10
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

- mzVsRtBalance:

  For `NearestPeaksParam`: `numeric(1)` representing the factor by which
  m/z values are multiplied before calculating the (euclician) distance
  between two peaks.

- absMz:

  For `NearestPeaksParam` and `MzClustParam`: `numeric(1)` maximum
  tolerated distance for m/z values.

- absRt:

  For `NearestPeaksParam`: `numeric(1)` maximum tolerated distance for
  retention times.

- kNN:

  For `NearestPeaksParam`: `integer(1)` representing the number of
  nearest neighbors to check.

## Value

A `list` with elements `"featureDefinitions"` and `"peakIndex"`.
`"featureDefinitions"` is a `matrix`, each row representing an (mz-rt)
feature (i.e. peak group) with columns:

- `"mzmed"`: median of the peaks' apex mz values.

- `"mzmin"`: smallest mz value of all peaks' apex within the feature.

- `"mzmax"`:largest mz value of all peaks' apex within the feature.

- `"rtmed"`: the median of the peaks' retention times.

- `"rtmin"`: the smallest retention time of the peaks in the feature.

- `"rtmax"`: the largest retention time of the peaks in the feature.

- `"npeaks"`: the total number of peaks assigned to the feature.

`"peakIndex"` is a `list` with the indices of all peaks in a feature in
the `peaks` input matrix.

## References

Katajamaa M, Miettinen J, Oresic M: MZmine: Toolbox for processing and
visualization of mass spectrometry based molecular profile data.
Bioinformatics 2006, 22:634-636. doi:
[10.1093/bioinformatics/btk039](https://doi.org/10.1093/bioinformatics/btk039)

## See also

Other core peak grouping algorithms:
[`do_groupChromPeaks_density()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.md),
[`do_groupPeaks_mzClust()`](https://sneumann.github.io/xcms/reference/do_groupPeaks_mzClust.md)
