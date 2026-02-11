# Core API function for peak density based chromatographic peak grouping

The `do_groupChromPeaks_density` function performs chromatographic peak
grouping based on the density (distribution) of peaks, found in
different samples, along the retention time axis in slices of
overlapping m/z ranges. By default (with parameter `ppm = 0`) these m/z
ranges have all the same (constant) size (depending on parameter
`binSize`). For values of `ppm` larger than 0 the m/z bins (ranges or
slices) will have increasing sizes depending on the m/z value. This
better models the m/z-dependent measurement error/precision seen on some
MS instruments.

## Usage

``` r
do_groupChromPeaks_density(
  peaks,
  sampleGroups,
  bw = 30,
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.25,
  maxFeatures = 50,
  sleep = 0,
  index = seq_len(nrow(peaks)),
  ppm = 0
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

- bw:

  For `PeakDensityParam`: `numeric(1)` defining the bandwidth (standard
  deviation ot the smoothing kernel) to be used. This argument is passed
  to the \[stats::density() method.

- minFraction:

  For `PeakDensityParam`: `numeric(1)` defining the minimum fraction of
  samples in at least one sample group in which the peaks have to be
  present to be considered as a peak group (feature).

- minSamples:

  For `PeakDensityParam`: `numeric(1)` with the minimum number of
  samples in at least one sample group in which the peaks have to be
  detected to be considered a peak group (feature).

- binSize:

  For `PeakDensityParam`: `numeric(1)` defining the size of the
  overlapping slices in m/z dimension.

- maxFeatures:

  For `PeakDensityParam`: `numeric(1)` with the maximum number of peak
  groups to be identified in a single mz slice.

- sleep:

  `numeric(1)` defining the time to *sleep* between iterations and plot
  the result from the current iteration.

- index:

  An optional `integer` providing the indices of the peaks in the
  original peak matrix.

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

## Value

A `data.frame`, each row representing a (mz-rt) feature (i.e. a peak
group) with columns:

- `"mzmed"`: median of the peaks' apex mz values.

- `"mzmin"`: smallest mz value of all peaks' apex within the feature.

- `"mzmax"`:largest mz value of all peaks' apex within the feature.

- `"rtmed"`: the median of the peaks' retention times.

- `"rtmin"`: the smallest retention time of the peaks in the group.

- `"rtmax"`: the largest retention time of the peaks in the group.

- `"npeaks"`: the total number of peaks assigned to the feature.

- `"peakidx"`: a `list` with the indices of all peaks in a feature in
  the `peaks` input matrix.

Note that this number can be larger than the total number of samples,
since multiple peaks from the same sample could be assigned to a
feature.

## Details

For overlapping slices along the mz dimension, the function calculates
the density distribution of identified peaks along the retention time
axis and groups peaks from the same or different samples that are close
to each other. See (Smith 2006) for more details.

## Note

The default settings might not be appropriate for all LC/GC-MS setups,
especially the `bw` and `binSize` parameter should be adjusted
accordingly.

## References

Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
Anal. Chem. 2006, 78:779-787. doi:
[10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

## See also

Other core peak grouping algorithms:
[`do_groupChromPeaks_nearest()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_nearest.md),
[`do_groupPeaks_mzClust()`](https://sneumann.github.io/xcms/reference/do_groupPeaks_mzClust.md)

## Author

Colin Smith, Johannes Rainer

## Examples

``` r
## Load the test file
library(xcms)
library(MsExperiment)
faahko_sub <- loadXcmsData("faahko_sub2")

## Disable parallel processing for this example
register(SerialParam())

## Extract the matrix with the identified peaks from the xcmsSet:
pks <- chromPeaks(faahko_sub)

## Perform the peak grouping with default settings:
res <- do_groupChromPeaks_density(pks, sampleGroups = rep(1, 3))

## The feature definitions:
head(res)
#>   mzmed mzmin mzmax    rtmed    rtmin    rtmax npeaks 1      peakidx
#> 1 279.0 279.0 279.0 2787.765 2787.765 2787.766      2 2      11, 199
#> 2 286.2 286.2 286.2 3254.904 3250.992 3258.815      2 2     115, 205
#> 3 300.2 300.2 300.2 3387.143 3379.317 3390.271      4 3 35, 125,....
#> 4 301.0 301.0 301.0 2787.766 2786.200 2792.459      3 3  10, 97, 198
#> 5 305.1 305.1 305.1 2994.338 2994.338 2994.339      2 2      15, 203
#> 6 305.1 305.1 305.1 2923.917 2923.916 2923.917      2 2      14, 202
```
