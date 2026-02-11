# Align spectrum retention times across samples using peak groups found in most samples

The function performs retention time correction by assessing the
retention time deviation across all samples using peak groups (features)
containg chromatographic peaks present in most/all samples. The
retention time deviation for these features in each sample is described
by fitting either a polynomial (`smooth = "loess"`) or a linear
(`smooth = "linear"`) model to the data points. The models are
subsequently used to adjust the retention time for each spectrum in each
sample.

## Usage

``` r
do_adjustRtime_peakGroups(
  peaks,
  peakIndex,
  rtime = list(),
  minFraction = 0.9,
  extraPeaks = 1,
  smooth = c("loess", "linear"),
  span = 0.2,
  family = c("gaussian", "symmetric"),
  peakGroupsMatrix = matrix(ncol = 0, nrow = 0),
  subset = integer(),
  subsetAdjust = c("average", "previous")
)
```

## Arguments

- peaks:

  a `matrix` or `data.frame` with the identified chromatographic peaks
  in the samples.

- peakIndex:

  a `list` of indices that provides the grouping information of the
  chromatographic peaks (across and within samples).

- rtime:

  a `list` of `numeric` vectors with the retention times per
  file/sample.

- minFraction:

  For `PeakGroupsParam`: `numeric(1)` between 0 and 1 defining the
  minimum required proportion of samples in which peaks for the peak
  group were identified. Peak groups passing this criteria will be
  aligned across samples and retention times of individual spectra will
  be adjusted based on this alignment. For `minFraction = 1` the peak
  group has to contain peaks in all samples of the experiment. Note that
  if `subset` is provided, the specified fraction is relative to the
  defined subset of samples and not to the total number of samples
  within the experiment (i.e., a peak has to be present in the specified
  proportion of subset samples).

- extraPeaks:

  For `PeakGroupsParam`: `numeric(1)` defining the maximal number of
  additional peaks for all samples to be assigned to a peak group
  (feature) for retention time correction. For a data set with 6
  samples, `extraPeaks = 1` uses all peak groups with a total peak count
  `<= 6 + 1`. The total peak count is the total number of peaks being
  assigned to a peak group and considers also multiple peaks within a
  sample that are assigned to the group. This parameter is ignored for
  [`adjustRtime()`](https://sneumann.github.io/xcms/reference/adjustRtime.md)
  on an
  [`XcmsExperimentHdf5()`](https://sneumann.github.io/xcms/reference/XcmsExperimentHdf5.md).

- smooth:

  For `PeakGroupsParam`: `character(1)` defining the function to be used
  to interpolate corrected retention times for all peak groups. Can be
  either `"loess"` or `"linear"`.

- span:

  For `PeakGroupsParam`: `numeric(1)` defining the degree of smoothing
  (if `smooth = "loess"`). This parameter is passed to the internal call
  to [`stats::loess()`](https://rdrr.io/r/stats/loess.html).

- family:

  For `PeakGroupsParam`: `character(1)` defining the method for loess
  smoothing. Allowed values are `"gaussian"` and `"symmetric"`. See
  [`stats::loess()`](https://rdrr.io/r/stats/loess.html) for more
  information.

- peakGroupsMatrix:

  optional `matrix` of (raw) retention times for peak groups on which
  the alignment should be performed. Each column represents a sample,
  each row a feature/peak group. If not provided, this matrix will be
  determined depending on parameters `minFraction` and `extraPeaks`. If
  provided, `minFraction` and `extraPeaks` will be ignored.

- subset:

  For `ObiwarpParam` and `PeakGroupsParam`: `integer` with the indices
  of samples within the experiment on which the alignment models should
  be estimated. Samples not part of the subset are adjusted based on the
  closest subset sample. See *Subset-based alignment* section for
  details.

- subsetAdjust:

  For `ObiwarpParam` and `PeakGroupsParam`: `character(1)` specifying
  the method with which non-subset samples should be adjusted. Supported
  options are `"previous"` and `"average"` (default). See *Subset-based
  alignment* section for details.

## Value

A `list` with `numeric` vectors with the adjusted retention times
grouped by sample.

## Details

The alignment bases on the presence of compounds that can be found in
all/most samples of an experiment. The retention times of individual
spectra are then adjusted based on the alignment of the features
corresponding to these *house keeping compounds*. The parameters
`minFraction` and `extraPeaks` can be used to fine tune which features
should be used for the alignment (i.e. which features most likely
correspond to the above mentioned house keeping compounds).

Parameter `subset` allows to define a subset of samples within the
experiment that should be aligned. All samples not being part of the
subset will be aligned based on the adjustment of the closest sample
within the subset. This allows to e.g. exclude blank samples from the
alignment process with their retention times being still adjusted based
on the alignment results of the *real* samples.

## Note

The method ensures that returned adjusted retention times are
increasingly ordered, just as the raw retention times.

## References

Colin A. Smith, Elizabeth J. Want, Grace O'Maille, Ruben Abagyan and
Gary Siuzdak. "XCMS: Processing Mass Spectrometry Data for Metabolite
Profiling Using Nonlinear Peak Alignment, Matching, and Identification"
*Anal. Chem.* 2006, 78:779-787.

## Author

Colin Smith, Johannes Rainer
