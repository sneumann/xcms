# Correspondence: group chromatographic peaks across samples

The `groupChromPeaks` method performs a correspondence analysis i.e., it
groups chromatographic peaks across samples to define the LC-MS
*features*. The correspondence algorithm can be selected, and
configured, using the `param` argument. See documentation of
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
and
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
for information on how to access and extract correspondence results.

The correspondence analysis can be performed on chromatographic peaks of
any MS level (if present and if chromatographic peak detection has been
performed for that MS level) defining features combining these peaks.
The MS level can be selected with the parameter `msLevel`. By default,
calling `groupChromPeaks` will remove any previous correspondence
results. This can be disabled with `add = TRUE`, which will add newly
defined features to already present feature definitions.

Supported `param` objects are:

- `PeakDensityParam`: correspondence using the *peak density* method
  (Smith 2006) that groups chromatographic peaks along the retention
  time axis within slices of (partially overlapping) m/z ranges. By
  default, these m/z ranges (bins) have a constant size. By setting
  `ppm` to a value larger than 0, m/z dependent bin sizes can be used
  instead (better representing the m/z dependent measurement error of
  some MS instruments). All peaks (from the same or from different
  samples) with their apex position being close on the retention time
  axis are grouped into a LC-MS feature. Only samples with non-missing
  sample group assignment (i.e. for which the value provided with
  parameter `sampleGroups` is different than `NA`) are considered and
  counted for the feature definition. This allows to exclude certain
  samples or groups (e.g. blanks) from the feature definition avoiding
  thus features with only detected peaks in these. Note that this
  affects only the **definition** of **new** features. Chromatographic
  peaks in these samples will still be assigned to features which were
  defined based on the other samples. See in addition
  [`do_groupChromPeaks_density()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_density.md)
  for the core API function.

- `NearestPeaksParam`: performs peak grouping based on the proximity of
  chromatographic peaks from different samples in the m/z - rt space
  similar to the correspondence method of *mzMine* (Katajamaa 2006). The
  method creates first a *master peak list* consisting of all
  chromatographic peaks from the sample with the most detected peaks and
  iteratively calculates distances to peaks from the sample with the
  next most number of peaks grouping peaks together if their *distance*
  is smaller than the provided thresholds. See in addition
  [`do_groupChromPeaks_nearest()`](https://sneumann.github.io/xcms/reference/do_groupChromPeaks_nearest.md)
  for the core API function.

- `MzClustParam`: performs high resolution peak grouping for **single
  spectrum** metabolomics data (Kazmi 2006). This method should **only**
  be used for such data as the retention time is not considered in the
  correspondence analysis. See in addition
  [`do_groupPeaks_mzClust()`](https://sneumann.github.io/xcms/reference/do_groupPeaks_mzClust.md)
  for the core API function.

For specific examples and description of the method and settings see the
help pages of the individual parameter classes listed above.

## Usage

``` r
groupChromPeaks(object, param, ...)

# S4 method for class 'XcmsExperiment,Param'
groupChromPeaks(object, param, msLevel = 1L, add = FALSE)

PeakDensityParam(
  sampleGroups = numeric(),
  bw = 30,
  minFraction = 0.5,
  minSamples = 1,
  binSize = 0.25,
  ppm = 0,
  maxFeatures = 50
)

MzClustParam(
  sampleGroups = numeric(),
  ppm = 20,
  absMz = 0,
  minFraction = 0.5,
  minSamples = 1
)

NearestPeaksParam(
  sampleGroups = numeric(),
  mzVsRtBalance = 10,
  absMz = 0.2,
  absRt = 15,
  kNN = 10
)

# S4 method for class 'PeakDensityParam'
as.list(x, ...)

# S4 method for class 'XCMSnExp,PeakDensityParam'
groupChromPeaks(object, param, msLevel = 1L, add = FALSE)

# S4 method for class 'XCMSnExp,MzClustParam'
groupChromPeaks(object, param, msLevel = 1L)

# S4 method for class 'XCMSnExp,NearestPeaksParam'
groupChromPeaks(object, param, msLevel = 1L, add = FALSE)
```

## Arguments

- object:

  The data object on which the correspondence analysis should be
  performed. Can be an
  [`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md),
  [`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
  object.

- param:

  The parameter object selecting and configuring the algorithm.

- ...:

  Optional parameters.

- msLevel:

  `integer(1)` defining the MS level on which the chromatographic peak
  detection should be performed.

- add:

  `logical(1)` (if `object` contains already chromatographic peaks, i.e.
  is either an `XCMSnExp` or `XcmsExperiment`) whether chromatographic
  peak detection results should be **added** to existing results. By
  default (`add = FALSE`) any additional `findChromPeaks` call on a
  result object will remove previous results.

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

- maxFeatures:

  For `PeakDensityParam`: `numeric(1)` with the maximum number of peak
  groups to be identified in a single mz slice.

- absMz:

  For `NearestPeaksParam` and `MzClustParam`: `numeric(1)` maximum
  tolerated distance for m/z values.

- mzVsRtBalance:

  For `NearestPeaksParam`: `numeric(1)` representing the factor by which
  m/z values are multiplied before calculating the (euclician) distance
  between two peaks.

- absRt:

  For `NearestPeaksParam`: `numeric(1)` maximum tolerated distance for
  retention times.

- kNN:

  For `NearestPeaksParam`: `integer(1)` representing the number of
  nearest neighbors to check.

- x:

  The parameter object.

## Value

For `groupChromPeaks`: either an
[`XcmsExperiment()`](https://sneumann.github.io/xcms/reference/XcmsExperiment.md)
or
[`XCMSnExp()`](https://sneumann.github.io/xcms/reference/XCMSnExp-class.md)
object with the correspondence result.

## References

Smith, C.A., Want E.J., O'Maille G., Abagyan R., and Siuzdak G. (2006)
"XCMS: Processing Mass Spectrometry Data for Metabolite Profiling Using
Nonlinear Peak Alignment, Matching, and Identification" *Anal. Chem.*
78:779-787. doi: [10.1021/ac051437y](https://doi.org/10.1021/ac051437y)

Katajamaa, M., Miettinen, J., Oresic, M. (2006) "MZmine: Toolbox for
processing and visualization of mass spectrometry based molecular
profile data". *Bioinformatics*, 22:634-636. doi:
[10.1093/bioinformatics/btk039](https://doi.org/10.1093/bioinformatics/btk039)

Kazmi S. A., Ghosh, S., Shin, D., Hill, D.W., and Grant, D.F. (2006)
"Alignment of high resolution mass spectra: development of a heuristic
approach for metabolomics. *Metabolomics* Vol. 2, No. 2, 75-83.

## Author

Colin Smith, Johannes Rainer
